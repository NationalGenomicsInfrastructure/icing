#!/usr/bin/env python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import click
import sys
import ssw
import operator

# it is a global variable - once it is refactored to be OO, it will disappear
# the whole point is to map HLA00005 to HLA*02:01:01:01
gHLAtypes = {}

def validate_locus(ctx,param,value):
    try:
        assert value in ["HLA-A","HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DRB6", "HLA-DRB7", "HLA-DRB8", "HLA-DRB9"]
    except:
        raise click.BadParameter('Please define locus as HLA-A, HLA-B, HLA-DRB1 ... as you can find in awk -F[\*\ ] \'/^DE/ && /HLA/ {print $4}\' hla.dat|sort -u')
    return value

@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--cons','-c', type=str, help='FASTA file containing consensuses (or a single consensus)')
@click.option('--dat','-d', type=str, help='the IMGT/HLA reference hla.dat file from ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla.dat')
@click.option('--locus','-l', type=str, help='the locus [either HLA-A, HLA-B, HLA-DRB1 ...]',default="HLA-A",callback=validate_locus)
def doGenotype(cons, dat, locus):
    fixedFile = fixIMGTfile(dat)
    ( primaryExons, secondaryExons, intronsAndUTRs) = getCompartmenstForAlleles(fixedFile,locus)
	# for each consensus
	#	preselect types considering only the important exons
	#	refine the preselected list by checking mismatch in the secondary exons
	#	final touches by looking at the introns/UTRs
    for seq_record in SeqIO.parse(cons,"fasta"):
		print "Processing consensus: " + seq_record.id
		# select genotypes considering only the important exons
		genotypes = preSelectTypes(primaryExons,seq_record)
		print genotypes
		if len(genotypes) == 0:		# we have failed for some reason
			print "Could not find a proper type for consensus"
		elif len(genotypes) == 1:	# there is a single genotype only: no need to shrink the candidate set
			print "Final HLA type for consensus: " + genotypes[0]
		else:	# there are more than one type candidates, go for exons
			genotypes = selectGenotypesConsideringCommonExons(genotypes, secondaryExons,seq_record)
			if len(genotypes) > 1:
				genotypes = selectGenotypesConsideringIntronsAndUTRs(genotypes, intronsAndUTRs,seq_record)
			print "Final HLA type for consensus:" 
			for gt in genotypes:
				print gHLAtypes[gt]


def getCompartmenstForAlleles(fixedIMGT,locus):
	"""
	Go through each entry in the IMGT/HLA EMBL file, and store out exons in a dictionary.
	For Class-I we are storing exons 2 & 3, for Class-II only exon 2. 
	The data is is like

	primary{}
	ex2='ACTGATCGATCGATACG'
	ex3='CCAGGCCTGGATCGCATTAGC'
	primary['HLA000101']=[ex2,ex3]
	{'HLA000101': ['ACTGATCGATCGATACG', 'CCAGGCCTGGATCGCATTAGC']}
	"""
	print "Processing reference IMGT file"
	primary = {}
	secondary = {}
	intronsAndUTRs = {}
	for seq_record in SeqIO.parse(fixedIMGT,"imgt"):
		gHLAtypes[seq_record.id] = seq_record.description
		# if it is the correct locus and there is a sequence record (not a deleted one)
		if seq_record.description.startswith(locus) and len(seq_record.seq) > 1:
			primary[seq_record.id] = getPrimaryExons(seq_record, locus)
			secondary[seq_record.id] = getSecondaryExons(seq_record, locus)
			intronsAndUTRs[seq_record.id] = getIntronsAndUTRs(seq_record, locus)
		
	print "ready"	
	return (primary,secondary,intronsAndUTRs)

def getPrimaryExons(sr,locus):
	"""
	Primary exons are exon 2 and 3 for HLA-A,B,C and exon 2 for all the rest.
	TODO: Note, this is actually wrong. There are quite a few other loci where the important polymorphic
	exons is not exon 2 only, but in the moment we are ignoring this fact.
	"""
	exonList = []
	for f in sr.features:
		if f.type == "exon":
			if f.qualifiers['number'] == ['2'] :
				exonList.append( sr.seq[f.location.start:f.location.end] )
			if locus in ["HLA-A","HLA-B", "HLA-C"]  and f.qualifiers['number'] == ['3']:
				exonList.append( sr.seq[f.location.start:f.location.end] )
	return exonList

def getSecondaryExons(sr,locus):
	"""
	We are adding all the other exons as secondary.
	Instead of a list, we have to add qualifiers as well, since for some alleles there are 5 exons defined, and only 3 for the other, etc.
	The data looks like:
	{ 
		'allele1': { 'ex1': seq, 'ex4': seq, 'ex5': seq},
		'allele2': { 'ex1': seq, 'ex4': seq, 'ex5': seq},
		'allele3': { 'ex1': seq, 'ex4': seq},
		...
	}
	For a set like that only ex1 and ex4 are in the smallest common set. 
	So, we are returning with the { 'ex1': seq, 'ex4': seq, 'ex5': seq,...} part
	"""
	exonDict = {}
	for f in sr.features:
		if f.type == "exon":	# consider only exons
			# treat Class-I and Class-II separately: it can be faster, but it is more readable this way
			# Class-I first
			if locus in ["HLA-A","HLA-B", "HLA-C"] and f.qualifiers['number'][0] not in ['2','3']:
				exonDict['ex'+f.qualifiers['number'][0] ] = str( sr.seq[f.location.start:f.location.end] )
			# Class-II and other (all non-Class-I entries)
			elif locus not in ["HLA-A","HLA-B", "HLA-C"] and f.qualifiers['number'] != ['2']:
				exonDict['ex'+f.qualifiers['number'][0] ] = str( sr.seq[f.location.start:f.location.end] )
		
	return exonDict

def getIntronsAndUTRs(sr,locus):
	"""
	All non-exonic compartments are added
	"""
	nonExonsList = []
	for f in sr.features:
		if f.type != "exon":
			nonExonsList.append( sr.seq[f.location.start:f.location.end] )	
	return nonExonsList


def getBestScoringAlleles( sortedTupleList ):
	"""
	We are expecting a sorted list of tuples, containing allele names and alignment scores like
	[('HLA10395.1', 545), ('HLA06895.1', 545), ('HLA14515.1', 540), ('HLA12165.1', 532), ('HLA03670.1', 532)]
	We are returning with a list of best matches:
	['HLA10395.1','HLA06895.1']
	"""	
	(firstAllele,bestScore) = sortedTupleList[0]
	bestAlleles = [firstAllele]
	for (allele,score) in sortedTupleList[1:]:
		if bestScore > score:
			break
		else:
			bestAlleles.append(allele)
	return bestAlleles	 


def preSelectTypes(primary,consensus):
	"""
	For each primary exon (or exon pair)
		make an alignment for exon 2, and put the result into a dictionary as alignmentEx2['allele'] = #score 
		# we are getting [score, matches, mismatches, inserts, deletions] if we want to 
		if there is exon 3 also, 
			do an alignmentEx3['allele'] = #mismatches
			merge the two in a way that sort both, keep the best for both, and make an intersect
		else
			sort, and keep the best alignments only
		
	"""
	print "Selecting best alleles using primary exons"
	# we are going to use https://github.com/vishnubob/ssw that is using
	# https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library for SW and SAM output
	# unfortunatelly it can not made multiprocess yet, but it is fast enough
	alignmentEx2 = {}
	alignmentEx3 = {}
	sw = ssw.Aligner()
	for allele,exons in primary.items():
		alignment = sw.align(str(consensus.seq),exons[0])
		alignmentEx2[allele] = alignment.score
		# ditto for exon 3 (TODO: no class-I/class-II checking here)
		alignment = sw.align(str(consensus.seq),exons[1])
		alignmentEx3[allele] = alignment.score
		# print allele + " scores: exon 2: " + str(alignmentEx2[allele]) + " exon 3 " + str(alignmentEx3[allele])
	# sort the dict by values and reverse the result
	bestEx2 = getBestScoringAlleles( sorted(alignmentEx2.items(),key=operator.itemgetter(1),reverse=True) )
	bestEx3 = getBestScoringAlleles( sorted(alignmentEx3.items(),key=operator.itemgetter(1),reverse=True) )
	# return with the intersect of the two sets, leaving only entries that are in the besty matching set for both exon 2 and exon 3
	print "done"
	return list(set(bestEx2) & set(bestEx3))

def getCommonExons(genotypes,exons):
	"""
	genotypes is a list of strings like
	['HLA10395.1','HLA06895.1']
	exons is a dict having dicts as values like
	{	 
        'HLA10395.1': { 'ex1': seq, 'ex4': seq, 'ex5': seq},
        'HLA10495.1': { 'ex1': seq, 'ex4': seq, 'ex5': seq},
        'HLA10355.1': { 'ex1': seq, 'ex4': seq},
        ...
    }
	"""
	# as the initial value of the common set, get exons of the first allele
	print "Searching for common exons"
	commonExons = set(exons[genotypes[0]].keys())
	print "starting from "
	print genotypes[0], commonExons
	for gt in genotypes[0:]:
		commonExons = commonExons & set( exons[gt].keys() )
		#print gt,exons[gt].keys(),commonExons
	print "Common exons: "
	print commonExons
	return commonExons

def selectGenotypesConsideringCommonExons(genotypes,secondary,consensus):
	# select a set of common exons
	numOfCandidates = len(genotypes)
	commonExons = getCommonExons(genotypes,secondary)
	# select the set of best matching alleles for the very first exon in the list
	# we have to do it for all exons, but initialize for the first
	newGenotypes = set(genotypes) & getBestScoringAllelesForExon(genotypes,commonExons.pop(),secondary,consensus)
	while commonExons:
		newGenotypes = set(newGenotypes) & getBestScoringAllelesForExon(newGenotypes,commonExons.pop(),secondary,consensus)
	print "new genotypes: "
	print newGenotypes
	if numOfCandidates == len(newGenotypes):
		# we were not able to decrease the number of candidates
		return list(newGenotypes)
	else:
		return selectGenotypesConsideringCommonExons(list(newGenotypes),secondary,consensus)

def getBestScoringAllelesForExon(genotypes,commonExon,secondary,consensus):
	exAlign = {}
	sw = ssw.Aligner()
	for allele in genotypes:
		print allele
		exonSeq = secondary[allele][commonExon]
		alignment = sw.align(str(consensus.seq),exonSeq)
		exAlign[allele] = alignment.score
	return set(getBestScoringAlleles( sorted(exAlign.items(),key=operator.itemgetter(1),reverse=True)) )

def selectGenotypesConsideringIntronsAndUTRs(genotypes,introns,consensus):
	return genotypes
	
def fixIMGTfile(hladat):
    """
    For some reason IMGT is not following the standard EMBL format, so we have to add entries to the ID line
    So, we will add an extra "IMGT;" as:
    ID   HLA00001; SV 1; standard; DNA; HUM; 3503 BP.
     going to be -> 
    ID   HLA00001; SV 1; standard; DNA; IMGT; HUM; 3503 BP.
    """
    newFileName = "fixed" + hladat
    fixedfh = open(newFileName,"w")
         
    with open(hladat,"r") as fh:
        for line in fh:
            if line.startswith("ID"):
                parts = line.split()
                line = parts[0] + "   " + \
                        parts[1] + " " +\
                        parts[2] + " " +\
                        parts[3] + " " +\
                        parts[4] + " " +\
                        parts[5] + " " +\
                        "IMGT; " +\
                        parts[6] + " " +\
                        parts[7] + " " +\
                        parts[8] + "\n"
            fixedfh.write(line)
    return newFileName


if __name__ == "__main__":
	doGenotype()

#		import pdb
#		pdb.set_trace()	
