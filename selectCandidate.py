# Selecting the best HLA consensus sequences from a samtools pileup consensus FASTQ
from __future__ import print_function
import argparse
import os.path
from Bio import SeqIO
import numpy
import matplotlib.pyplot as plt
import scipy.stats.stats as stats
import copy

def parse_args():
    parser = argparse.ArgumentParser(description='Selecting genomic sequence candidates from consensus FASTQ generated from samtools pileup')
    parser.add_argument('-c', help='The FASTQ consensus file from the pileup',required=True, dest="pileupFASTQ")
    parser.add_argument('-l', help='Length of consensus fragments to retain, shorter ones are discarded. Default is 3000 bps', required=False, dest="shortestContig", default="3000")
    parser.add_argument('-a', help='Averge PHRED qualities in the consensus FASTQ. Entries with lower average are discarded. Default is 60', required=False, dest="averageLimit", default="60")
    return parser.parse_args()

def selectLongCandidates(fastqRecords, shortestContig):
    lc = []
    for rec in fastqRecords:
        #print("%s %d" % (rec.id, len(rec)))
        if len(rec) >= shortestContig:
            lc.append(rec)
    return lc
    #return (rec for rec in fastqRecords if len(rec) >= shortestContig )

def selectHQCandidates(fastqRecords, averageLimit):
    hqR = []
    for record in fastqRecords:
        #print(record.id + " ###########################################################################")
        phred = record.letter_annotations["phred_quality"]
#        print "average: " + str(numpy.average(phred))
#        print "mean: " + str(numpy.mean(phred))
#        print "median: " + str(numpy.median(phred))
#        print "std: " + str(numpy.std(phred))
#        print "var: " + str(numpy.var(phred))
        if numpy.average(phred) > averageLimit:
            hqR.append(record)
    return hqR


def selectSkewedCandidates(fastqRecords):
    skR = []
    for record in fastqRecords:
        phred = record.letter_annotations["phred_quality"]
        #print ("skewness: " + str(stats.skew(phred))
        if stats.skew(phred) < 0.0:
            skR.append(record)
#        print "kurtosis: " + str(stats.kurtosis(phred))
#        num_bins=15
#        n, bins, patches = plt.hist(phred, num_bins, normed=1, facecolor='green', alpha=0.5)
#        plt.title("FASTQ histogram " + record.id)
#        plt.show(True)
    return skR

def getRecords(pileupFASTQ):
    fqRec = []
    for rec in SeqIO.parse(pileupFASTQ, "fastq"):
        #print("%s %f" % (rec.id, numpy.median(rec.letter_annotations["phred_quality"]) ))
        try:
            if numpy.median(rec.letter_annotations["phred_quality"]) >= 25:
                fqRec.append(rec)
        except AttributeError:
           pass 
    return fqRec

    #return (rec for rec in SeqIO.parse(pileupFASTQ, "fastq") if min(rec.letter_annotations["phred_quality"]) >= 30)

def makePairs(candidates):
    pairs = []
    if len(candidates) == 1:
        clone = copy.deepcopy(candidates[0])
        clone.id = clone.id + "clone"
        clone.name = clone.name + "clone"
        clone.description = clone.description + "clone"
        candidates.append(clone)
        for c in candidates:
            print(c)
        
    return candidates
    

def main():
    # parse command line
    args = parse_args()
    fastqRecords = getRecords(args.pileupFASTQ)
    # select long candidates: get the top 10% of length distribution (regarding length)
    longCandidates = selectLongCandidates(fastqRecords, int(args.shortestContig))
    # select candidates with high quality bases (based on median base quality value and stddev)
    hqCandidates = selectHQCandidates(longCandidates, int(args.averageLimit))
    # select sequences where the base quality distribution is right skewed (are skewness is negative)
    skewedCandidates = selectSkewedCandidates(hqCandidates)
    # in some cases there is only a single candidate (i.e. homozygous or hemizygous cases, like DRB1,4,7 etc)
    # we have to be sure we are emitting a pair, not a single candidate (the same allele in those cases)
    pairs = makePairs(skewedCandidates)
    for record in pairs:
        # have to change ":"s in filenames
        seqID = record.id.replace(":","_")
        fileName = seqID + "_" + os.path.basename(args.pileupFASTQ) + ".candidate.fasta"
        SeqIO.write(record,fileName, "fasta")

if __name__ == '__main__':
    main()
