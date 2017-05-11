#!/usr/bin/env python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import click
import sys


def validate_locus(ctx,param,value):
    try:
        assert value in ["HLA-A","HLA-B", "HLA-C", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DPB2", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-DRB2", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DRB6", "HLA-DRB7", "HLA-DRB8", "HLA-DRB9", "HLA-E", "HLA-F", "HLA-G", "HLA-H", "HLA-HFE", "HLA-J", "HLA-K", "HLA-L", "HLA-P", "HLA-V", "HLA-Y"]
    except:
        raise click.BadParameter('Please define locus as HLA-A, HLA-B, HLA-DRB1 ... as you can find in awk -F[\*\ ] \'/^DE/ && /HLA/ {print $4}\' hla.dat|sort -u')
    return value

@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--dat','-d', type=str, help='the IMGT/HLA reference hla.dat file from ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla.dat')
@click.option('--locus','-l', type=str, help='the locus [either HLA-A, HLA-B, HLA-DRB1 ...]',default="HLA-A",callback=validate_locus)
def printExons(dat,locus):
    fixedFile = fixIMGTfile(dat)
    printExonsOnly(fixedFile,locus)

def printExonsOnly(fh,locus):
    """ Go through each entry, and print out exons 1 and 2"""

    for seq_record in SeqIO.parse(fh,"imgt"):
        if seq_record.description.startswith(locus) and len(seq_record.seq) > 1:
            # the new exon record we are extracting
            newSeq = ""
            for f in seq_record.features:
                if f.type == "exon":
                    if f.qualifiers['number'] == ['2'] :
                        newSeq += seq_record.seq[f.location.start:f.location.end ]
                    if locus in ["HLA-A","HLA-B", "HLA-C"]  and f.qualifiers['number'] == ['3']:
                        newSeq += seq_record.seq[f.location.start:f.location.end ]
            sys.stdout.write( SeqRecord(newSeq, id=seq_record.id, description=seq_record.description).format("fasta") )

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
    printExons()


#a=SeqIO.read("test.dat","embl")
##print a
#print a.id
#print a.description
#newSeq = ""
#for f in a.features:
#    if f.type == "exon" and (f.qualifiers['number'] == ['2'] or f.qualifiers['number'] == ['3'] ):
#        newSeq += a.seq[f.location.start:f.location.end ]
#print newSeq
#record = SeqRecord(newSeq,id=a.id,description=a.description)
#print record
#print record.format("fasta")
