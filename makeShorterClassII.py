#!/usr/bin/env python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import click
import sys
import os


def validate_locus(ctx,param,value):
    try:
        assert value in ["HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DPB2", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-DRB2", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DRB6", "HLA-DRB7", "HLA-DRB8", "HLA-DRB9", "HLA-E", "HLA-F", "HLA-G", "HLA-H", "HLA-HFE", "HLA-J", "HLA-K", "HLA-L", "HLA-P", "HLA-V", "HLA-Y"]
    except:
        raise click.BadParameter('Please define locus as HLA-A, HLA-B, HLA-DRB1 ... as you can find in awk -F[\*\ ] \'/^DE/ && /HLA/ {print $4}\' hla.dat|sort -u')
    return value

@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--dat','-d', type=str, help='the IMGT/HLA reference hla.dat file from ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla.dat')
@click.option('--locus','-l', type=str, help='the locus [HLA-DRB1 ...]',default="HLA-DRB1",callback=validate_locus)
def printShortClassII(dat,locus):
    fixedFile = fixIMGTfile(dat)
    printShortenedFASTA(fixedFile,locus)

def fixIMGTfile(hladat):
    """
    For some reason IMGT is not following the standard EMBL format, so we have to add entries to the ID line
    So, we will add an extra "IMGT;" as:
    ID   HLA00001; SV 1; standard; DNA; HUM; 3503 BP.
     going to be -> 
    ID   HLA00001; SV 1; standard; DNA; IMGT; HUM; 3503 BP.
    """
    newFileName = "fixed" + os.path.basename(hladat)
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

def printShortenedFASTA(fixed,locus):
    """
    Prints out a only genomic FASTA sequences with ID, but only  -200 bases before exon 2
    """
    for seq_record in SeqIO.parse(fixed,"imgt"):
        if seq_record.description.startswith(locus) and len(seq_record.seq) > 1:
            # the new exon record we are extracting
            newSeq = ""
            extracted = False   # we will not be able to extract all parts as many sequences has missing introns
            for f in seq_record.features:   # we hope the features are in order
                if f.type == "intron" and f.qualifiers['number'] == ['1']:
                    newSeq = seq_record.seq[f.location.end-200:]
                    extracted = True
                    break   # we have the whole sequence now
            if extracted:
                sys.stdout.write( SeqRecord(newSeq, id=seq_record.id, description=seq_record.description).format("fasta") )

if __name__ == "__main__":
    printShortClassII()


