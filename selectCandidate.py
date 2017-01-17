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
    parser = argparse.ArgumentParser(description='Selecting long genomic sequence candidates from consensus FASTA generated from samtools pileup')
    parser.add_argument('-c', help='The FASTA consensus file from the pileup',required=True, dest="pileupFASTA")
    parser.add_argument('-l', help='Length of consensus fragments to retain, shorter ones are discarded. Default is 3000 bps', required=False, dest="shortestContig", default="3000")
    return parser.parse_args()

def getRecords(pileupFASTA, shortestContig):
    return (rec for rec in SeqIO.parse(pileupFASTA, "fasta") if len(rec) >= shortestContig )

def makePairs(candidates):
    pairs = []
    if len(list(candidates)) == 1:
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
    fastaRecords = getRecords(args.pileupFASTA, args.shortestContig)
    # in some cases there is only a single candidate (i.e. homozygous or hemizygous cases, like DRB1,4,7 etc)
    # we have to be sure we are emitting a pair, not a single candidate (the same allele in those cases)
    pairs = makePairs(fastaRecords)
    for record in pairs:
        # have to change ":"s in filenames
        seqID = record.id.replace(":","_")
        fileName = seqID + "_" + os.path.basename(args.pileupFASTA) + ".candidate.fasta"
        SeqIO.write(record,fileName, "fasta")

if __name__ == '__main__':
    main()
