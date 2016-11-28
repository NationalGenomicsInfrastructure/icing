# Selecting the best HLA consensus sequences from a samtools pileup consensus FASTQ
import argparse
from Bio import SeqIO
import numpy
import matplotlib.pyplot as plt
import scipy.stats.stats as stats

def parse_args():
    parser = argparse.ArgumentParser(description='Selecting genomic sequence candidates from consensus FASTQ generated from samtools pileup')
    parser.add_argument('-c', help='The FASTQ consensus file from the pileup',required=True, dest="pileupFASTQ")
    return parser.parse_args()

def selectLongCandidates(aSeq):
    return aSeq

def selectHQCandidates(aSeq):
    return aSeq

def selectSkewedCandidates(aSeq):
    return aSeq

def reduceCloseRelatives(aSeq):
    return aSeq

def getRecords(pileupFASTQ):
    return (rec for rec in SeqIO.parse(pileupFASTQ, "fastq") if min(rec.letter_annotations["phred_quality"]) >= 30)

def main():
    # parse command line
    args = parse_args()
    fastqRecords = getRecords(args.pileupFASTQ)
    # select long candidates: get the top 10% of length distribution (regarding length)
    longCandidates = selectLongCandidates(fastqRecords)
    # select candidates with high quality bases (based on median base quality value and stddev)
    hqCandidates = selectHQCandidates(longCandidates)
    # select sequences where the base quality distribution is right skewed (are skewness is negative)
    skewedCandidates = selectSkewedCandidates(hqCandidates)
    # reduce candidates if there are still too many similar ones
    bestOnes = reduceCloseRelatives(skewedCandidates)
    for ca in bestOnes:
        print ca

if __name__ == '__main__':
    main()
