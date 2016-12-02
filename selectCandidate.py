# Selecting the best HLA consensus sequences from a samtools pileup consensus FASTQ
from __future__ import print_function
import argparse
import os.path
from Bio import SeqIO
import numpy
import matplotlib.pyplot as plt
import scipy.stats.stats as stats

def parse_args():
    parser = argparse.ArgumentParser(description='Selecting genomic sequence candidates from consensus FASTQ generated from samtools pileup')
    parser.add_argument('-c', help='The FASTQ consensus file from the pileup',required=True, dest="pileupFASTQ")
    parser.add_argument('-l', help='Length of consensus fragments to retain, shorter ones are discarded. Default is 3000 bps', required=False, dest="shortestContig", default="3000")
    parser.add_argument('-a', help='Averge PHRED qualities in the consensus FASTQ. Entries with lower average are discarded. Default is 60', required=False, dest="averageLimit", default="60")
    return parser.parse_args()

def selectLongCandidates(fastqRecords, shortestContig):
    return (rec for rec in fastqRecords if len(rec) >= shortestContig )

def selectHQCandidates(fastqRecords, averageLimit):
    hqR = []
    for record in fastqRecords:
        phred = record.letter_annotations["phred_quality"]
#        print record.id + " ###########################################################################"
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
        #print "skewness: " + str(stats.skew(phred))
        if stats.skew(phred) < 0.0:
            skR.append(record)
#        print "kurtosis: " + str(stats.kurtosis(phred))
#        num_bins=15
#        n, bins, patches = plt.hist(phred, num_bins, normed=1, facecolor='green', alpha=0.5)
#        plt.title("FASTQ histogram " + record.id)
#        plt.show(True)
    return skR

def getRecords(pileupFASTQ):
    return (rec for rec in SeqIO.parse(pileupFASTQ, "fastq") if min(rec.letter_annotations["phred_quality"]) >= 30)

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
    for record in skewedCandidates:
        # print out without the extra newline
        SeqIO.write(record,record.id + "_" +os.path.basename(args.pileupFASTQ) + ".candidate.fasta", "fasta")

if __name__ == '__main__':
    main()
