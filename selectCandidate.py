# Selecting the best HLA consensus sequences from a samtools pileup consensus FASTA
import click
import copy
from Bio import SeqIO


@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--fasta','-f', type=str, help='multi FASTA input file to select candidates from')
@click.option('--length','-l', type=int, help='length of the shortest contig')
def getFastaRecords(fasta,length):
    """Selecting only particular entries from a multi-FASTA (right now only considering length)"""
    frs = list()
    for rec in SeqIO.parse(fasta, "fasta"):
        if len(rec) >= length:
            frs.append(rec)
    printSeqs(makePairs(frs))

def printSeqs(records):
    for seq in records:
        print seq.format("fasta")

def makePairs(candidates):
    # we are expecting a list here, not a generator
    pairs = []
    if len(candidates) == 1:
        clone = copy.deepcopy(candidates[0])
        clone.id = clone.id + "clone"
        clone.name = clone.name + "clone"
        clone.description = clone.description + "clone"
        candidates.append(clone)
        
    return candidates
 

if __name__ == '__main__':
    getFastaRecords()
