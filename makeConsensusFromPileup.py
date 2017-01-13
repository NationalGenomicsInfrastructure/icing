import operator
import click


@click.command()
@click.option('--pileup','-p', type=str, help='The pileup file to use for consensus generation')
@click.option('--depth', '-d', type=int, help='Minimal expected depth. With coverage lower than this N will be inserted (default = 10).', default=10)
def makeConsensus(pileup,depth):
    """
    Utility to make a consensus from samtools pileup. Use samtools as
    samtools mpileup -B -d 1000 -Q 10 -A  file.bam > file.pileup
    and then feed the pileup to this script
    makeConsensusFromPileup.py -p file.pileup
    It is a naiive consensus generator, indels are ignored.
    """
    seqStr = ''
    seqId = ''
    fh = open(pileup,"r")
    for line in fh:
        (seqId,seqChr) = getConsensusChar(line, depth)
        seqStr += seqChr 
        if len(seqStr) % 80 == 0:
            seqStr += '\n'

    print ">" + seqId
    print seqStr
    
def getMostCommonBase(aPileup):
    pl = aPileup.upper()
    cA = pl.count('A')
    cC = pl.count('C')
    cG = pl.count('G')
    cT = pl.count('T')
    letters = ['A','C','G','T']
    counts = [ cA, cC, cG, cT]
    base = ''
    if sum(counts) > 0:
        index, value = max(enumerate(counts), key=operator.itemgetter(1))
        base = letters[index]
    return base

def getConsensusChar(aLine, minDepth):
    sl = aLine.split()
    (seqId, baseDepth, pl ) = (sl[0], sl[3],sl[4])
    commonBase = ''
    if minDepth <= baseDepth:
        commonBase = getMostCommonBase(pl)
    else:
        commonBase = 'N'
    return seqId, commonBase

if __name__ == "__main__":
    makeConsensus()
