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
    prevId = ''
    fh = open(pileup,"r")
    for line in fh:
        (seqId,seqChr) = getConsensusChar(line, depth)
        # if it is the very first sequence, initialize the "previous" seqId as the first (we want to see seqID changes)
        if len(seqStr) == 0:
            prevId = seqId
        # add the next base to the growing sequence string
        seqStr += seqChr 
        # add a newline if too long (aka UNIX fold)
        if len(seqStr) % 80 == 0:
            seqStr += '\n'
        # start a new sequence if ID has changed
        if prevId != seqId:
            printOutFASTA(prevId,seqStr)
            prevId = seqId
            seqStr = ''
    printOutFASTA(seqId,seqStr)

def printOutFASTA(seqId, seqStr):
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
    seqId = ''
    baseDepth = 0
    if len(sl) == 4:
        (seqId, baseDepth) = (sl[0], int(sl[3]) )
    else:
        try:
            (seqId, baseDepth, pl ) = (sl[0], int(sl[3]), sl[4])
        except:
            print "offending line:"
            print aLine
            raise 
    commonBase = ''
    if minDepth <= baseDepth and baseDepth > 0:
        commonBase = getMostCommonBase(pl)
    else:
        commonBase = 'N'
    return seqId, commonBase

if __name__ == "__main__":
    makeConsensus()
