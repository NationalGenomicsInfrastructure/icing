from __future__ import print_function
import operator
import click

@click.command()
@click.option('--pileup','-p', type=str, help='The pileup file to use for consensus generation')
@click.option('--depth', '-d', type=int, help='Minimal expected depth. With coverage lower than this N will be inserted (default = 10).', default=10)
@click.option('--nonN', '-n', type=float, help='Minimal expected percentage of non-N bases. (default = 90.0).', default=90.0)
def makeConsensus(pileup, depth, nonn):
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
            printOutFASTA(prevId,seqStr, nonn)
            prevId = seqId
            seqStr = ''
    printOutFASTA(seqId,seqStr, nonn)

def printOutFASTA(seqId, seqStr, nonn):
	"""
	We are writing out sequences that are containing at least nonN% of non-N parts. 
	I. e. if the sequence is ACCTGANNNN and we have a 90% limit, it will not be printed out.
	"""
        percentN = 1.0 - float(len(seqStr) - seqStr.count("N"))/float(len(seqStr))
	if(percentN <= 1.0 - nonn/100.0):
            print(">" + seqId + " Ns: " + str(percentN))
            print(seqStr)#,end='')
    
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
            print("offending line:")
            print(aLine)
            raise 
    commonBase = ''
    if minDepth <= baseDepth and baseDepth > 0:
        commonBase = getMostCommonBase(pl)
    else:
        commonBase = 'N'
    return seqId, commonBase

if __name__ == "__main__":
    makeConsensus()
