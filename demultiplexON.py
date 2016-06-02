import swalign
import argparse
import yaml
from Bio.Seq import Seq
from Bio import SeqIO

"""
To separate (demultiplex) pooled HLA reads, we have to generate handle sequences (and their reverse complement parts)
to align to the FASTQ reads. For Oxford Nanopore we are expecting ~90% reliability, meaning that in the ~50bps long 
handle we will have approx 5bps errors, indels mostly.
Processing steps:
- generate handles from yaml
- get reads from the FASTQ, and try to align a handle to the beginning of the read. Note, there are usually a ~15bps
  chunk extra prefix that precedes the prefix
- if a read has only a single best scoring handle, store as it is
- if the handle alignment is ambiguous, store in a separate file (we will find out later what to do, maybe use 
  different alignment scoring)
- if no handle found, store as junk
"""

def parse_args():
    parser = argparse.ArgumentParser(description='Demultiplexing Oxford Nanopore reads')
    parser.add_argument('-y', help='the YAML definiton of handle prefix,postfix and indexes',required=True, dest="indexes")
    parser.add_argument('-r', help='FASTQ file containing reads',required=True, dest="reads")
    parser.add_argument('-m', help='number of mismatches expected in the handle sequence',required=True, dest="mm")

    return parser.parse_args()


def readYAMLConf(yamlfile):
    config = {}
    with open(yamlfile, 'r') as ymlfile:
        cfg = yaml.load(ymlfile)

        for section in cfg:
            config[section] = cfg[section]
    return config

def generateSeqHandles(anIndexCfg):
    """
        The YAML config file to parse is like:
    
    handles:
        prefix: "TTAGTCTCCGACGGCAGGCTTCAAT"
        postfix: "ACGCACCCACCGGGACTCAG"
    indexes: [
        "ACAGTC",
        "TGATGC",
        "TCTCAG"
    ]

    There is a handle at one end of each sequence which is as follows:
    TTAGTCTCCGACGGCAGGCTTCAAT-ACAGTC-ACGCACCCACCGGGACTCAG
              prefix         -index -      postfix
    """
    forwardIdx= []     # the result array to collect handle sequence strings
    handlePrefix = anIndexCfg["handles"]["prefix"]
    handlePostfix = anIndexCfg["handles"]["postfix"]
    for index in anIndexCfg["indexes"]:
        forwardIdx.append(handlePrefix + index + handlePostfix)
    
    reverseIdx = []       # to collect reverse complements
    for handle in forwardIdx:
        seq = Seq(handle)
        rc = str(seq.reverse_complement())
        reverseIdx.append(rc)

    return (forwardIdx,reverseIdx)

def getBestMatches(aSeq, aHandleList, anSWAligner, aMismatch):
    bestMatches = []
    for handle in aHandleList:
        al = anSWAligner.align(aSeq,handle)
        if al.mismatches < aMismatch:
            bestMatches.append(handle)
    return bestMatches

def swFactory():
    match = 2
    mismatch = -1
    gap_penalty = -1
    gap_extension_penalty = -1
    gap_extension_decay = 0.0
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    return swalign.LocalAlignment( (scoring), 
                                    gap_penalty, 
                                    gap_extension_penalty, 
                                    gap_extension_decay=gap_extension_decay, 
                                    verbose=False, 
                                    globalalign=False, 
                                    full_query=False)

def fileFactory(aHandleList):
    filesToSave = {}
    for seq in aHandleList:
        fh = open(seq+".dmx.fastq","w")
        filesToSave[seq] = fh
    return filesToSave

def writeFASTQRecord(aRecord,aFASTQFile):
    readId = aRecord.id
    seqStr = aRecord.__dict__['_seq'].__dict__['_data']
    quals = aRecord.__dict__['_per_letter_annotations']['phred_quality']
    qualsStr = ""
    for c in quals:
        qualsStr += chr(c+33)

    #import pdb;pdb.set_trace()
    aFASTQFile.write("@" + readId+" "+ str(len(seqStr))+ " "+ str(len(qualsStr)) +"\n")
    aFASTQFile.write(seqStr+"\n")
    aFASTQFile.write("+\n")
    aFASTQFile.write(qualsStr+"\n")   
    aFASTQFile.flush()


def demultiplexFastqByBestMatch(aFASTQFile,aHandleList,aMismatch,isForward=True):
    # we are exporting each handle into a different file
    # this dictionary has the sequence handles as keys and the files as values
    exportFiles = fileFactory(aHandleList)
    sw = swFactory()
    fh = open(aFASTQFile,'rU')
    for idx, record in enumerate(SeqIO.parse(fh, "fastq")):
        seqStr = str(record.seq)
        # if we are looking for reverse sequence, get it from the end
        if isForward:
            seqStr = seqStr[0:100]
        else:
            #import pdb; pdb.set_trace()
            seqStr = seqStr[len(seqStr)-99:]
        # bestMatches is a list of handles having the same alignment score
        # if there is only one, save it, else ignore ambiguities
        bestMatches = getBestMatches(seqStr, aHandleList, sw, aMismatch)   # get the best matches for all the handles
        if len(bestMatches) == 1:                                       # there is a single best match: store it
            # unfortunately FASTQ export for very long reads looks to be buggy.
            # So we have to export records ourselves
            #SeqIO.write(record,exportFiles[bestMatches[0]],"fastq")
            writeFASTQRecord(record,exportFiles[bestMatches[0]])
            print "Wrote record " +str(idx)+" "+ record.id + " to " + (exportFiles[bestMatches[0]]).name
    fh.close()
    # be nice and close the exported files
    for seq in aHandleList:
        exportFiles[seq].close()
    print "ready "

def main():
    args = parse_args()                         # parse commmand-line
    indexes = readYAMLConf(args.indexes)        # read sequences from yaml config file
    (forwardHnd,reverseHnd)= generateSeqHandles(indexes)    # generate forward and reverse complement indexes
    demultiplexFastqByBestMatch(args.reads, reverseHnd, int(args.mm),isForward=False)
    demultiplexFastqByBestMatch(args.reads, forwardHnd, int(args.mm),isForward=True)

if __name__ == '__main__':
    main()
