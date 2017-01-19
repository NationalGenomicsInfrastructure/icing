if (!params.candidates) {
  exit 1, "Please specify the FASTA file containing consensus candidates"
}
if(!params.refDir) {
  exit 1, "Please specify a directory, where the genomic references are for the locus."
}

refChannel = Channel.fromPath( params.refDir + '/HLA*.fasta')

process makeSeparateConsensuses {
    input:
    file conss from file(params.candidates)

    output:
    file "cons*.fasta" into singleConsensuses

    """
    python ${workflow.projectDir}/separateFASTA.py -f ${conss}
    """
}
// make a Cartesian product of reference and consensus sequences
singleConsensuses = singleConsensuses.flatten()
refQueryPairs = refChannel.spread(singleConsensuses)
// calculate the delta files with mummer/nucmer
process calcDeltas {
    tag {params.locus + " " + ref + " " + query + " " + params.candidates}
    module 'mummer/3.23'

    input:
    set file(ref),file(query) from refQueryPairs

    output:
    file '*.delta' into deltas

    """
    rn=`basename ${ref}`
    qn=`basename ${query}`
    PREFIX=\${rn%.fasta}_\${qn%.fasta}
    nucmer --prefix=\${PREFIX} ${ref} ${query}
    """
}

// now calculate the mismatches (see http://mummer.sourceforge.net/manual/#nucmer for the delta file format)
// and sort deltas according to mismatches

deltas = deltas.flatten().toList()//.view{"D " + it}
process sortByMismatches {
    input:
    file fl from deltas

    output:
    file "sorted.deltas" into sortedDeltas

    script:
    """
    for f in ${fl}; do awk 'BEGIN{sum=0}{if(NF==7)sum+=\$5}END{print sum,FILENAME}' \$f; done| sort -n > sorted.deltas 
    """

}
