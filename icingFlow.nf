#!/usr/bin/env nextflow
/*
    Nextflow script for HLA typing from Oxford Nanopore amplicon reads
 */
if (!params.sample) {
  exit 1, "Please specify the sample FASTQ file containing 2D reads only."
}

if(!params.reference) {
  exit 1, "Please specify the genomic references."
}

fastq = file(params.sample)
base = fastq.getBaseName()

reference = file(params.reference)
referenceIdx = file(params.reference+".fai")
referenceAmb= file(params.reference+".amb")
referenceAnn= file(params.reference+".ann")
referenceBwt= file(params.reference+".bwt")
referenceFai= file(params.reference+".fai")
referencePac= file(params.reference+".pac")
referenceSa= file(params.reference+".sa")

process mapBWA {
    module 'bioinfo-tools'
    module 'samtools/1.3'
    module 'bwa'

    cpus 4

    input:
    file reference
    file referenceIdx
    file referenceAmb
    file referenceAnn
    file referenceBwt
    file referenceFai
    file referencePac
    file referenceSa

    output:
    file '*.bam' into HLAalignment
    file '*.bam' into HLAbam

    """
    bwa mem -x ont2d -t ${task.cpus} ${reference} ${fastq} | samtools view -bS -t ${reference} - | samtools sort - > ${base}.bam
    """
}


process selectCandidate {
    input:
    file HLAalignment

    output:
    file '*.candidates' into candidates

    """
    python ${workflow.launchDir}/selectCandidate.py -b ${HLAalignment} -r ${reference} > ${base}.candidates
    """
}

process extractReads {
    
    module 'samtools'

    input:
    file HLAbam
    file candidates

    output:
    file '*.cnd.fastq' into readsets

    """
    for f in `awk '{print \$1}' ${candidates}`; do echo \$f; samtools index ${HLAbam} ; samtools view -bh ${HLAbam} \${f}":" -o ${base}_\${f}.bam; picard-tools SamToFastq I=${base}_\${f}.bam F=${base}"_"\${f}".cnd.fastq"; done
    """
}
/*
    TODO:
    - add demultiplexing to the very beginning
    - Extract reads mapped to the candidate genomic references separately
    - feed each read set to cuna and generate consensuses
    - select best consensuses
    - black magic starts when you are finding the best match from IMGT/HLA to your consensus 
 */
