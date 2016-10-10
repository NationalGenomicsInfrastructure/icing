#!/usr/bin/env nextflow
/*
    Nextflow script for HLA typing from Oxford Nanopore amplicon reads
		run as 
		nextflow run icingFlow.nf --reference A_gen.fasta --sample reads.fastq
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
    module 'bwa/0.7.15'

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

    /*
        This relatively long line of script 
        - generates index for the resulting BAM
        - gets the candidate names from the previous result file
        - extracts reads into a BAM that are mapped to the candidate
        - extracts reads into a FASTQ from the BAM
        We should get a FASTQ for each candidate
     */
    """
    samtools index ${HLAbam}; for f in `awk '{print \$1}' ${candidates}`; do echo \$f;  samtools view -bh ${HLAbam} \${f}":" -o ${base}_\${f}.bam; picard-tools SamToFastq I=${base}_\${f}.bam F=${base}"_"\${f}".cnd.fastq"; done
    """
}

process generateConsensuses {
    input:
    file reads from readsets

    output:
    file '*fasta' into consensuses 

    shell:
    """
    for f in ${reads}; do (/home/szilva/dev/canu/Linux-amd64/bin/canu useGrid=false -p \${f}.canu -d \${f}.canu genomeSize=15000 -nanopore-raw \${f} && cp \${f}.canu/\${f}.canu.unassembled.fasta \${f}.fasta ) || echo \$f" failed to assemble"; done
    """
}
/*
    TODO:
    - feed each read set to cuna and generate consensuses
    - select best consensuses
    - black magic starts when you are finding the best match from IMGT/HLA to your consensus 
    - add demultiplexing to the very beginning
 */
