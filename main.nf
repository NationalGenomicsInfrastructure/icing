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
referenceAmb = file(params.reference+".amb")
referenceAnn = file(params.reference+".ann")
referenceBwt = file(params.reference+".bwt")
referenceFai = file(params.reference+".fai")
referencePac = file(params.reference+".pac")
referenceSa = file(params.reference+".sa")
outDirName = reference.getBaseName()

process mapBWA {
    module 'bioinfo-tools'
    module 'samtools/1.3'
    module 'bwa/0.7.15'

    cpus 6

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
    file '*.pileup' into HLApileup
    file '*.bam' into HLAbam

    """
    bwa mem -x ont2d -t ${task.cpus} ${reference} ${fastq} | \
    samtools view -bS -t ${reference} - | \
    samtools sort - | \
    tee ${base}.bam | \
    samtools mpileup -uf ${reference} - > ${base}.pileup
    """
}

process getConsensuses {
    publishDir "rawCons"
    module 'bcftools/1.3'

    input:
    file pileup from HLApileup 

    output:
    file "cons.fastq" into consensusFASTQ

    script: 
    """
    bcftools call -c ${pileup} --ploidy 1| vcfutils.pl vcf2fq > cons.fastq
    """
} 

process selectConsensusCandidate {
    // export LD_LIBRARY_PATH=${HOME}/miniconda2/pkgs/libpng-1.6.22-0/lib/

    publishDir "candidates"
    input:
    file cfq from consensusFASTQ

    output:
    file "${base}.candidates.fasta" into candidates

    """
    echo "Hajra"
    python ${workflow.launchDir}/selectCandidate.py -c ${cfq} -l 2000 -a 50 > ${base}.candidates.fasta
    """
}
//
//process extractReads {
//    
//    module 'samtools'
//
//    input:
//    file HLAbam
//    file candidates
//
//    output:
//    file '*.cnd.fastq' into readsets
//
//    /*
//        This relatively long line of script 
//        - generates index for the resulting BAM
//        - gets the candidate names from the previous result file
//        - extracts reads into a BAM that are mapped to the candidate
//        - extracts reads into a FASTQ from the BAM
//        We should get a FASTQ for each candidate
//     */
//    """
//    samtools index ${HLAbam}; for f in `awk '{print \$1}' ${candidates}`; do echo \$f;  samtools view -bh ${HLAbam} \${f}":" -o ${base}_\${f}.bam; picard-tools SamToFastq I=${base}_\${f}.bam F=${base}"_"\${f}".cnd.fastq"; done
//    """
//}
//
//// generating consensuses, and saving parts containing reads
//// the funny bash part is to get rid of failed consensus building runs
//process generateConsensuses {
//    publishDir "ContigsAndReads_" + outDirName
//
//    input:
//    file reads from readsets
//
//    output:
//    file '*fasta' into consensuses 
//
//    shell:
//    """
//    samtools mpileup -uf ../ON0526/IMGT/A_gen.fasta 1_All_BC01.bam > BC01.mpileup
//    """
//}
//
//// Now we are selecting consensuses, where there are more reads assigned to the result
//process pruneContigs {
//    publishDir "ContigsAndReads_" + outDirName
//
//    input:
//    file cont from consensuses
//    
//    output:
//    file '*.pruned.fasta' into prunedContigs
//
//    """
//    # 
//    awk '/>/ && \$3 !~/reads=[1-9]\$/{print FILENAME,\$0}' ${consensuses} |awk -F"[:.]" '{print \$3}'| sort -u    
//
//    """
//}
//
/*
    TODO:
    - select best consensuses
    - black magic starts when you are finding the best match from IMGT/HLA to your consensus 
    - add demultiplexing to the very beginning
 */
