#!/usr/bin/env nextflow
workflow.onError { // Display error message
    log.info "ICING: MinION HLA Genotyping"
    log.info "Usage (default values in [], aiming for Class-I typing ):"
    log.info "nextflow run main.nf --steps align --genDir dir --locus locus --sample file.fastq --minReadLength [2000] --minContigLength [3000]"
    log.info "or"
    log.info "nextflow run main.nf --steps genotype --cdsDir dir --bam file.bam --minCoverage [70] --minContigLength [3000] --locus locus --sampleSuffix mySample"
    log.info "steps can be merged like"
    log.info "nextflow run main.nf --steps align,genotype --genDir dir --cdsDir dir --locus locus --sample file.fastq --minReadLength [2000] --minCoverage [70] --minContigLength [3000]"
    log.info "Optionally number of working threads can be assigned as --threads [8]"
    log.info "Workflow execution stopped with the following message: " + workflow.errorMessage
}

fastq = ''
base = ''
resultSuffix = ''
dcBam = Channel.create()
filteredPileup = Channel.create()
genChannel = Channel.create()
cdsChannel = Channel.create()
stepCount = 0

workflowSteps = params.steps.split(',').collect {it.trim()}
if ('align' in workflowSteps) {
    if(!params.sample) { exit 1, "Please specify the sample FASTQ file containing 2D reads only." }
    if(!params.genDir) { exit 1, "Please specify a directory, where the genomic references are for the locus." }
    fastq = file(params.sample)
    base = fastq.getBaseName()//.replaceFirst(/.fastq/, "")
    resultSuffix = "_"+base.replaceFirst(/.fastq/, "")+"_"+params.locus
    genChannel = Channel.fromPath( params.genDir + '/HLA*.fasta')
} else if ('genotype' in workflowSteps ) {
    if(!params.cdsDir) { exit 1, "Please specify a directory, where the CDS references are for the locus." }
    base = file(params.bam).getBaseName()
    resultSuffix = "_"+params.sampleSuffix +"_"+params.locus
    cdsChannel = Channel.fromPath( params.cdsDir + '/cdsHLA*.fasta')
} else {
    exit 1, 'Please provide \'align\', or \'genotype\' as one of the workflow steps (--steps align,genotype)'
}

process mapBWA {
    tag {params.locus + " " + ref + " " + fastq}

    module 'bioinfo-tools'
    module 'samtools/1.3'
    module 'bwa/0.7.15'

    cpus params.threads

    input:
    file ref from genChannel

    output:
    file 'HLAXX*.bam' into HLAbam

    when: 'align' in workflowSteps
    script:
    """
    set -eo pipefail
    bwa index ${ref}
    bwa mem -t ${task.cpus} -a -k 70 -W100 -r10 -A1 -B1 -O1 -E1 -L0 ${ref} ${fastq}|  samtools view --threads ${task.cpus} -bS -T ${ref} -m ${params.minReadLength} - |
    samtools sort --threads ${task.cpus} -  > HLAXX`tr -cd '[:alnum:]' < /dev/urandom | fold -w30 | head -n1`.bam
    """
}
// we are not using genChannel anymore: close it or it will hang and wait at genotyping
genChannel.close()

process mergeBAMs {
    publishDir "ICING_" + incrementSteps() + "_rawAlignment"+resultSuffix, mode: 'copy'

    module 'samtools/1.3'
   
    input:
    file singleBAM from HLAbam.toList()

    output:
    file "merged.bam" into mergedBAM 

    when: 'align' in workflowSteps
    script:
    """
    samtools merge --threads ${task.cpus} merged.bam HLAXX*bam
	samtools index merged.bam
    """
}
// here we are getting rid of reads that are matching the reference, but are aligned with many mismatches (high edit distance)
// generally we are ignoring reads that are containing more mismatches then the 5% of the expected minimal contig size
editDistance = params.minContigLength.toInteger()*0.05

process filterBestMatchingReads {
    publishDir "ICING_" + incrementSteps() + "_bestMatchingReads"+resultSuffix, mode:'copy'
    module 'samtools/1.3'
    module 'bamtools/2.4.1'

    input:
    file bam from mergedBAM

    output:
    file "filtered_${bam}" into filteredBam
    file "filtered_${bam}.bai" into filteredBai

    when: 'align' in workflowSteps
    script:
    """
    samtools view -h -bS ${bam} | bamtools filter -tag "NM:<${editDistance}" -in - -out "filtered_"${bam} 
	samtools index "filtered_"${bam}
    """
}

// Here we are running some stats on the bam files and trying to get only alleles
// with nice coverage. The rest of the crap is thrown away
process selectAllelesWithDecentCoverage {
    publishDir "ICING_" + incrementSteps() + "_goodCoverage"+resultSuffix
    module 'samtools/1.3'

    input:
    file fbam from filteredBam

    output:
    file "candidate.stats" into stats
    file "mergedDC.bam" into dcBam
    file "mergedDC.bam.bai"
    
    when: 'align' in workflowSteps
    script:
    """
    # - get the allele names from the BAM header, 
    # - calculate coverage depth for each position of each allele, but only using reads that are at least minReadLength long
    # - calculate statistics, and select the best 16 alleles (no point to get more)
    # - create a BAM file that contains only long reads and the best alleles
    SCRIPTDIR=`dirname ${workflow.scriptFile}`
    samtools index ${fbam}
    for allele in `samtools view -H ${fbam} | awk '/^@SQ/{print \$2}'| sed 's/SN://g'`; do 
                samtools depth -a -l ${params.minReadLength} ${fbam} -r \${allele}: > \${allele}.coverage; 
                python \${SCRIPTDIR}/binCoverage.py -f \${allele}.coverage; 
    done | sort -k3n | tail -16 > candidate.stats
    for allele in `awk '{print \$1}' candidate.stats`; do 
        samtools view -hb ${fbam} \${allele}: > FHB\${allele}.bam
    done
    samtools merge mergedDC.bam FHB*bam
	samtools index mergedDC.bam
    """
}
// END OF ALIGNMENT
//////////////////////////////////////    expecting a BAM        ////////////////////////////////////////
if('genotype' in workflowSteps) {
    if(params.bam) {
        dcBam = Channel.fromPath(params.bam)
    }
} else {
    dcBam.close()
}

process makePileUp {
    publishDir "ICING_" + incrementSteps() + "_pileup"+resultSuffix, mode:'copy'
    module 'samtools/1.3'

    input: 
    file fb from dcBam 

    output:
    file "${base}.pileup" into filteredPileup

    when: 'genotype' in workflowSteps
    script:
    """
    samtools mpileup -B -d 1000 -Q 10 -A ${fb} > ${base}.pileup
    """
}

//////////////////////////////////////    expecting a pileup      ////////////////////////////////////////
if('genotype' in workflowSteps) {
    if(params.pileup) {
        filteredPileup = Channel.fromPath(params.pileup)
    }
} else {
    filteredPileup.close()
}
filteredPileup = filteredPileup.view{"Pileup from filtered reads: " + it}
process getConsensuses {
    publishDir "ICING_" + incrementSteps() + "_rawCons"+resultSuffix

    input:
    file pileup from filteredPileup

    output:
    file "cons.fasta" into consensusFASTA

    when: 'genotype' in workflowSteps
    script: 
    """
    python ${workflow.projectDir}/makeConsensusFromPileup.py -p ${pileup} -d ${params.minCoverage} > cons.fasta 
    """
} 

process selectConsensusCandidates {
    publishDir "ICING_" + incrementSteps() + "_candidates"+resultSuffix

    input:
    file cf from consensusFASTA

    output:
    file "*.candidate.fasta" into candidatesFASTA

    when: 'genotype' in workflowSteps
    script:
    """
    python ${workflow.projectDir}/selectCandidate.py -f ${cf} -l ${params.minContigLength} > ${base}.candidate.fasta
    """
}

//////////////////////////////////////    expecting a consensus FASTA   ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// Now we have a list of candidates, try to find the most similar ones.
// Note, we are using the CDS only FASTAs, as we want to get the genotypes,
// but this also means that the process provides only 3 fields precision for the candidate list.
// To get 4 fields precision you have to go through the candidate list in an alignment browser like IGV one-by-one
// We are assuming the cdsHLA*.fasta files were generated from the ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/*_nuc.fasta
// files, containing CDSs
// ditto, we are assigning IDs and filename to consensuses as "cons"
process makeSeparateConsensuses {
    input:
    file conss from candidatesFASTA

    output:
    file "cons*.fasta" into singleConsensuses

    when: 'genotype' in workflowSteps
    script:
    """
    python ${workflow.projectDir}/separateFASTA.py -f ${conss} -p cons
    """
}
singleConsensuses = singleConsensuses.view{"Consensuses to use: $it"}

process doGenotyping {
    tag {params.locus + " " + ref + " " + query }

    input:
    set file(ref),file(query) from consPairs 

    output:

    when: 'genotype' in workflowSteps
    script:
    """
    """
}

def incrementSteps() {
	return ++stepCount

}

// making a sam with SW from exons 2&3
// szilva@galatea ~/sr/minion/hlapool161222/prep/redA/ICING_6_candidates_bc1_HLA-A $ /home/szilva/dev/Complete-Striped-Smith-Waterman-Library/src/pyssw.py -f 1000 -rcs -l /home/szilva/dev/Complete-Striped-Smith-Waterman-Library/src/  -e 0 candidates.fasta ../exons.fasta | tee CSSWL.sam
// awk '$3 ~/consHLA00005/{print }' CSSWL.sam > consHLA00005.sam
// 
 
