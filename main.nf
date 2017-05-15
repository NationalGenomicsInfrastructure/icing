#!/usr/bin/env nextflow
workflow.onError { // Display error message
    log.info "ICING: MinION HLA Genotyping"
    log.info "Usage (default values in [], aiming for Class-I typing ):"
    log.info "nextflow run main.nf --reducedRefDir /dir/to/HLA*fasta  --ref hla.dat --locus locus --sample file.fastq [ --minReadLength [2000] --minContigLength [3000] --threads [8] ]"
    log.info "Workflow execution stopped with the following message: " + workflow.errorMessage
}

fastq = file(params.sample)
base = fastq.getBaseName()
reducedRefDir = params.reducedRefDir
resultSuffix = "_"+base.replaceFirst(/.fastq/, "")+"_"+params.locus

//TODO: get rid of it dcBam = Channel.create()
//TODO: groi filteredPileup = Channel.create()
//TODO: groi genChannel = Channel.create()
//TODO: groi cdsChannel = Channel.create()
stepCount = 0

process mapWithALTcontigs {
    tag {params.locus + " " + params.sample + " ALT"}
	publishDir "ICING_" + incrementSteps() + "_ALT"+resultSuffix, mode: 'copy'
	module 'bwa/0.7.15'
    module 'samtools/1.3'

	cpus params.threads

	input: 
	set reads from fastq
	
	output:
	file 'rawALTmaps.bam' into rawALTbam
	file 'rawALTmaps.bam.bai' into rawALTbai

	script:
	"""
    set -eo pipefail
	LOCUS=${params.locus}
	
	# first make an ALT-aware index
	FIRSTALT=`ls ${reducedRefDir}/HLA*fasta|head -1`
	cat ${reducedRefDir}/HLA*fasta > alt\${LOCUS}.fasta
	bwa index \$FIRSTALT
	bwa mem -a \$FIRSTALT alt\${LOCUS}.fasta > alt\${LOCUS}.fasta.64.alt 
	bwa index -6 alt\${LOCUS}.fasta

	# now we have an index with ALT contigs; do mapping
    bwa mem -t ${task.cpus} -a -k 70 -W100 -r10 -A1 -B1 -O1 -E1 -L0 alt\${LOCUS}.fasta ${reads}|\
	samtools view --threads ${task.cpus} -bS -T alt\${LOCUS}.fasta -m ${params.minReadLength} - |\
    samtools sort --threads ${task.cpus} -  > rawALTmaps.bam
	samtools index rawALTmaps.bam
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
    file bam from rawALTbam

    output:
    file "filtered_${bam}" into filteredBam
    file "filtered_${bam}.bai" into filteredBai

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
process makePileUp {
    publishDir "ICING_" + incrementSteps() + "_pileup"+resultSuffix, mode:'copy'
    module 'samtools/1.3'

    input: 
    file fb from dcBam 

    output:
    file "${base}.pileup" into filteredPileup

    script:
    """
    samtools mpileup -B -d 1000 -Q 10 -A ${fb} > ${base}.pileup
    """
}

//////////////////////////////////////    expecting a pileup      ////////////////////////////////////////
//filteredPileup = filteredPileup.view{"Pileup from filtered reads: " + it}
process getConsensuses {
    publishDir "ICING_" + incrementSteps() + "_rawCons"+resultSuffix

    input:
    file pileup from filteredPileup

    output:
    file "cons.fasta" into consensusFASTA

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

    script:
    """
    python ${workflow.projectDir}/selectCandidate.py -f ${cf} -l ${params.minContigLength} > ${base}.candidate.fasta
    """
}

//////////////////////////////////////    expecting a consensus FASTA   ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// Now we have a list of candidates, try to find the most similar ones.

//process doGenotyping {
//    tag {params.locus + " " + ref + " " + query }
//
//    input:
//    set file(query) from candidatesFASTA
//
//    output:
//
//    when: 'genotype' in workflowSteps
//    script:
//    """
//	python genotypeHLA.py -c ${query} -d ${ref} -l ${params.locus}
//    """
//}

def incrementSteps() {
	return ++stepCount

}

