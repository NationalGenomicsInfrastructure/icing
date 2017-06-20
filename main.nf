#!/usr/bin/env nextflow
workflow.onError { // Display error message
    log.info "ICING: MinION HLA Genotyping"
    log.info "Usage (default values in [], aiming for Class-I typing):"
    log.info "nextflow run main.nf --genomicRef A_gen.fasta  --ref hla.dat --locus locus --sample file.fastq [ --minPileupVote [40] --minReadLength [2000] --minContigLength [3000] --threads [8] --bam <alignments.bam>]"
    log.info "Workflow execution stopped with the following message: " + workflow.errorMessage
}


if(!params.minPileupVote) {params.minPileupVote = 40}
if(!params.minReadLength) {params.minReadLength = 2000}
if(!params.minContigLength) {params.minContigLength = 3000}
if(!params.editDistance ) {params.editDistance = params.minContigLength.toInteger()*0.05}
println params

fastq = file(params.sample)
base = fastq.getBaseName()
genomicRef = file(params.genomicRef)
resultSuffix = "_"+base.replaceFirst(/.fastq/, "")+"_"+params.locus
ref = file(params.ref)

// to save consecutive step statuses (alignment, consensuses, stats) in ICING_... directories
stepCount = 0
rawALTbam = Channel.create()
hasBam = params.bam
hasEditDistance = params.editDistance

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

	// if we do not have a BAM with aligned reads, run alignment
	when: !hasBam
	script:
	"""
    set -eo pipefail
	LOCUS=${params.locus}
	
	# first make an ALT-aware index
	cp ${params.genomicRef} alt\${LOCUS}.fasta
	# get the very first genomic sequence
	FIRSTALT=`awk '/>/{print \$1; exit 0}' alt\${LOCUS}.fasta |sed 's/>//'`".fasta"
	awk '/>/{print;getline;while(\$1!~/>/){getline;print}exit 0}' alt\${LOCUS}.fasta> \$FIRSTALT
	bwa index \$FIRSTALT
	bwa mem -a \$FIRSTALT alt\${LOCUS}.fasta > alt\${LOCUS}.fasta.64.alt 
	bwa index -6 alt\${LOCUS}.fasta

	# now we have an index with ALT contigs; do mapping
    #bwa mem -t ${task.cpus} -x ont2d alt\${LOCUS}.fasta ${reads}|
    bwa mem -t ${task.cpus} -a -k 70 -W100 -r10 -A1 -B1 -O1 -E1 -L0 alt\${LOCUS}.fasta ${reads}|\
	samtools view --threads ${task.cpus} -bS -T alt\${LOCUS}.fasta -m ${params.minReadLength} - |\
    samtools sort --threads ${task.cpus} -  > rawALTmaps.bam
	samtools index rawALTmaps.bam
	"""
}

// here we are getting rid of reads that are matching the reference, but are aligned with many mismatches (high edit distance)
// generally we are ignoring reads that are containing more mismatches then the 5% of the expected minimal contig size
if(hasBam) {
	rawALTbam = Channel.fromPath(params.bam)
}

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
    samtools view -h -bS ${bam} | bamtools filter -tag "NM:<${params.editDistance}" -in - -out "filtered_"${bam} 
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
    #done | sort -k3n > candidate.stats
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
    python ${workflow.projectDir}/makeConsensusFromPileup.py -p ${pileup} -d ${params.minPileupVote} > cons.fasta 
    """
} 

process selectConsensusCandidates {
    publishDir "ICING_" + incrementSteps() + "_candidates"+resultSuffix

    input:
    file cf from consensusFASTA

    output:
    file "*.candidates.fasta" into candidatesFASTA

    script:
    """
    python ${workflow.projectDir}/selectCandidate.py -f ${cf} -l ${params.minContigLength} > ${base}.candidates.fasta
    """
}

//////////////////////////////////////    expecting a consensus FASTA   ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// Now we have a list of candidates, try to find the most similar ones.

process doGenotyping {
    tag {params.locus + " " + ref }
    publishDir "ICING_" + incrementSteps() + "_genotypes"+resultSuffix

    input:
    set file(query) from candidatesFASTA

    output:
	file "typing.log"
	file 'genotypes.txt' into types

    script:
    """
	python ${workflow.projectDir}/genotypeHLA.py -c ${query} -d ${ref} -l ${params.locus} > typing.log
	GT=genotypes.txt
	echo "#################################################################################################" > \$GT
	echo "">>\$GT
	echo "HLA types from sample ${params.sample} and locus ${params.locus}: ">> \$GT
	echo "">>\$GT
	grep ^HLA typing.log | sort -u >> \$GT
	echo "">>\$GT
	echo "#################################################################################################" >> \$GT
    """
}
types.subscribe {println it.text}

workflow.onComplete { // Display complete message
		log.info "ICING - MinION HLA typing" 
		log.info "Command Line: $workflow.commandLine"
		log.info "Project Dir : $workflow.projectDir"
		log.info "Launch Dir  : $workflow.launchDir"
		log.info "Work Dir    : $workflow.workDir"
		log.info "Completed at: $workflow.complete"
		log.info "Duration    : $workflow.duration"
		log.info "Success     : $workflow.success"
		log.info "Exit status : $workflow.exitStatus"
}

def incrementSteps() {
	return ++stepCount
}

