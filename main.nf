#!/usr/bin/env nextflow
/*
    Nextflow script for HLA typing from Oxford Nanopore amplicon reads
		run as 
		nextflow run main.nf --reference A_gen.fasta --sample reads.fastq --minReadLength 1600
 */
if (!params.sample) {
  exit 1, "Please specify the sample FASTQ file containing 2D reads only."
}

if(!params.refDir) {
  exit 1, "Please specify a directory, where the genomic references are for the locus."
}

if(!params.minReadLength) {
    exit 1, "Please provide minimal read length to consider at final consensus build (usually the half of the length of the expected reference is a good guess)"
}

if(!params.minContigLength ) {
    exit 1, "Please provide minimal contig length (use the length of the locus as an initial guess) "
}

fastq = file(params.sample)
base = fastq.getBaseName()//.replaceFirst(/.fastq/, "")
resultSuffix = "_"+base.replaceFirst(/.fastq/, "")+"_"+params.locus
refChannel = Channel.fromPath( params.refDir + '/HLA*.fasta')

process mapBWA {
    tag {params.locus + " " + ref + " " + fastq}

    module 'bioinfo-tools'
    module 'samtools/1.3'
    module 'bwa/0.7.15'

    cpus 6

    input:
    file ref from refChannel

    output:
    file 'HLAXX*.bam' into HLAbam

    """
    set -eo pipefail
    bwa index ${ref}
    bwa mem -x ont2d -a -t ${task.cpus} ${ref} ${fastq}|  samtools view --threads ${task.cpus} -bS -T ${ref} -m ${params.minReadLength} - |\
    samtools sort --threads ${task.cpus} -  > HLAXX`tr -cd '[:alnum:]' < /dev/urandom | fold -w30 | head -n1`.bam
    """
}

process mergeBAMs {
    publishDir "ICING_rawAlignment"+resultSuffix, mode: 'copy'

    module 'samtools/1.3'
   
    input:
    file singleBAM from HLAbam.toList()

    output:
    file "merged.bam" into mergedBAM 

    """
    samtools merge merged.bam HLAXX*bam
    """
}
//
//// when using pre-calculated BAM, start by adding --bam FILE.bam
//
//fHLAbam = Channel.create()
//rawHLAbam = Channel.create()
//
////HLAbam = file(params.bam)
////HLAbam.separate(fHLAbam, rawHLAbam){ a -> [a,a] }
//
//
//// here we are getting rid of reads that are matching the reference, but are aligned with many mismatches (high edit distance)
//// generally we are ignoring reads that are containing more mismatches then the 10% of the expected minimal contig size
editDistance = params.minContigLength.toInteger() - params.minContigLength.toInteger()/10
process filterBestMatchingReads {
    publishDir "ICING_bestMatchingReads"+resultSuffix, mode:'copy'
    module 'samtools/1.3'

    input:
    file bam from mergedBAM

    output:
    file "filtered_${bam}" into filteredBam

    script:
    """
    samtools view -h -bS ${bam} | bamtools filter -tag "AS:>${editDistance}" -in - -out "filtered_"${bam} 
    """
}

// Here we are running some stats on the bam files and trying to get only alleles
// with nice coverage. The rest of the crap is thrown away
process selectAllelesWithDecentCoverage {
    publishDir "ICING_goodCoverage"+resultSuffix
    module 'samtools/1.3'

    input:
    file fbam from filteredBam

    output:
    file "candidate.stats" into stats
    file "mergedDC.bam" into dcBam
    
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
    samtools merge -@8 mergedDC.bam FHB*bam
    """
}

process makePileUp {
    publishDir "ICING_pileup"+resultSuffix, mode:'copy'
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

process getConsensuses {
    publishDir "ICING_rawCons"+resultSuffix

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
    publishDir "ICING_candidates"+resultSuffix

    cpus 6

    input:
    file cf from consensusFASTA

    output:
    file "*.candidate.fasta" into candidatesFASTA

    """
    python ${workflow.projectDir}/selectCandidate.py -f ${cf} -l ${params.minContigLength} > ${base}.candidate.fasta
    """
}

// Now we have a list of candidates, try to find the most similar ones.
// Note, we are using the CDS only FASTAs, as we want to get the genotypes,
// but this also means that the process provides only 3 fields precision for the candidate list.
// To get 4 fields precision you have to go through the candidate list in an alignment browser like IGV one-by-one
// We are assuming the cdsHLA*.fasta files were generated from the ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/*_nuc.fasta
// files, containing CDSs
// ditto, we are assigning IDs and filename to consensuses as "cons"
cdsChannel = Channel.fromPath(params.refDir + '/cdsHLA*.fasta')
process makeSeparateConsensuses {
    input:
    file conss from candidatesFASTA

    output:
    file "cons*.fasta" into singleConsensuses

    """
    python ${workflow.projectDir}/separateFASTA.py -f ${conss} -p cons
    """
}
// make a Cartesian product of reference and consensus sequences
singleConsensuses = singleConsensuses.flatten()
refQueryPairs = cdsChannel.spread(singleConsensuses)
// calculate the delta files with mummer/nucmer
process calcDeltas {
    tag {params.locus + " " + ref + " " + query }
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
    publishDir "ICING_final_candidate_list"+resultSuffix, mode: 'copy'
    input:
    file fl from deltas

    output:
    file "sorted.deltas" into sortedDeltas

    script:
    """
    for f in ${fl}; do awk 'BEGIN{sum=0}{if(NF==7)sum+=\$5}END{print sum,FILENAME}' \$f; done| sort -n > sorted.deltas 
    """
}
