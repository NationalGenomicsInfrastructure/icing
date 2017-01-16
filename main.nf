#!/usr/bin/env nextflow
/*
    Nextflow script for HLA typing from Oxford Nanopore amplicon reads
		run as 
		nextflow run main.nf --reference A_gen.fasta --sample reads.fastq --minReadLength 1600
 */
if (!params.sample) {
  exit 1, "Please specify the sample FASTQ file containing 2D reads only."
}

if(!params.reference) {
  exit 1, "Please specify the genomic references."
}

if(!params.minReadLength) {
    exit 1, "Please provide minimal read length to consider at final consensus build (usually the half of the length of the expected reference is a good guess)"
}

if(!params.minContigLength || !params.minAvgQuality) {
    exit 1, "Please provide minimal contig length and minimal average quality (use the length of the locus and 40 as an initial guess) "
}

fastq = file(params.sample)
base = fastq.getBaseName()//.replaceFirst(/.fastq/, "")
locus = file(params.reference).getBaseName().replaceFirst(/.fasta/, "")

resultSuffix = "_"+base.replaceFirst(/.fastq/, "")+"_"+locus

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
    publishDir "rawAlignment"+resultSuffix

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
    file '*.bam' into HLAbam
    file '*.bai' into HLAbai

    """
    set -eo pipefail
    bwa mem -x ont2d -t ${task.cpus} ${reference} ${fastq} -a | \
    samtools view -bS -t ${reference} - | \
    samtools sort -  > ${base}.bam 
    samtools index ${base}.bam
    """
}

// when using pre-calculated BAM, start by adding --bam FILE.bam

fHLAbam = Channel.create()
rawHLAbam = Channel.create()

//HLAbam = file(params.bam)
//HLAbam.separate(fHLAbam, rawHLAbam){ a -> [a,a] }


// here we are getting rid of reads that are matching the reference, but are aligned with many mismatches (high edit distance)
// generally we are ignoring reads that are containing more mismatches then the 10% of the expected minimal contig size
editDist = params.minContigLength.toInteger()/10
process filterBestMatchingReads {
    publishDir "bestMatchingReads"+resultSuffix, mode:'copy'
    module 'samtools/1.3'

    input:
    file bam from HLAbam

    output:
    file "filtered_${bam}" into filteredBam

    script:
    """
    samtools view -h -m ${params.minReadLength} -bS ${bam} | bamtools filter -tag "AS:>1000" -in - -out "filtered_"${bam} 
    """
}

// Here we are running some stats on the bam files and trying to get only alleles
// with nice coverage. The rest of the crap is thrown away
process selectAllelesWithDecentCoverage {
    publishDir "goodCoverage"+resultSuffix
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
    publishDir "pileup"+resultSuffix, mode:'copy'
    module 'samtools/1.3'

    input: 
    file fb from dcBam 

    output:
    file "${base}.pileup" into filteredPileup

    script:
    """
    samtools mpileup -uf ${reference} ${fb} > ${base}.pileup
    """
}

process getConsensuses {
    publishDir "rawCons"+resultSuffix
    module 'bcftools/1.3'

    input:
    file pileup from filteredPileup

    output:
    file "cons.fastq" into consensusFASTQ

    script: 
    """
    #bcftools call -c ${pileup} --ploidy 1| vcfutils.pl vcf2fq > cons.fastq
    bcftools call -c ${pileup} | vcfutils.pl vcf2fq > cons.fastq
    """
} 

//consensusFASTQ = file("results/cons.fastq")

process selectConsensusCandidate {
    // export LD_LIBRARY_PATH=${HOME}/miniconda2/pkgs/libpng-1.6.22-0/lib/

    publishDir "candidates"+resultSuffix
    input:
    file cfq from consensusFASTQ

    output:
    file "*.candidate.fasta" into candidates

    """
    python ${workflow.projectDir}/selectCandidate.py -c ${cfq} -l ${params.minContigLength} -a ${params.minAvgQuality} 
    """
}
//candidates = candidates.view{"Candidates "+it}
// What if there is a single candidate? 



// All this magic is to calculate only the upper triangle of the distance matrix (since 
// it is symetrical and the main diagonals are all zeros).
//
// make two lists of candidates so we can have a Cartesian product and
// we will be able to calculate all the distances
(s1,s2,toSelectFrom) = candidates.flatten().into(3)
// make the cartesian, but leave out entries where the indexes are equal
doublePairs = s1.spread(s2).filter{it[0] != it[1]}
// now sort by filenames (so ["a", "b"] and ["b","a"] both will be ["a", "b"] )
// and store only unique entries
singlePairs = doublePairs.map { x -> x.sort() }.unique()
//singlePairs = singlePairs.view{"Single pairs to compare: "+ it}
process calculateSimilarity {
    publishDir "distanceMatrix"+resultSuffix, mode: 'copy'

    input:
    set file(seq1), file(seq2) from singlePairs 

    output:
    file "*_*.dist" into distances

    script:
    """
    needle -aformat score -datafile EDNAFULL -outfile stdout -gapopen 10.0 -gapextend 0.5 \
        -asequence $seq1 \
        -bsequence $seq2 |\
     awk '/HLA/{print }'|\
     tr -d "()" > $seq1"_"$seq2".dist"
    """
}

allDistances = distances.collectFile(){it -> "${it} "}

process selectMostDistantSequences {
    publishDir "mostDistant"+resultSuffix

    input:
    set file(dist) from allDistances

    output:
    file "sorted.candidates.txt" into sortedCandidateDistances
    file "most.distant.candidates.txt" into mostDistantCandidates 

    script:
    """
    for f in `cat $dist`; do
        cat \$f 
    done |\
    sort -k4n > sorted.candidates.txt
    tail -1 sorted.candidates.txt |awk '{printf("%s\\n%s",\$1,\$2)}' > most.distant.candidates.txt
    """
}
// now in the very last line of the "sorted.candidates.txt" file we have 
// the name of the two most distant candidates. Providing the consensuses are
// commensurable (they have similar length and quality) we are extracting reads that are
// - able to align to one of these most distant candidates (this should be OK 
//   by filtering reads only to aligned to the candidate reference - bwa reports only 3 alignments by default, so have to be careful)
// - are at least $minReadLength long (samtools view -m length(ref)/2)
// - alignment score is higher than 1000 (AS tag in bwa alignment) 

mostDistantCandidates = mostDistantCandidates.view{"mostDistantCandidates $it"}

//process extractReadsForBestCandidates {
//    publishDir "candidateReads"+resultSuffix, mode: 'copy'
//    module 'samtools/1.3'
//
//    input:
//    set file(alignment) from rawHLAbam
//    set file(cIds) from mostDistantCandidates
//
//    output:
//    file "*bam" into candidateReads
//    file "*bai" into candidateReadIndexes
//
//    script:
//    """
//    samtools index $alignment
//    for id in `cat $cIds`; do
//        samtools view -h $alignment "HLA:\${id}:" | samtools view -h -m ${params.minReadLength} -bS - | bamtools filter -tag "AS:>1000" -in - -out \${id}.bam
//        samtools index \${id}.bam
//    done
//    """
//}
