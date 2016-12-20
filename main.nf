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
    publishDir "rawAlignment"

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

//consensusFASTQ = file("rawCons/cons.fastq")
//HLAbam = file("rawAlignment/BC01.fastq.bam")
process selectConsensusCandidate {
    // export LD_LIBRARY_PATH=${HOME}/miniconda2/pkgs/libpng-1.6.22-0/lib/

    publishDir "candidates"
    input:
    file cfq from consensusFASTQ

    output:
    file "*.candidate.fasta" into candidates

    """
    python ${workflow.launchDir}/selectCandidate.py -c ${cfq} -l 2000 -a 50 
    """
}

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
process calculateSimilarity {
    publishDir "distanceMatrix", mode: 'copy'

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
    publishDir "mostDistant"

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
//   by filtering reads only to aligned to the candidate reference)
// - are at least $minReadLength long
// - are not soft-clipped (alignedLength > readLength * $minAlignedRatio where minAlignedRatio is between [0,1]
mostDistantCandidates = mostDistantCandidates.view{"mostDistantCandidates $it"}

process extractReadsForBestCandidates {
    publishDir "candidateReads"
    module 'samtools/1.3'

    input:
    set file(alignment) from HLAbam
    set file(cIds) from mostDistantCandidates

    output:
    file "*sam" into candidateReads

    script:
    """
    samtools index $alignment
    for id in `cat $cIds`; do
        samtools view -h $alignment "HLA:\${id}:" > \${id}.sam
    done
    # TODO: process SAM files, get rid of low QC reads
    """
}
