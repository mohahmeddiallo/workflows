#!/usr/bin/env nextflow

version = '0.10'

params.help = false
params.resume = false


log.info """

----------------------------------------------------------------------------
Pipeline for merging paired Sanger reads. ${version}
----------------------------------------------------------------------------
trace files                             : ${params.dataFolder}
unique identifier of forward read       : ${params.forwardReadID}
unique identifier of reverse read       : ${params.reverseReadID}
quality threshold for trimming          : ${params.qualityThreshold}
is there only one fragment of this gene : ${params.singleFragment}
merge two fragments together            : ${params.mergeFragments}
email                                   : ${params.email}
"""

// this ensures that the parameters are listed before the pipeline exits
if (params.help) exit 1

// a common mistake is to pass resume as a parameter of the workflow instead of a nextflow flag
// --resume is a parameter while -resume is the actual flag that activates resume
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

params.pairs   = "${params.dataFolder}/*_{${params.forwardReadID},${params.reverseReadID}}_*.ab1"


Channel
    .fromFilePairs(params.pairs, flat: true)
    .set {read_ch}

// Checking if the reads being merged are a fragment of a longer
// marker as is often the case with full-length 18S rDNA of nematodes
// the default is "YES", meaning the reads cover the whole region
if (params.singleFragment == "NO") {
   println("There is another fragment of this region")
   if (params.fragment != 1 && params.fragment != 2 ) exit 1, "Please specify --fragment=1  or --fragment=2"
} else if (params.singleFragment == "YES") {
   println("There is only one fragment of this region")
} 
else {
   exit 1, "Please specify YES or NO for single or paired ends"
}

// Convert ab1 files to fastq using seqret 
// for both forward and reverse reads
process ab1ToFastq {
    publishDir "$baseDir/01_fastqFiles"

    input:
    tuple val (sample_id), path (read1), path (read2) from read_ch

    output:
    path ("*${params.forwardReadID}_premixed.fastq") into fastq1_ch
    path ("*${params.reverseReadID}_premixed_tmp.fastq") into fastq2_ch

    script:
    """
    seqret -sformat abi -osformat fastq $read1 "${sample_id}_${params.forwardReadID}_premixed.fastq"
    seqret -sformat abi -osformat fastq $read2 "${sample_id}_${params.reverseReadID}_premixed_tmp.fastq"
    """
}

// Reverse complement the reverse read before merging using 
// seqkit, but first take the sample name from the read name
process reverseComplement {
    publishDir "$baseDir/01_fastqFiles"

    input:
    path (reverseRead) from fastq2_ch

    output:
    path("*${params.reverseReadID}_premixed.fastq") into fastq2_rc_ch

    script:
    """
    sample_id=\$(cut -d_ -f1 <<< $reverseRead)
    seqkit seq -p -r $reverseRead > "\${sample_id}_${params.reverseReadID}_premixed.fastq"
    """
}

// Trim the reads to high quality using seqtk with quality 
// threshold (-q) of 0.01, equivalent to Phred score of Q20
process qualityTrim {
    publishDir "$baseDir/02_qualityTrimmed"

    input:
    path(inputFoward) from  fastq1_ch
    path(inputReverse) from fastq2_rc_ch

    output:
    path("*${params.forwardReadID}*.fastq") into trimmed_forward
    path("*${params.reverseReadID}*.fastq") into trimmed_reverse

    script:
    """
    sample_id=\$(cut -d_ -f1 <<< $inputFoward)
    seqtk trimfq -q 0.01  $inputFoward > "\${sample_id}_${params.forwardReadID}.fastq"
    seqtk trimfq -q 0.01  $inputReverse > "\${sample_id}_${params.reverseReadID}.fastq"
    """
}


process mergeReads {
    publishDir "$baseDir/03_mergedReads"

    input:
    path(trimmed_read1) from trimmed_forward
    path(trimmed_read2) from trimmed_reverse

    output:
    path("*.fasta") into merged_fasta
    path("*merger") into merged_merger

    script:
    """
    sample_id=\$(cut -d_ -f1 <<< $trimmed_read1)
    merger -asequence $trimmed_read1 -bsequence $trimmed_read2 \
        -outseq "\${sample_id}_${params.fragment}.fasta" -outfile "\${sample_id}_${params.fragment}.megamerger"
    """
}


if (params.mergeFragments == "YES") {
   println("Continuing with merging two fragments")
}


Channel
    .fromFilePairs("$baseDir/03_mergedReads/*_{1,2}.fasta", flat: true)
    .set {fragment_ch}


process mergeFrags {
    publishDir "$baseDir/03_mergedReads"

    input:
    tuple val (sample_id), path (fragment1), path (fragment2) from fragment_ch

    output:
    path("*.fasta") into merged_frag_fasta
    path("*merger") into merged_frag_merger

    when:
    params.mergeFragments == "YES"

    script:
    """
    merger -asequence $fragment1 -bsequence $fragment2 \
        -outseq "${sample_id}.fasta" -outfile "${sample_id}.megamerger"
    """
}