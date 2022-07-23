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
is there another fragment of this gene  : ${params.singleFragment}
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

if (params.singleFragment == "NO") {
   println("There is another fragment of this region")
   if (params.fragment != 1 && params.fragment != 2 ) exit 1, "Please specify --fragment=1  or --fragment=2"
} else if (params.singleFragment == "YES") {
   println("There is only one fragment of this region")
} 
else {
   exit 1, "Please specify YES or NO for single or paired ends"
}


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

    script:
    """
    sample_id=\$(cut -d_ -f1 <<< $trimmed_read1)
    merger -asequence $trimmed_read1 -bsequence $trimmed_read2 \
        -outseq "\${sample_id}_${params.fragment}.fasta" -outfile "\${sample_id}_${params.fragment}.megamerger"
    """
}


