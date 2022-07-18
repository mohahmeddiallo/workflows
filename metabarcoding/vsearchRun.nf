#!/usr/bin/env nextflow


version = '0.1'


params.help = false
params.resume = false


log.info """

---------------------------------------------------------------
Amplicon Sequence Analysis Pipeline Using VSEARCH v. ${version}
---------------------------------------------------------------
reads                           : ${params.rawreads}
reference database              : ${params.ref}
minimum overlap length          : ${params.minOvLen}
maximum overlap difference      : ${params.maxDiff}
maximum expected error          : ${params.maxEE}
fastq minimus length            : ${params.fastqMinLen}
fastq maximum length            : ${params.fastqMaxLen}
fastq maximum number of Ns      : ${params.fastqmaxNs}
minimum unique size             : ${params.minUnSize}
precluster ID                   : ${params.PreClusterID}
cluster ID                      : ${params.clusterID}
sintax cut-off                  : ${params.sintaxCutoff}
email                           : ${params.email}
"""



Channel
    .fromFilePairs(params.rawreads, flat: true)
    .set{read_ch}


Channel
    .fromPath(params.ref)
    .into{ref_for_chimera_ch; ref_for_sintaxTax_ch}



process mergeReads {
    label "insideContainer"
    

    input:
    tuple val(sample_id), path(read1), path(read2) from read_ch


    output:
    path ("*merged.fastq") into (merged_for_statistics_ch, merged_for_filtering_ch)


    script:
    """
    echo 
    echo ====================================
    echo Processing sample $sample_id
    echo ====================================
    vsearch --threads $params.cpus \
        --fastq_mergepairs $read1 \
        --reverse $read2 \
        --fastq_minovlen $params.minOvLen \
        --fastq_maxdiffs $params.maxDiff \
        --fastqout "$sample_id".merged.fastq \
        --fastq_eeout
    """
}



process mergeStats {
    label "insideContainer"


    input:
    file (merged_fastq) from merged_for_statistics_ch


    output:
    path ("*.stats")

    script:
    """
    sample_id=`cut -d_ -f1 <<< $merged_fastq`
    echo 
    echo ====================================
    echo Calculate quality statistics
    echo ====================================
    vsearch --threads $params.cpus \
        --fastq_eestats $merged_fastq \
        --output "\${sample_id}".stats
    """
}

process qualityFiltering {
    label "insideContainer"


    input:
    file(merged_fastq2) from merged_for_filtering_ch


    output:
    path ("*.filtered.fasta") into filtered_ch


    script:
    """
    sample_id=`cut -d. -f1 <<< $merged_fastq2`
    echo 
    echo ====================================
    echo Quality filtering
    echo ====================================
    vsearch --threads $params.cpus \
        --fastq_filter $merged_fastq2 \
        --fastq_maxee $params.maxEE \
        --fastq_minlen $params.fastqMinLen \
        --fastq_maxlen $params.fastqMaxLen \
        --fastq_maxns $params.fastqmaxNs \
        --fastaout "\${sample_id}".filtered.fasta \
        --fasta_width $params.fastaWidth
        """

}

process derepSample {
    label "insideContainer"


    input:

    file(filtered_fasta) from filtered_ch


    output:
    path("*.derep.fasta") into derep_ch
    path("*.derep.uc") into derep_uc_ch


    script:
    """
    sample_id=\$(cut -d. -f1 <<< $filtered_fasta)
    echo =====================================================
    echo Dereplicating relable with sample_n
    echo =====================================================
    vsearch --threads $params.cpus \
        --derep_fulllength $filtered_fasta \
        --strand plus \
        --output "\${sample_id}".derep.fasta \
        --sizeout \
        --uc "\${sample_id}".derep.uc \
        --relabel "\${sample_id}". \
        --fasta_width $params.fastaWidth
    """
}


process mergeAllSamples {
    label "insideContainer"

    input:
    path "*" from derep_ch.collect()

    output:
    file "all.fasta" into (concat_all_ch, concat_all_for_map_ch)

    script:
    """
    cat * >> all.fasta
    """
}


process derepAllSamples {
    label "insideContainer"

    input:
    file (all_fasta) from concat_all_ch

    output:
    path("all_derep.fasta") into (derep_all_ch, derep_all_for_map_ch)
    path("all_derep.uc") into all_derep_uc_ch

    script:
    """
    vsearch --threads $params.cpus \
        --derep_fulllength $all_fasta \
        --minuniquesize $params.minUnSize \
        --sizein \
        --sizeout \
        --fasta_width $params.fastaWidth \
        --uc all_derep.uc \
        --output all_derep.fasta
    """
}


process preCluster {
    label "insideContainer"

    input:
    file(all_derep_fasta) from derep_all_ch

    output:
    path("all.preclustered.fasta") into preclustered_ch
    path("all.preclustered.uc") into all_preclustered_uc_ch

    script:
    """
    vsearch --threads $params.cpus \
        --cluster_size $all_derep_fasta \
        --id $params.PreClusterID \
        --strand plus \
        --sizein \
        --sizeout \
        --fasta_width $params.fastaWidth \
        --uc all.preclustered.uc \
        --centroids all.preclustered.fasta
    """
}


process denovoChimDetection {
    label "insideContainer"

    input:
    file(preclustered) from preclustered_ch

    output:
    path("all.denovo.nonchimeras.fasta") into denovo_non_chimera_ch

    script:
    """
    vsearch --threads $params.cpus \
        --uchime_denovo $preclustered \
        --sizein \
        --sizeout \
        --fasta_width $params.fastaWidth \
        --nonchimeras all.denovo.nonchimeras.fasta
    """
}


process referenceChimDetection {
    label "insideContainer"

    input:
    file(nonchimera) from denovo_non_chimera_ch
    file(ref_database) from ref_for_chimera_ch

    output:
    path("all.ref.nonchimeras.fasta") into reference_non_chimera_ch


    script:
    """
    vsearch --threads $params.cpus \
        --uchime_ref all.denovo.nonchimeras.fasta \
        --db $ref_database \
        --sizein \
        --sizeout \
        --fasta_width $params.fastaWidth \
        --nonchimeras all.ref.nonchimeras.fasta
    """
}

process extractAllNonChimerasStep1 {
    label "locally"

    input:
    file(all_derep) from derep_all_for_map_ch
    file(all_preclustered) from all_preclustered_uc_ch
    file(all_ref_non_chimera) from reference_non_chimera_ch
    
    
    output:
    path("all.nonchimeras.derep.fasta") into all_nonchim_derep_ch


    script:
    """
    map.pl $all_derep $all_preclustered $all_ref_non_chimera > all.nonchimeras.derep.fasta
    """

}


process extractAllNonChimerasStep2 {
    label "locally"

    input:
    file(all_fasta) from concat_all_for_map_ch
    file(all_derep_uc) from all_derep_uc_ch
    file(all_nonChime_derep) from all_nonchim_derep_ch

    output:
    path("all.nonchimeras.fasta") into all_non_chimeras_ch

    script:
    """
    map.pl $all_fasta $all_derep_uc $all_nonChime_derep > all.nonchimeras.fasta
    """
}


process cluster {
    label "insideContainer"

    input:
    file(all_non_chimera) from all_non_chimeras_ch

    output:
    path("all.otus.fasta") into otus_ch
    path("all.clustered.uc") into clustered_uc_ch
    path("all.otutab.txt") into otutab_ch

    script:
    """
    vsearch --threads $params.cpus \
        --cluster_size $all_non_chimera \
        --id $params.clusterID \
        --strand plus \
        --sizein \
        --sizeout \
        --fasta_width $params.fastaWidth \
        --uc all.clustered.uc \
        --relabel OTU_ \
        --centroids all.otus.fasta \
        --otutabout all.otutab.txt
    """
}

process sintaxTaxonomy {
    label "insideContainer"

    input:
    file(otus) from otus_ch
    file(ref_database) from ref_for_sintaxTax_ch
    

    output:
    path("all.otus.tax.txt")

    script:
    """
    vsearch --sintax $otus \
        --db $ref_database \
        --tabbedout all.otus.tax.txt \
        --sintax_cutoff $params.sintaxCutoff
    """
}
