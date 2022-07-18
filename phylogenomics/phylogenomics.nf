#!/usr/bin/env nextflow


version = '0.1'

// params.phypyPrunerMinTax="$(($(ls -1 $(dirname $params.proteins) | wc -l)/2))"
params.help = false
params.resume = false

// calculating the default value for minTax (Half the number of taxa)
// number of proteomes in the input directory
proteinCount = new File(params.proteomesPath).listFiles().count { it.name ==~ /.*\.faa/ }

// interger value of the proteome number divided by 2
params.minTax = proteinCount.intdiv(2)

log.info """

----------------------------------------------------------------------------
Phylogenetic analysis pipeline. ${version}
----------------------------------------------------------------------------
protein fasta                           : ${params.proteins}
proteinortho min. reciprocal similarity : ${params.proteinortho_sim}
proteinortho min. percent identity      : ${params.proteinortho_id}
proteinortho E value                    : ${params.proteinortho_evalue}
proteinortho results directory          : ${params.proteinorthoDir}
max number of iterative refinement      : ${params.mafftmaxIterate}
iqtree bootstrap                        : ${params.iqtreeBS}
phylopypruner min. length               : ${params.phypyPrunerMinLen}
phylopypruner min. support              : ${params.phypyPrunerMinSupport}
phylopypruner min. number of taxa       : ${params.minTax}
phylopypruner outgroup                  : ${params.phypyPrunerOutgroup}
raxml analysis name                     : ${params.raxMLName}
raxml model                             : ${params.raxMLModel}
raxml seed (-p)                         : ${params.raxMLP}
raxml asc-corr                          : ${params.raxMLascCorr}
raxml analysis name                     : ${params.raxMLName}
iqtree bootstrap                        : ${params.iqTrBS}
iqtree model                            : ${params.iqTrModel}
iqtree seed                             : ${params.iqTrSeed}
iqtree outgroup taxon                   : ${params.iqTrOutgroup}
astral minimum branch support           : ${params.astralCollapse}
email                                   : ${params.email}
"""

// this ensures that the parameters are listed before the pipeline exits
if (params.help) exit 1

// a common mistake is to pass resume as a parameter of the workflow instead of a nextflow flag
// --resume is a parameter while -resume is the actual flag that activates resume
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"



Channel
    .fromPath(params.proteins)
    .collect()
    .into {proteins_for_orthology; proteins_for_grab_script}

process proteinOrtho {
    label "proteinOrtho"

    publishDir './01_proteinOrtho'
    
    input:
    path(proteins_fasta) from proteins_for_orthology

    output:
    path("orthology_*.tsv") into proteinorthoOuttsv
    path("orthology_*blast-graph") into proteinorthoOutblastgraph
    path("orthology_*.info") into proteinorthoOutinfo
    path("orthology_*.proteinortho-graph") into proteinorthoOutprotgraph
    path("orthology_*.summary") into proteinorthoOutputsummary

    script:
    """
    proteinortho -cpus=$params.cpus -project=orthology_$params.proteinortho_evalue \
        -clean -sim=$params.proteinortho_sim -identity=$params.proteinortho_id $params.proteins
    """
}


process filteringProteins {
    label "protOrthPerl"
    debug = true

    publishDir './01_proteinOrtho'
    
    input:
    path(protein_tsv) from proteinorthoOuttsv

    output:
    file("orthology.reduced.tsv") into reduced_tsv_ch

    shell:
    """
    awk -v minTaxa=$params.minTax '/^#/ {header=\$0; print \$0}; \$1 > minTaxa' $protein_tsv > orthology.reduced.tsv
    """
    // Only worked with 'shell'. Did not work with 'script' 
}


process gettingProteins {
    label "protOrthPerl"
    debug = true
    publishDir './02_orthoGroups'
    
    input:
    path(reduced_tsv) from reduced_tsv_ch
    path(proteins_fasta) from proteins_for_grab_script

    output:
    path("*.fasta") into orthogroup_ch
    

    script:
    """
    proteinortho_grab_proteins.pl -t -exact $reduced_tsv $proteins_fasta
    """
}


process mafftAlignment {
    label "mafftAlign"

    publishDir './03_alignments'

    input:
    file (orthogroup_fasta) from orthogroup_ch.flatten()
    //Only worked with the 'flatten()' operator as opposed to 'collect()'

    output:
    path("*.fasta") into (aligned_orthogroup_ch, aligned_orthogroup2_ch)

    script:
    """
    mafft --anysymbol --auto --maxiterate 1000 $orthogroup_fasta > aligned_"$orthogroup_fasta"
    """
}


process aliscore {
    label "aliscore"
    errorStrategy 'ignore'

    time '1m'

    publishDir './04_aliscore'

    input:
    file(aligned_fasta) from aligned_orthogroup_ch

    output:
    path("*.txt") into aliscore_ch
    path("*.fasta") into alicut_ch


    script:
    """
    aliscore.pl -i $aligned_fasta
    if [ -s ${aligned_fasta}_List_random.txt]; then
        ALICUT_V2.31.pl -s
    else
        cat $aligned_fasta > ALICUT_$aligned_fasta
    fi
    """
}


process renameTrimmedAlignments {
    label "rename"

    publishDir './05_iqtreePhylogeny'

    input:
    path(alicut_trimmed_fasta) from alicut_ch

    output:
    path("*.fasta") into (renamed_align_ch, renamed_align_ch2)

    script:
    """
    rename 's/ALICUT_aligned_orthology.reduced.tsv.//g' $alicut_trimmed_fasta
    """
}


process iqtreePhylogeny_1 {
    errorStrategy 'ignore'
    label "iqtree_1"

    publishDir './05_iqtreePhylogeny'

    input:
    path (renamed_trimmed_fasta) from renamed_align_ch

    output:
    path("*") into iqtree_1_ch

    script:
    """
    iqtree2 -s $renamed_trimmed_fasta -st AA -m MFP -msub nuclear -B $params.iqtreeBS -T $params.cpus
    """
}


process phyloPruner {
    label "phylopypruner"
    
    publishDir './06_phyloPyPruner'

    input:
    path all from iqtree_1_ch.mix(renamed_align_ch2).collect()

    output:
    path("*") into (phyprune_ch, phyprune_ch_1, phyprune_ch_2)

    script:
    """
    phypyPrunerMinTax=\$((\$(ls -1 $params.proteins | wc -l)/2))
    mkdir input
    mv $all \$PWD/input/
    for fasta in ./input/*fasta; do gsed -i '/^>/ s/_/|/g' \$fasta; done
    for tree in ./input/*.treefile; do gsed -i 's/_/|/g' \$tree; done
    phylopypruner --dir ./input \
            --min-len $params.phypyPrunerMinLen \
            --trim-lb 5 \
            --min-support $params.phypyPrunerMinSupport \
            --prune MI --mask pdist \
            --min-taxa "\$phypyPrunerMinTax" \
            --outgroup $params.phypyPrunerOutgroup
    """
}


process raxmlPhylogeny {
    label "raxmlphylogeny"
    publishDir './07_raxmlPhyogeny'

    input:
    path phyprune_ch

    output:
    path("*") into raxml_ch

    script:
    """
    sed 's/AUTO/LG/g' "$phyprune_ch"/phylopypruner_output/partition_data.txt > "$phyprune_ch"/phylopypruner_output/partition_data_edit.txt
    cat "$phyprune_ch"/phylopypruner_output/partition_data_edit.txt

    raxmlHPC -s "$phyprune_ch"/phylopypruner_output/supermatrix.fas \
        -n $params.raxMLName \
        -m $params.raxMLModel \
        -o $params.phypyPrunerOutgroup \
        -p $params.raxMLP \
        -q "$phyprune_ch"/phylopypruner_output/partition_data_edit.txt \
        -T $params.cpus \
    """
}


process iqtreePhylogeny_2 {
    label "iqtree_2"

    publishDir "./08_iqtreePhylogeny"

    input:
    path phyprune_ch_1
    path raxml_ch

    output:
    path("*") into iqtree_2_ch

    script:
    """
    gsed -i 's/p78 = 4344/p79 = 4344/g' "$phyprune_ch_1"/phylopypruner_output/partition_data_edit.txt
    gsed -i 's/p5 = 16188/p5_2 = 16188/g' "$phyprune_ch_1"/phylopypruner_output/partition_data_edit.txt

    iqtree2 -nt $params.cpus -bb $params.iqTrBS \
        -bnni -safe -st AA \
        -s "$phyprune_ch_1"/phylopypruner_output/supermatrix.fas \
        -m $params.iqTrModel -fmax -seed $params.iqTrSeed \
        -ft "$raxml_ch"/RAxML_bestTree.raxmlOutput \
        -o $params.phypyPrunerOutgroup --prefix output \
        -q "$phyprune_ch_1"/phylopypruner_output/partition_data_edit.txt
    """ 
}


Channel
    .fromPath("$baseDir/06_phyloPyPruner/input/phylopypruner_output/output_alignments/*.fasta")
    .set {alignmentPath}

process iqtreePhylogeny_3 {
    label "iqtree_3"
    errorStrategy 'ignore'

    publishDir "./09_iqtreePhylogeny"

    input:
    path alignmentPath

    output:
    path("*") into iqtree_3_ch

    script:
    """
    iqtree2 -s $alignmentPath -st AA \
        -m MFP -msub nuclear \
        -B $params.iqtreeBS -nt $params.cpus
    """
}


Channel
    .fromPath("$baseDir/09_iqtreePhylogeny/*.treefile")
    .set {tree_ch}

process branchCollapse {
    label "branchcollapse"
    publishDir "./10_collapseBranches"

    input:
    path(treefile) from tree_ch

    output:
    path("*") into branchcollapse_ch

    script:
    """
    nw_ed $treefile 'i & b<=75' o > "$treefile"_BS75.fa.treefile 
    """
}


process astralPhylogeny {
    label "astral"
    debug = true
    publishDir "./11_astralPhylogeny"

    input:
    path(all) from branchcollapse_ch.collect()

    output:
    path("*") into astral_ch

    script:
    """
    cat $all | gsed 's/[|_].\\{5,7\\}:/:/g' >> multiple_trees.tre
    java -jar -Xmx6000M -D"java.library.path=lib/" ~/Astral/astral.5.7.3.jar \
        -i multiple_trees.tre -o consensus_tree.tre 
    """  
}