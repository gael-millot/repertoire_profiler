
process Tree {
    publishDir path: "${out_path}/phylo/aa", mode: 'copy', pattern: "{*_aligned_aa.fasta.treefile}", overwrite: false
    publishDir path: "${out_path}/phylo/nuc", mode: 'copy', pattern: "{*_aligned_nuc.fasta.treefile}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{tree.log}", overwrite: false

    label 'iqtree'

    cache 'true'

    input:
    tuple path(fasta_nuc_alignments), path(fasta_aa_alignments)  // parallelization expected (by clonal groups over align_clone_nb sequences)
    path phylo_tree_model_file

    output:
    path "*_aligned_aa.fasta.treefile", emit: tree_aa_ch, optional: true
    path "*_aligned_nuc.fasta.treefile", emit: tree_nuc_ch, optional: true
    path "tree.log", emit: tree_log_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    if (( \$(wc -l < ${fasta_aa_alignments}) / 2 >= 3 )) ; then
        # see https://iqtree.github.io/doc/Command-Reference
        iqtree -nt AUTO -s ${fasta_aa_alignments} -m ${phylo_tree_model_file}+I+R6 --seed 123456789 |& tee -a tree.log
    else
        echo -e "\\n\\nNUMBER OF AA ALIGNED SEQUENCES < 3: NO TREE RETURNED.\\n\\n" |& tee -a tree.log
    fi
    if (( \$(wc -l < ${fasta_nuc_alignments}) / 2 >= 3 )) ; then
        iqtree -nt AUTO -s ${fasta_nuc_alignments} -m ${phylo_tree_model_file}+GTR+I+R6 --seed 123456789 -T AUTO |& tee -a tree.log
    else
        echo -e "\\n\\nNUMBER OF NUC ALIGNED SEQUENCES < 3: NO TREE RETURNED.\\n\\n" |& tee -a tree.log
    fi
    """
}
