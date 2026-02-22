
process Tree_nuc {
    publishDir path: "${out_path}/phylo/nuc", mode: 'copy', pattern: "{*_aligned_nuc.fasta.treefile}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{tree_nuc.log}", overwrite: false

    label 'iqtree'

    cache 'true'

    input:
    path fasta_nuc_alignments  // parallelization expected (by clonal groups over align_clone_nb sequences)
    path phylo_tree_model_file

    output:
    path "*_aligned_nuc.fasta.treefile", emit: tree_nuc_ch, optional: true
    path "tree_nuc.log", emit: tree_log_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    if (( \$(wc -l < ${fasta_nuc_alignments}) / 2 >= 3 )) ; then
        # see https://iqtree.github.io/doc/Command-Reference
        iqtree -nt AUTO -s ${fasta_nuc_alignments} -m ${phylo_tree_model_file}+GTR+I+R6 --seed 123456789 |& tee -a tree_nuc.log
    else
        echo -e "\\n\\nNUMBER OF NUC ALIGNED SEQUENCES < 3: NO TREE RETURNED.\\n\\n" |& tee -a tree_nuc.log
    fi
    """
}
