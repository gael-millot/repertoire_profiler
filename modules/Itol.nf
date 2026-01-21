
process Itol{
    label 'gotree'

    input:
    path tree
    path itol_files
    val phylo_tree_itolkey

    output:
    path "*_aligned_aa.fasta_itol_url.txt", emit: itol_aa_ch, optional: true
    path "*_aligned_nuc.fasta_itol_url.txt", emit: itol_nuc_ch, optional: true

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    # Check if itol_files is made of iTOL.empty # -s means exists and non-empty
    if [ ! -s "iTOL.empty" ]; then 
        gotree upload itol --project gotree_uploads --user-id $phylo_tree_itolkey -i $tree > ${tree.baseName}_itol_url.txt 2>&1
    else
        # add metadata
        gotree upload itol --project gotree_uploads --user-id $phylo_tree_itolkey -i $tree $itol_files > ${tree.baseName}_itol_url.txt 2>&1 
    fi
    """
}
