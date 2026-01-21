

process Meta2Itol  {
    label 'tabletoitol'

    input:
    path meta_file // no parallelization
    val meta_seq_names

    output:
    path "iTOL*", emit: itol_out_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    if [ "${meta_file}" == "NULL" ]; then
        echo "No metadata provided; skipping iTOL annotation" > skip.iTOL.txt
        touch iTOL.empty # create the empty file iTOL.empty
        exit 0
    else
        table2itol.R -i ${meta_seq_names} ${meta_file}
        # create three files: iTOL_colorstrip-Isotype.txt, iTOL_colorstrip-Name.txt, iTOL_gradient-KD.txt
    fi
    """
}
