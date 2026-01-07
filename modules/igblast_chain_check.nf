
process Igblast_chain_check {
    label 'r_ig_clustering'
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{igblast_chain_check.log}", overwrite: false
    cache 'true'

    input:
    path data_assembly_ch // no parallelization
    path allele_names_tsv_all_ch // to have all the files in the work dir
    val igblast_v_ref_files 
    val igblast_d_ref_files 
    val igblast_j_ref_files 
    val igblast_constant_ref_files 
    path cute_file

    output:
    path "selected.tsv", emit: checked_tsv_ch, optional: true
    path "not_selected.tsv", emit: not_checked_tsv_ch, optional: true
    path "igblast_chain_check.log", emit: igblast_chain_check_log

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    # igblast_data_check_ch not required here
    igblast_chain_check.R \
"${data_assembly_ch}" \
"${igblast_v_ref_files}" \
"${igblast_d_ref_files}" \
"${igblast_j_ref_files}" \
"${igblast_constant_ref_files}" \
"${cute_file}" \
"igblast_chain_check.log"
    """
}