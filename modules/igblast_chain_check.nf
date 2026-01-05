
process Igblast_chain_check {
    label 'r_ig_clustering'
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{igblast_chain_check.log}", overwrite: false
    cache 'true'

    input:
    path data_assembly_ch // no parallelization
    path igblast_data_check_ch
    path cute_file

    output:
    path "productive_seq.tsv", emit: igblast_chain_check_ch
    path "igblast_chain_check.log"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    igblast_chain_check.R \
"${data_assembly_ch}" \
"${igblast_data_check_ch}" \
"${cute_file}" \
"igblast_chain_check.log"
    """
}