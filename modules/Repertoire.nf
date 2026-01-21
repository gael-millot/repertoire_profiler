
process Repertoire {
    label 'r_ig_clustering'
    publishDir path: "${out_path}/pdf", mode: 'copy', pattern: "{*_repertoire.pdf}", overwrite: false
    publishDir path: "${out_path}/figures/png", mode: 'copy', pattern: "{*.png}", overwrite: false
    publishDir path: "${out_path}/figures/svg", mode: 'copy', pattern: "{*.svg}", overwrite: false
    publishDir path: "${out_path}/repertoires", mode: 'copy', pattern: "{rep_*.tsv}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{repertoire.log}", overwrite: false
    cache 'true'

    input:
    path seq_name_replacement_ch2 // no parallelization
    path allele_names_tsv_all_ch
    path cute_file

    output:
    path "*_repertoire.pdf"
    path "*.svg"
    path "*.png", emit: repertoire_png_ch
    path "rep_*.tsv"
    path "repertoire.log"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    repertoire.R \
"${seq_name_replacement_ch2}" \
"${allele_names_tsv_all_ch}" \
"${cute_file}" \
"repertoire.log"
    """
}


