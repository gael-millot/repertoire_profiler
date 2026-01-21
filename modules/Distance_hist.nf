process Distance_hist {
    label 'immcantation'
    publishDir path: "${out_path}/figures/png", mode: 'copy', pattern: "{*.png}", overwrite: false
    publishDir path: "${out_path}/figures/svg", mode: 'copy', pattern: "{*.svg}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{distance_hist.log}", overwrite: false
    cache 'true'

    input:
    path distToNearest_ch // no parallelization
    path cute_file
    val clone_model
    val clone_normalize
    val clone_distance

    output:
    path "seq*.pdf", emit: histogram_pdf_ch
    path "*.png", emit: distance_hist_ch // png plot (but sometimes empty) sustematically returned
    path "*.svg"
    path "distance_hist.log"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    histogram.R \
"${distToNearest_ch}" \
"${clone_model}" \
"${clone_normalize}" \
"${clone_distance}" \
"${cute_file}" \
"distance_hist.log"
    """
}