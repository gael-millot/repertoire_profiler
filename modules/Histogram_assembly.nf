

process Histogram_assembly {
    label 'r_ig_clustering'
    publishDir path: "${out_path}/pdf", mode: 'copy', pattern: "{seq_distance.pdf}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{histogram_assembly.log}", overwrite: false
    cache 'true'

    input:
    path histogram_pdf_ch

    output:
    path "seq_distance.pdf"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    Rscript -e '
    # assignation to prevent a returned element
        tempo <- qpdf::pdf_combine(input = list.files(path = ".", pattern = "^seq.*.pdf\$"), output = "./seq_distance.pdf")
    ' |& tee -a histogram_assembly.log
    """
}
