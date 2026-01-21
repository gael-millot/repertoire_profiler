
process Donut_assembly {
    label 'r_ig_clustering'
    publishDir path: "${out_path}/pdf", mode: 'copy', pattern: "{donuts.pdf}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{donut_assembly.log}", overwrite: false
    cache 'true'

    input:
    path donut_pdf_ch2

    output:
    path "donuts.pdf", emit: donut_assembly_ch
    path "donut_assembly.log"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    Rscript -e '
        files <- list.files(path = ".", pattern = ".pdf\$")

        sorted_files <- files[order(
            !grepl("allele", files),  # Alleles will be displayed before genes
            grepl("gene", files),     
            !grepl("vj_", files)     # Amongst alleles and genes, vj will be displayed before c
        )]

        qpdf::pdf_combine(input = sorted_files, output = "./donuts.pdf")
    ' |& tee -a donut_assembly.log
    """
}
