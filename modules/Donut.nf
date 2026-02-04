// Creates donut plots for the constant and variable region ; grouped by same allele or genes
// Inputs:
//      - kind: kind of donut plot to display. can be "all", "annotated" or "clone"
//      - data: tsv file containing the sequences to plot 
//      - col: columns in the data file to take into account when plotting the donut (can be c_call if the donut should be grouped by same constant region alleles for instance)
//      - donut_*: parameters for the display of the donut. These are entered in the nextflow.config file
//      - cute_file: file containing R functions used in the donut.R script
//      - igblast_variable_ref_files & igblast_constant_ref_files: only used to be certain that all studied loci are displayed in the donut's title (if we included IGL and IGK in the study but only one was found, both loci still need to be in the title)
process Donut {
    label 'r_ext'
    //publishDir path: "${out_path}", mode: 'copy', pattern: "{*.tsv}", overwrite: false
    //publishDir path: "${out_path}", mode: 'copy', pattern: "{*.pdf}", overwrite: false
    publishDir path: "${out_path}/figures/png", mode: 'copy', pattern: "{*.png}", overwrite: false
    publishDir path: "${out_path}/figures/svg", mode: 'copy', pattern: "{*.svg}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{*.log}", overwrite: false
    cache 'true'

    input:
    tuple val(kind), path(data), val(col) // 12 parallelization expected
    val donut_palette
    val donut_hole_size
    val donut_hole_text
    val donut_hole_text_size
    val donut_border_color
    val donut_border_size
    val donut_annotation_distance
    val donut_annotation_size
    val donut_annotation_force
    val donut_annotation_force_pull
    val donut_legend_width
    val donut_legend_text_size
    val donut_legend_box_size
    val donut_legend_box_space
    val donut_legend_limit
    path cute_file
    path allele_names_tsv_all_ch

    output:
    path "*.tsv", emit: donut_tsv_ch, optional: true
    path "*.pdf", emit: donut_pdf_ch, optional: true
    path "*.png", emit: donuts_png_ch
    path "*.svg"
    path "*.log"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    FILENAME=\$(basename -- ${data}) # recover a file name without path
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\nKIND: ${kind}\\nCOL: ${col}\\n\\n################################\\n\\n" |& tee -a ${kind}_donut.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a ${kind}_donut.log
    donut.R \
"${data}" \
"${kind}" \
"${col}" \
"${donut_palette}" \
"${donut_hole_size}" \
"${donut_hole_text}" \
"${donut_hole_text_size}" \
"${donut_border_color}" \
"${donut_border_size}" \
"${donut_annotation_distance}" \
"${donut_annotation_size}" \
"${donut_annotation_force}" \
"${donut_annotation_force_pull}" \
"${donut_legend_width}" \
"${donut_legend_text_size}" \
"${donut_legend_box_size}" \
"${donut_legend_box_space}" \
"${donut_legend_limit}" \
"${cute_file}" \
"${allele_names_tsv_all_ch}" \
"${kind}_donut.log"
    """
}
