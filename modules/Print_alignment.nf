process PrintAlignment{
    label 'goalign'
    cache 'true'
    publishDir path: "${out_path}/alignments/nuc/imgt", mode: 'copy', pattern: "*.html", overwrite: false
    publishDir path: "${out_path}/alignments/aa/imgt", mode: 'copy', pattern: "*.html", overwrite: false

    input:
    tuple path(align_ch), path(other), val(tag)  // parallelization expected (by clonal groups over clone_nb_seq sequences)


    output:
    tuple path("*.html"), val(tag), emit : alignment_html, optional: true

    script:
    """
    goalign draw biojs -i ${align_ch} -o ${align_ch.baseName}.html
    """
}