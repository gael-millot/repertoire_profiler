process PrintAlignment{
    label 'goalign'

    input:
    path align_ch // parallelization expected (by clonal groups over clone_nb_seq sequences)

    output:
    path "*.html", emit : alignment_html

    script:
    """
    goalign draw biojs -i ${align_ch} -o ${align_ch.baseName}.html
    """
}