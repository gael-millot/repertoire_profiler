// This process takes the germline_alignment_d_mask column, removes the gaps and adds a column to the tsv
// It then adds to the tsv the translation in amino-acids of the germline_d_mask without gaps
process TranslateGermline {
    label 'r_ig_clustering'

    input:
    path add_germ_ch // parallelization expected (by clonal groups)
    path cute_file

    output:
    path "*-trans_germ-pass.tsv", emit: translate_germ_ch
    path "clonal_germline_sequence_no_gaps_problems.tsv", emit: translate_problem_ch
    path "*.log", emit: translate_germ_log_ch
    path "warnings.txt", emit: translate_germ_warn_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    FILENAME=\$(basename -- ${add_germ_ch}) # recover a file name without path
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a translateGermline.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a translateGermline.log
    TranslateGermline.R \
    "${add_germ_ch}" \
    "${cute_file}" \
    "translateGermline.log"
    """
}

