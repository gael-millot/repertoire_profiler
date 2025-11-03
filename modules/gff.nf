// Makes a gff file from the coordinates info contained inside the tsv and pairs it with the fastas
process Gff {
    label 'r_ig_clustering'
    cache 'true'

    input:
    tuple path(all_files_ch), val(seq_kind) // parallelization expected (by clonal groups over align_clone_nb sequences)
    val align_seq
    val align_clone_nb
    path cute_file

    output:
    path "*.gff", emit: gff_ch, optional: true
    path "Gff.log", emit: gff_log_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    FILENAME=\$(basename -- ${all_files_ch}) # recover a file name without path
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a Gff.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a Gff.log
    Gff.R \
    "${all_files_ch}" \
    "sequence_id" \
    "${align_seq}" \
    "${align_clone_nb}" \
    "${cute_file}" \
    "${seq_kind}" \
    "Gff.log"
    """
}