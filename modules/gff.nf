// Makes a gff file from the coordinates info contained inside the tsv and pairs it with the fastas
process Gff {
    label 'r_ig_clustering'
    cache 'true'

    input:
    tuple path(fasta), val(tag), path(tsv) // parallelization expected (by clonal groups over align_clone_nb sequences)
    val align_seq
    val clone_germline_kind
    val align_clone_nb
    path cute_file

    output:
    path "*.gff", emit: gff_ch, optional: true
    path "Gff.log", emit: gff_log_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    FILENAME=\$(basename -- ${fasta}) # recover a file name without path
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a Gff.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a Gff.log
    Gff.R \
    "${fasta}" \
    "${tsv}" \
    "sequence_id" \
    "${align_seq}" \
    "${clone_germline_kind}" \
    "${align_clone_nb}" \
    "${cute_file}" \
    "${tag}" \
    "Gff.log"
    """
}