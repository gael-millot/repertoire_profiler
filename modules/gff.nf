// Makes a gff file from the coordinates info contained inside the tsv and pairs it with the fastas
process Gff {
    label 'r_ig_clustering'
    cache 'true'

    input:
    tuple path(fasta), path (other), val(tag), path(tsv) // parallelization expected (by clonal groups over align_clone_nb sequences) // warning: fasta is either nuc of aa, while other is the other one
    val align_seq
    val align_clone_nb
    path cute_file

    output:
    tuple path("*.gff"), val(tag), emit: gff_ch
    path "*_sequence_alignment_with_gaps_imgt_*.tempo", emit: imgt_gff_ch, optional: true
    path "Gff.log", emit: gff_log_ch
    path "warnings.txt", emit: gff_warn_ch

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
    "${align_clone_nb}" \
    "${cute_file}" \
    "${tag}" \
    "Gff.log"
    """
}