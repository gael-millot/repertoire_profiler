// Makes a gff file from the coordinates info contained inside the tsv and pairs it with the fastas
process Gff_imgt {
    label 'r_ig_clustering'
    publishDir path: "${out_path}/alignments/nuc/imgt", mode: 'copy', pattern: "{*_imgt_nuc.gff}", overwrite: false
    publishDir path: "${out_path}/alignments/aa/imgt", mode: 'copy', pattern: "{*_imgt_aa.gff}", overwrite: false 
    cache 'true'

    input:
    tuple path(fasta), path(aa), val(tag) // tag not required here 
    path tsv  // no parallelization
    path cute_file

    output:
    path "*_imgt_nuc.gff"
    path "*_imgt_aa.gff"
    path "Gff_imgt.log", emit: gff_log_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    FILENAME=\$(basename -- ${fasta}) # recover a file name without path
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a Gff.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a Gff.log
    Gff_imgt.R \
    "${fasta}" \
    "${tsv}" \
    "sequence_id" \
    "${cute_file}" \
    "Gff_imgt.log"
    """
}