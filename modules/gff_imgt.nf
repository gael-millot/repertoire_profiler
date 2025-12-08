// Makes a gff file from the coordinates info contained inside the tsv and pairs it with the fastas
process Gff_imgt {
    label 'r_ig_clustering'
    publishDir path: "${out_path}/alignments/nuc/imgt", mode: 'copy', pattern: "{*_imgt_nuc_biojs.gff}", overwrite: false
    publishDir path: "${out_path}/alignments/aa/imgt", mode: 'copy', pattern: "{*_imgt_aa_biojs.gff}", overwrite: false 
    //publishDir path: "${out_path}/alignments/nuc/imgt", mode: 'copy', pattern: "{*_imgt_nuc_jalview2.gff}", overwrite: false
    //publishDir path: "${out_path}/alignments/aa/imgt", mode: 'copy', pattern: "{*_imgt_aa_jalview2.gff}", overwrite: false 
    cache 'true'

    input:
    tuple path(fasta), path(aa), val(tag) // even if tag not required here 
    val align_seq
    path tsv  // no parallelization
    path cute_file

    output:
    path "*_imgt_nuc_biojs.gff"
    path "*_imgt_aa_biojs.gff"
    path "*_imgt_nuc_jalview.gff", emit: nuc_imgt_gff_ch, optional: true
    path "*_imgt_aa_jalview.gff", emit: aa_imgt_gff_ch, optional: true
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
    "${align_seq}" \
    "${cute_file}" \
    "Gff_imgt.log"
    """
}