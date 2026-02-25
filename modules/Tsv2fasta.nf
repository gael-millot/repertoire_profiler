// Makes a fasta file with several sequences based on the sequences present in the tsv input in BOTH nucleotidic format and amino-acidic format
// The tsv input is expected to already only contain sequences of one same clonal group. For each clonal group, both nuc and aa fastas will be emitted paired up.
process Tsv2fasta {
    label 'r_ig_clustering'
    cache 'true'
    // pushished files are dealt after the process but these work:
    // publishDir path: "${out_path}/fasta/for_alignment_nuc", mode: 'copy', pattern: "for_alignment_nuc/*.fasta", overwrite: false
    // publishDir path: "${out_path}/fasta/for_alignment_aa", mode: 'copy', pattern: "for_alignment_aa/*.fasta", overwrite: false
    // publishDir path: "${out_path}/alignments/nuc/imgt", mode: 'copy', pattern: "sequence_alignment_with_gaps_imgt_nuc.fasta", overwrite: false
    // publishDir path: "${out_path}/alignments/aa/imgt", mode: 'copy', pattern: "sequence_alignment_with_gaps_imgt_aa.fasta", overwrite: false
    

    input:
    tuple path(all_files_ch), val(seq_kind) // parallelization expected (by clonal groups over align_clone_nb sequences)
    val align_seq
    val clone_germline_kind
    val align_clone_nb
    path cute_file

    output:
    tuple path("for_alignment_nuc/*.fasta"), path("for_alignment_aa/*.fasta"), val(seq_kind), emit: fasta_align_ch, optional: true // // already aligned fasta file with seq_kind == "IMGT"
    tuple path("sequence_alignment_with_gaps_imgt_nuc.fasta"), path("sequence_alignment_with_gaps_imgt_aa.fasta"), val(seq_kind), emit: fasta_align_imgt_ch, optional: true // // already aligned fasta file with seq_kind == "IMGT"
    path "Tsv2fasta.log", emit: tsv2fasta_log_ch
    path "warnings.txt", emit: warning_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    FILENAME=\$(basename -- ${all_files_ch}) # recover a file name without path
    echo -e "\\n\\n################################\\n\\nFILE \${FILENAME}\\n\\nKIND OF SEQ  ${seq_kind}\\n\\n################################\\n\\n" |& tee -a Tsv2fasta.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a Tsv2fasta.log
    Tsv2fasta.R \
    "${all_files_ch}" \
    "sequence_id" \
    "${align_seq}" \
    "${clone_germline_kind}" \
    "${align_clone_nb}" \
    "${cute_file}" \
    "${seq_kind}" \
    "Tsv2fasta.log"
    """
}
