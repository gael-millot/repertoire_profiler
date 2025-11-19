// Trim the 5' end of the sequence if required (to remove the leader peptide) and Translate the nucleotidic sequences into amino acid sequences
process TrimTranslate {
    label 'seqkit'
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "productive_nuc/trimmed/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "productive_nuc/removed/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "productive_nuc/query/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "productive_nuc/align/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "productive_nuc/align_with_gaps/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "productive_aa/trimmed/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "productive_aa/igblast/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "productive_aa/query/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "productive_aa/align/*.fasta", overwrite: false
    cache 'true'

    input:
    path select_ch // parallelization expected

    output:
    path "trimtranslate.tsv", emit: trimtranslate_ch // productive file with column sequence_alignment_aa added
    path "productive_nuc/trimmed/*.*"
    path "productive_nuc/removed/*.*"
    path "productive_nuc/query/*.*"
    path "productive_nuc/align/*.*"
    path "productive_nuc/align_with_gaps/*.*"
    path "productive_aa/trimmed/*.*"
    path "productive_aa/igblast/*.*"
    path "productive_aa/query/*.*"
    path "productive_aa/align/*.*"
    path "trimtranslate.log", emit: trimtranslate_log_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    trimtranslate.sh ${select_ch} # |& tee -a trimtranslate.log not used because in trimtranslate.sh
    """
}