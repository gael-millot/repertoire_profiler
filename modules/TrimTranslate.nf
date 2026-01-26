// Trim the 5' end of the sequence if required (to remove the leader peptide) and Translate the nucleotidic sequences into amino acid sequences
process TrimTranslate {
    label 'seqkit'
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "wanted_nuc/trimmed/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "wanted_nuc/removed/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "wanted_nuc/query/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "wanted_nuc/align/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "wanted_nuc/align_with_gaps/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "wanted_aa/trimmed/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "wanted_aa/igblast/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "wanted_aa/query/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "wanted_aa/align/*.fasta", overwrite: false
    cache 'true'

    input:
    tuple path(checked_tsv_ch), val(nlines) // parallelization
    val nb_wanted

    when:
    nb_wanted > 0
    nlines > 1

    output:
    path "trimtranslate.tsv", emit: trimtranslate_ch // productive file with column sequence_alignment_aa added
    path "wanted_nuc/trimmed/*.*"
    path "wanted_nuc/removed/*.*"
    path "wanted_nuc/query/*.*"
    path "wanted_nuc/align/*.*"
    path "wanted_nuc/align_with_gaps/*.*"
    path "wanted_aa/trimmed/*.*"
    path "wanted_aa/igblast/*.*"
    path "wanted_aa/query/*.*"
    path "wanted_aa/align/*.*"
    path "trimtranslate.log", emit: trimtranslate_log_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    trimtranslate.sh ${checked_tsv_ch} # |& tee -a trimtranslate.log not used because in trimtranslate.sh
    """
}