
// This process aligns the nucleotidic fasta files by transfering amino-acidic alignments onto the nucleotidic ones
process Abalign_align_nuc {
    label 'goalign'
    //publishDir path: "${out_path}/alignments/nuc", mode: 'copy', pattern: "{*_aligned_nuc.fasta}", overwrite: false

    input:
    tuple path(fasta_nuc), path(aligned_aa), val(seq_kind) // parallelization expected (by clonal groups over align_clone_nb sequences)

    output:
    tuple path("*_aligned_nuc.fasta"), path(aligned_aa), val(seq_kind), emit: aligned_all_ch
    path "abalign_align_nuc.log", emit: abalign_align_nuc_log_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    FILENAME=\$(basename -- ${fasta_nuc}) # recover a file name without path
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a abalign_align_nuc.log

    goalign codonalign -i ${aligned_aa} -f ${fasta_nuc} -o ${fasta_nuc.baseName}_aligned_nuc_tempo.fasta |& tee -a abalign_align_nuc.log
    awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' ${fasta_nuc.baseName}_aligned_nuc_tempo.fasta > ${fasta_nuc.baseName}_aligned_nuc.fasta # remove \\n in seq
    """
}

