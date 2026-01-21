
process  Mafft_align {
    label 'mafft'
    //publishDir path: "${out_path}/alignments/aa", mode: 'copy', pattern: "{*_aligned_aa.fasta}", overwrite: false
    //publishDir path: "${out_path}/alignments/nuc", mode: 'copy', pattern: "{*_aligned_nuc.fasta}", overwrite: false

    input:
    tuple path(fasta_nuc), path(fasta_aa), val(seq_kind) // parallelization expected (by clonal groups over align_clone_nb sequences)
    val align_mafft_all_options
    val align_mafft_clonal_options

    output:
    tuple path("*_aligned_nuc.fasta"), path("*_aligned_aa.fasta"), val(seq_kind), emit: aligned_all_ch
    path "mafft_align.log", emit: mafft_align_log_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    # --op N	Increase gap opening penalty (fewer gaps).
    # --ep N	Gap extension penalty (usually leave default).
    # --leavegappyregion	Don’t over-align ends or regions with many gaps.
    # --keeplength	Preserve sequence end lengths.
    # --localpair --maxiterate N	Use accurate local iterative refinement.
    if [[ "${seq_kind}" == "ALL" ]] ; then
        mafft ${align_mafft_all_options} ${fasta_nuc} | awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' > ${fasta_nuc.baseName}_aligned_nuc_tempo.fasta
        awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print toupper(\$0) ; next}}END{print "\\n"}' ${fasta_nuc.baseName}_aligned_nuc_tempo.fasta > ${fasta_nuc.baseName}_aligned_nuc.fasta # remove \\n in seq
        mafft ${align_mafft_all_options} ${fasta_aa} | awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' > ${fasta_aa.baseName}_aligned_aa_tempo.fasta
        awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' ${fasta_aa.baseName}_aligned_aa_tempo.fasta > ${fasta_aa.baseName}_aligned_aa.fasta # remove \\n in seq
    elif [[ "${seq_kind}" == "CLONE" ]] ; then
        mafft ${align_mafft_clonal_options} ${fasta_nuc} | awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' > ${fasta_nuc.baseName}_aligned_nuc_tempo.fasta
        awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print toupper(\$0) ; next}}END{print "\\n"}' ${fasta_nuc.baseName}_aligned_nuc_tempo.fasta > ${fasta_nuc.baseName}_aligned_nuc.fasta # remove \\n in seq
        mafft ${align_mafft_clonal_options} ${fasta_aa} | awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' > ${fasta_aa.baseName}_aligned_aa_tempo.fasta
        awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' ${fasta_aa.baseName}_aligned_aa_tempo.fasta > ${fasta_aa.baseName}_aligned_aa.fasta # remove \\n in seq
    fi
    echo "" > mafft_align.log
    """
}

