
process  Mafft_align {
    label 'mafft'

    cache 'false'
    //publishDir path: "${out_path}/alignments/aa", mode: 'copy', pattern: "{*_aligned_aa.fasta}", overwrite: false
    //publishDir path: "${out_path}/alignments/nuc", mode: 'copy', pattern: "{*_aligned_nuc.fasta}", overwrite: false

    input:
    tuple path(fasta_nuc), path(fasta_aa), val(seq_kind) // parallelization expected (by clonal groups over align_clone_nb sequences)
    val align_mafft_all_options
    val align_mafft_clonal_options

    output:
    tuple path("*_aligned_nuc.fasta"), path("*_aligned_aa.fasta"), val(seq_kind), emit: aligned_all_ch
    path "mafft_align.log", emit: mafft_align_log_ch
    path "warnings.txt", emit: mafft_align_warn_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    # --op N	Increase gap opening penalty (fewer gaps).
    # --ep N	Gap extension penalty (usually leave default).
    # --leavegappyregion	Don’t over-align ends or regions with many gaps.
    # --keeplength	Preserve sequence end lengths.
    # --localpair --maxiterate N	Use accurate local iterative refinement.
    # file truncation
    # Indeed Mafft has problem with more than 1355 sequences. Thus I fixed a limit at 1350 * 2 = 2700 fata lines (1350 sequences)
    # nuc_seq_nb=\$((\$(wc -l ${fasta_nuc} | cut -f1 -d' ') / 2))
    nuc_seq_nb=\$(wc -l ${fasta_nuc} | cut -f1 -d' ')
    if (("\${nuc_seq_nb}" <= "2700")) ; then 
        cp ${fasta_nuc} tempo_nuc.fasta
        echo "" >> warnings.txt
    else
        cat ${fasta_nuc} | head -2700 > tempo_nuc.fasta
        echo -e "\\n\\nWARNING:\\nMORE THAN 1,350 SEQUENCES (\$((nuc_seq_nb / 2))) DETECTED IN ${fasta_nuc}.\\n FILE IS TRUNCATED TO 1,350 SEQUENCES FOR MAFFT ALIGNMENT.\\n\\n" >> warnings.txt
    fi
    aa_seq_nb=\$(wc -l ${fasta_aa} | cut -f1 -d' ')
    if (("\${aa_seq_nb}" <= "2700")) ; then 
        cp ${fasta_aa} tempo_aa.fasta
        echo "" >> warnings.txt
    else
        cat ${fasta_aa} | head -2700 > tempo_aa.fasta
        echo -e "\\n\\nWARNING:\\nMORE THAN 1,350 SEQUENCES (\$((aa_seq_nb / 2))) DETECTED IN ${fasta_aa}.\\n FILE IS TRUNCATED TO 1,350 SEQUENCES FOR MAFFT ALIGNMENT.\\n\\n" >> warnings.txt
    fi
    # end file truncation
    if [[ "${seq_kind}" == "ALL" ]] ; then
        mafft ${align_mafft_all_options} tempo_nuc.fasta | awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' > ${fasta_nuc.baseName}_aligned_nuc_tempo.fasta
        awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print toupper(\$0) ; next}}END{print "\\n"}' ${fasta_nuc.baseName}_aligned_nuc_tempo.fasta > ${fasta_nuc.baseName}_aligned_nuc.fasta # remove \\n in seq
        mafft ${align_mafft_all_options} tempo_aa.fasta | awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' > ${fasta_aa.baseName}_aligned_aa_tempo.fasta
        awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' ${fasta_aa.baseName}_aligned_aa_tempo.fasta > ${fasta_aa.baseName}_aligned_aa.fasta # remove \\n in seq
    elif [[ "${seq_kind}" == "CLONE" ]] ; then
        mafft ${align_mafft_clonal_options} tempo_nuc.fasta | awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' > ${fasta_nuc.baseName}_aligned_nuc_tempo.fasta
        awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print toupper(\$0) ; next}}END{print "\\n"}' ${fasta_nuc.baseName}_aligned_nuc_tempo.fasta > ${fasta_nuc.baseName}_aligned_nuc.fasta # remove \\n in seq
        mafft ${align_mafft_clonal_options} tempo_aa.fasta | awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' > ${fasta_aa.baseName}_aligned_aa_tempo.fasta
        awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' ${fasta_aa.baseName}_aligned_aa_tempo.fasta > ${fasta_aa.baseName}_aligned_aa.fasta # remove \\n in seq
    fi
    echo "" > mafft_align.log
    """
}

