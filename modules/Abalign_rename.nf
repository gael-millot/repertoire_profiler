
// Abalign puts fasta headers in all caps. next script is meant to put those headers back to how they originally were
process Abalign_rename {
    label 'r_ext'
    //publishDir path: "${out_path}/alignments/aa", mode: 'copy', pattern: "{*_aligned_aa.fasta}", overwrite: false
    
    input:
    tuple path(fasta_nuc), path(fasta_aa), path(fasta_aa_align), val(seq_kind) // parallelization expected (by clonal groups over align_clone_nb sequences)
    
    output:
    tuple path(fasta_nuc), path("*_aligned_aa.fasta"), val(seq_kind), emit: renamed_aligned_aa_ch
    path "*_failed_abalign_align.tsv", emit: failed_abalign_align_ch
    path "abalign_rename.log", emit: abalign_rename_log_ch
    path "warnings.txt", emit: abalign_rename_warn_ch
    
    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    Rscript -e '
        options(warning.length = 8170)
        warn <- ""
        args = commandArgs(trailingOnly=TRUE)
        df_ini <- readLines(args[1])
        df_align <- readLines(args[2])
        ini_header_pos <- which(grepl(x = df_ini, pattern = "^>"))
        align_header_pos <- which(grepl(x = df_align, pattern = "^>"))
        if( ! all((ini_header_pos + 1) %% 2 == 0)){
            stop(paste0("\\n\\n========\\n\\nERROR IN Abalign_rename PROCESS\\n\\nHEADER POSITIONS CANNOT BE EVEN NUMBERS.\\nHERE THEY ARE:\\n", ini_header_pos, "\\n\\n========\\n\\n"), call. = FALSE)
        }
        if( ! all((align_header_pos + 1) %% 2 == 0)){
            stop(paste0("\\n\\n========\\n\\nERROR IN Abalign_rename PROCESS\\n\\nHEADER POSITIONS CANNOT BE EVEN NUMBERS.\\nHERE THEY ARE:\\n", align_header_pos, "\\n\\n========\\n\\n"), call. = FALSE)
        }
        failed_abalign_align <- c("fasta_name\\theader")
        df_ini_low <- tolower(df_ini)
        df_align_low <- tolower(df_align)
        tempo_log <-  ! df_ini_low[ini_header_pos] %in% df_align_low[align_header_pos]
        if(any(tempo_log)){
            tempo_txt <- paste0(paste0("${fasta_aa.baseName}", "\\t", sub(x = df_ini[ini_header_pos][tempo_log], pattern = "^>", replacement = "")), collapse = "\\n")
            tempo_warn <- paste0("\\n\\nWARNING:\nIN THE Abalign_rename PROCESS, ALIGNMENT FAILED FOR ", tempo_txt, "\\n\\n")
            warn <- base::paste0(base::ifelse(test = warn == "", yes = tempo_warn, no = base::paste0(warn, "\n\n", tempo_warn, collapse = NULL, recycle0 = FALSE)), collapse = NULL, recycle0 = FALSE)
            print(tempo_warn)
            failed_abalign_align <- paste0("fasta_name\\theader\\n", paste0(tempo_txt, collapse = "\\n"))
        }
        writeLines(failed_abalign_align, con = "${fasta_aa.baseName}_failed_abalign_align.tsv")
        for(i2 in align_header_pos){
            tempo_header_align_low <- tolower(df_align[i2])
            tempo_log <- df_ini_low %in% tempo_header_align_low
            if(sum(tempo_log) != 1){
                stop(paste0("\\n\\n========\\n\\nERROR IN Abalign_rename PROCESS\\n\\ntempo_header_align_low MUST EXIST IN df_ini_low.\\ntempo_header_align_low:\\n", tempo_header_align_low, "\\ndf_ini_low:\\n", df_ini_low, "\\n\\n========\\n\\n"), call. = FALSE)
            }
            df_align[i2] <- df_ini[tempo_log]
        }
        writeLines(warn, con = "warnings.txt")
        writeLines(df_align, con = "${fasta_aa.baseName}_aligned_aa_tempo.fasta")
    ' ${fasta_aa} ${fasta_aa.baseName}_align_aa.fasta |& tee -a abalign_rename.log
    awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' ${fasta_aa.baseName}_aligned_aa_tempo.fasta > ${fasta_aa.baseName}_aligned_aa.fasta # remove \\n in seq
    """
}
