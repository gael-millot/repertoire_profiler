// This process takes the germline_alignment_d_mask column, removes the gaps and adds a column to the tsv
// It then adds to the tsv the translation in amino-acids of the germline_d_mask without gaps
process TranslateGermline {
    label 'r_ig_clustering'

    input:
    path add_germ_ch // parallelization expected (by clonal groups)

    output:
    path "*-trans_germ-pass.tsv", emit: translate_germ_ch
    path "*.log", emit: translate_germ_log_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    FILENAME=\$(basename -- ${add_germ_ch}) # recover a file name without path
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a translateGermline.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a translateGermline.log
    Rscript -e '
        file_name <- "${add_germ_ch}"
        df <- read.table(file_name, sep = "\\t", header = TRUE)
        germ_nuc <- df\$clonal_germline_sequence_no_gaps
        # Make sure all germline sequences are the same (which they are supposed to be since this is a single clonal group)
        if(any(germ_nuc != germ_nuc[1])){
            stop(paste0("\\n\\n================\\n\\nERROR IN TranslateGermline PROCESS.\\nTHE VALUES INSIDE THE Germline COLUMN SHOULD ALL BE THE SAME, BUT THEY ARE NOT.\\nHERE THEY ARE:\\n", paste0(germ_nuc, collapse = "\\n"),"\\n\\n================\\n\\n"), call. = FALSE)
        }
        length_no_gaps <- nchar(germ_nuc[1])
        if(length_no_gaps %% 3 != 0){
            cat(paste0("\\nWARNING:\nTHE clonal_germline_sequence_no_gaps COLUMN CONTAINS ", length_no_gaps, " CHARACTERS WHEN GAPS ARE REMOVED, WHICH IS NOT A MULTIPLE OF 3. \\n"), file = "translateGermline.log", append = TRUE)
        }
        germ_dna <- Biostrings::DNAString(germ_nuc[1])
        # Catch a warning in the log file if raised
        withCallingHandlers(
            expr = {
                germ_aa <- Biostrings::translate(germ_dna, if.fuzzy.codon="X")
            },
            warning = function(w) {
                cat("WARNING OF Biostrings::translate FUNCTION: ", conditionMessage(w), "\\n", file = "translateGermline.log", append = TRUE)
                invokeRestart("muffleWarning")
            }
        )
        tempo_name <- "clonal_germline_sequence_aa"
        df[[tempo_name]] <- toString(germ_aa)
        file_base <- tools::file_path_sans_ext(basename(file_name))
        new_file_name <- paste0(file_base, "-trans_germ-pass.tsv")
        # add controls
        df <- data.frame(df, clonal_germline_identical = df\$clonal_germline_sequence_no_gaps == df\$clonal_germline_alignment_igblast_airr, clonal_germline_aa_identical = df\$clonal_germline_sequence_aa == df\$clonal_germline_alignment_aa_igblast_airr)
        # end add controls
        write.table(df, file = paste0("./", new_file_name), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\\t")
    '|& tee -a translateGermline.log
    """

}

