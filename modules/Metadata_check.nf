process Metadata_check { // cannot be in germ_tree_vizu because I have to use the clone_assigned_seq.tsv file for the check
    label 'immcantation'
    publishDir "${out_path}/reports", mode: 'copy', pattern: "{metadata_check.log}", overwrite: false
    cache 'true'

    input:
    path data_assembly_ch // no parallelization
    path meta_file
    val meta_seq_names
    val meta_name_replacement
    val meta_legend

    output:
    path "metadata_check.log"
    path "warnings.txt", emit: metadata_check_warn_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    if [[ "${meta_file}" != "NULL" ]] ; then
        Rscript -e '
            options(warning.length = 8170)
            warn <- ""
            df <- read.table("${data_assembly_ch}", header = TRUE, sep = "\\t")
            meta <- read.table("./${meta_file}", sep = "\\t", header = TRUE)
            meta_seq_names <- "${meta_seq_names}" # conversion here because if NULL, block the code
            meta_name_replacement <- "${meta_name_replacement}" # conversion here because if NULL, block the code
            meta_legend <- "${meta_legend}" # conversion here because if NULL, block the code
            if( ! meta_seq_names %in% names(meta)){
                stop(paste0("\\n\\n============\\n\\nERROR IN THE Metadata_check PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS NOT \\"NULL\\", THEN THE meta_seq_names PARAMETER MUST BE A COLUMN NAME OF THE METADATA FILE.\\n\\n============\\n\\n"), call. = FALSE)
            }
            if(meta_name_replacement == "NULL" & meta_legend == "NULL"){
                stop(paste0("\\n\\n============\\n\\nERROR IN THE Metadata_check PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS NOT \\"NULL\\", THEN THE meta_name_replacement AND meta_legend PARAMETERS CANNOT BE BOTH \\"NULL\\".\\n\\n============\\n\\n"), call. = FALSE)
            }
            if(names(df)[1] != "sequence_id"){
                stop(paste0("\\n\\n============\\n\\nINTERNAL ERROR IN THE Metadata_check PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS NOT \\"NULL\\", THEN THE TABLE GENERATED USING THE sample_path PARAMETER MUST CONTAIN sequence_id AS FIRST COLUMN NAME.\\n\\n============\\n\\n"), call. = FALSE)
            }
            if(meta_name_replacement != "NULL"){
                if( ! meta_name_replacement %in% names(meta)){
                    stop(paste0("\\n\\n============\\n\\nERROR IN THE Metadata_check PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS NOT \\"NULL\\" AND IF THE meta_name_replacement PARAMETER IS NOT \\"NULL\\", THEN meta_name_replacement MUST BE A COLUMN NAME OF THE METADATA FILE.\\n\\n============\\n\\n"), call. = FALSE)
                }
                if(names(df)[2] != "initial_sequence_id"){
                    stop(paste0("\\n\\n============\\n\\nINTERNAL ERROR IN THE Metadata_check PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS NOT \\"NULL\\" AND IF THE meta_name_replacement PARAMETER IS NOT \\"NULL\\", THEN THE TABLE GENERATED USING THE sample_path PARAMETER MUST CONTAIN initial_sequence_id AS FIRST COLUMN NAME.\\n\\n============\\n\\n"), call. = FALSE)
                }
                if( ! any(meta[ , meta_name_replacement] %in% df[ , 1])){
                   tempo_warn <- paste0("\\n\\n\\n\\nWARNING IN THE Metadata_check PROCESS OF NEXTFLOW\\nTHE meta_file AND meta_name_replacement PARAMETERS OF THE nextflow.config FILE ARE NOT NULL BUT NO NAME REPLACEMENT PERFORMED\\nPROBABLY THAT THE ${meta_seq_names} COLUMN OF THE FILE INDICATED IN THE meta_path PARAMETER IS NOT MADE OF NAMES OF FASTA FILES INDICATED IN THE sample_path PARAMETER\\nFIRST ELEMENTS OF THE ${meta_seq_names} COLUMN (meta_seq_names PARAMETER) OF THE META DATA FILE ARE:\\n", paste(head(meta[ , meta_seq_names], 20), collapse = "\\n"), "\\nFIRST FASTA FILES NAMES ARE:\\n", paste(head(df[ , 1], 20), collapse = "\\n"), "\\n\\n\\n\\n")
                    warn <- base::paste0(base::ifelse(test = warn == "", yes = tempo_warn, no = base::paste0(warn, "\n\n", tempo_warn, collapse = NULL, recycle0 = FALSE)), collapse = NULL, recycle0 = FALSE)
                    print(tempo_warn)
                }
            }
            if(meta_legend != "NULL"){
                if( ! meta_legend %in% names(meta)){
                    stop(paste0("\\n\\n============\\n\\nERROR IN THE Metadata_check PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS NOT \\"NULL\\" AND IF THE meta_legend PARAMETER IS NOT \\"NULL\\", THEN THE meta_legend PARAMETER MUST BE A COLUMN NAME OF THE METADATA FILE.\\n\\n============\\n\\n"), call. = FALSE)
                }
                if( ! meta_legend %in% names(df)){
                    stop(paste0("\\n\\n============\\n\\nINTERNAL ERROR IN THE Metadata_check PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS NOT \\"NULL\\" AND IF THE meta_legend PARAMETER IS NOT \\"NULL\\", THEN meta_legend MUST BE A COLUMN NAME OF clone_assigned_seq.tsv.\\n\\n============\\n\\n"), call. = FALSE)
                }
            }
            writeLines(warn, con = "warnings.txt")
        ' |& tee -a metadata_check.log
    else
        echo "" > metadata_check.log
        echo "" > warnings.txt
    fi
    """
}
