


// Rename the sequence names in the data file
// The sequence names are replaced by the values of the column meta_name_replacement (specified in .config) of the metadata file
process Seq_name_replacement {
    label 'r_ig_clustering'
    cache 'true'

    input:
    path add_aa_imgt_ch // parallelization expected
    path meta_file
    val meta_seq_names
    val meta_name_replacement
    val meta_legend

    output:
    path "*_renamed_seq.tsv", emit: seq_name_replacement_ch
    path "seq_name_replacement.log", emit: seq_name_replacement_log_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    FILENAME=\$(basename -- ${add_aa_imgt_ch}) # recover a file name without path
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a seq_name_replacement.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a seq_name_replacement.log
    # check first that the data file does not have the second column name starting by "initial_". Otherwise, with create unproper behavior in donut
    if [[ "${meta_file}" == "NULL" ]] ; then
        rm NULL # remove the initial file to avoid to send it into the channel
        echo -n "" > NULL # new hard file that can be sent into the channel
        chmod 777 NULL
        Rscript -e '
            seq <- read.table("./${add_aa_imgt_ch}", sep = "\\t", header = TRUE)
            if(grepl(x = names(seq)[2], pattern = "^initial_")){
                stop(paste0("\\n\\n============\\n\\nERROR IN THE Seq_name_replacement PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS \\"NULL\\", THEN THE SECOND COLUMN OF THE DATA IN THE sample_path PARAMETER CANNOT HAVE THE NAME OF THE SECOND COLUNM STARTING BY \\"initial_\\"\\n\\n============\\n\\n"), call. = FALSE)
            }
        ' |& tee -a seq_name_replacement.log
        IFS='_' read -r -a TEMPO <<< "\${FILENAME}" # string split into array
        cat ${add_aa_imgt_ch} > ./\${TEMPO[0]}_renamed_seq.tsv |& tee -a seq_name_replacement.log
    else
        # if [[ "${meta_file}" != "NULL" && "${meta_name_replacement}" != "NULL" ]] ; then # or [[ "${meta_file}" -ne "NULL" && "${meta_name_replacement}" -ne "NULL" ]], but not !=
        Rscript -e '
            meta <- read.table("./${meta_file}", sep = "\\t", header = TRUE)
            seq <- read.table("./${add_aa_imgt_ch}", sep = "\\t", header = TRUE)
            id <- seq[1, 1] # Extract the name of the colum one of seq
            if( ! "${meta_seq_names}" %in% names(meta)){
                stop(paste0("\\n\\n============\\n\\nERROR IN THE Seq_name_replacement PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS NOT \\"NULL\\", THEN THE meta_seq_names PARAMETER MUST BE A COLUMN NAME OF THE METADATA FILE.\\n\\n============\\n\\n"), call. = FALSE)
            }
            if(names(seq)[1] != "sequence_id"){
                stop(paste0("\\n\\n============\\n\\nINTERNAL ERROR IN THE Seq_name_replacement PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS NOT \\"NULL\\", THEN THE TABLE GENERATED USING THE sample_path PARAMETER MUST CONTAIN sequence_id AS FIRST COLUMN NAME.\\n\\n============\\n\\n"), call. = FALSE)
            }
            if("${meta_name_replacement}" != "NULL"){
                if( ! "${meta_name_replacement}" %in% names(meta)){
                    stop(paste0("\\n\\n============\\n\\nERROR IN THE Seq_name_replacement PROCESS OF NEXTFLOW\\nIF NOT \\"NULL\\", THE meta_name_replacement PARAMETER MUST BE A COLUMN NAME OF THE meta_path PARAMETER (METADATA FILE): ", "${meta_name_replacement}", "\\n\\n============\\n\\n"), call. = FALSE)
                }
                seq <- data.frame(seq[1], tempo = seq[1], seq[2:length(seq)]) # the second column is created to keep the initial sequence names, before replacement
                for(i2 in 1:nrow(meta)){
                    if(sum(seq[ , 2] %in% meta[i2, "${meta_seq_names}"]) > 1){
                        stop(paste0("\\n\\n============\\n\\nERROR IN THE Seq_name_replacement PROCESS OF NEXTFLOW\\nIN THE METADATA FILE, A SEQUENCE NAME CANNOT BELONG TO SEVERAL VALUES OF THE meta_name_replacement PARAMETER COLUMN NAME OF THE meta_path PARAMETER\\nTHE METAFILE IS: ${meta_file}\\nTHE COLUM NAME IS: ", "${meta_name_replacement}", "\\nTHE PROBLEMATIC REPLACEMENT NAME IN THE METAFILE IS: ", paste(meta[i2, "${meta_seq_names}"], collapse = " "), "\\n\\n============\\n\\n"), call. = FALSE)
                    }else if(any(seq[ , 2] == meta[i2, "${meta_seq_names}"])){
                        if( ! (meta[i2, "${meta_name_replacement}"] == "" | is.na(meta[i2, "${meta_name_replacement}"]))){
                            seq[seq[ , 2] == meta[i2, "${meta_seq_names}"], 1] <- meta[i2, "${meta_name_replacement}"] # replacement of the name in column 1
                        }
                    }
                }
                names(seq)[2] <- paste0("initial_", names(seq)[1])
            }
            if("${meta_legend}" != "NULL"){
                if( ! "${meta_legend}" %in% names(meta)){
                    stop(paste0("\\n\\n============\\n\\nERROR IN THE Seq_name_replacement PROCESS OF NEXTFLOW\\nIF THE meta_legend PARAMETER IS NOT \\"NULL\\", THEN IT MUST BE A COLUMN NAME OF THE METADATA FILE.\\n\\n============\\n\\n"), call. = FALSE)
                }
                seq <- data.frame(seq, TEMPO = as.character(NA))
                names(seq)[length(seq)] <- "${meta_legend}"
                for(i2 in 1:nrow(meta)){
                    # ifelse(test = "${meta_name_replacement}" == "NULL", yes = 1, no = 2) because 1st column to use or 2nd column
                    if(sum(seq[ , ifelse(test = "${meta_name_replacement}" == "NULL", yes = 1, no = 2)] %in% meta[i2, "${meta_seq_names}"]) > 1 | sum(meta[i2, "${meta_seq_names}"] %in% seq[ , ifelse(test = "${meta_name_replacement}" == "NULL", yes = 1, no = 2)]) > 1){
                        stop(paste0("\\n\\n============\\n\\nERROR IN THE Seq_name_replacement PROCESS OF NEXTFLOW\\nIN THE METADATA FILE, A SEQUENCE NAME CANNOT MATCH SEVERAL NAMES IN THE TABLE PROVIDED BY THE sample_path PARAMETER.\\nTHE METAFILE IS: ${meta_file}\\nTHE COLUM NAME IS: ", "${meta_name_replacement}", "\\nTHE PROBLEMATIC REPLACEMENT NAME IN THE METAFILE IS: ", paste(meta[i2, "${meta_seq_names}"], collapse = " "), "\\n\\n============\\n\\n"), call. = FALSE)
                    }else if(any(seq[ , ifelse(test = "${meta_name_replacement}" == "NULL", yes = 1, no = 2)] == meta[i2, "${meta_seq_names}"])){
                        if( ! (meta[i2, "${meta_legend}"] == "" | is.na(meta[i2, "${meta_legend}"]))){
                            seq[seq[ , ifelse(test = "${meta_name_replacement}" == "NULL", yes = 1, no = 2)] == meta[i2, "${meta_seq_names}"], "${meta_legend}"] <- meta[i2, "${meta_legend}"]
                        }
                    }
                }
            }
            write.table(seq, file = paste0("./", id, "_renamed_seq.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\\t")
            # modification of the metadata file for the correct use of ggtree::"%<+%" in germ_tree_vizu.R that uses the column name meta_seq_names for that 
            # meta <- data.frame(meta, initial_label = meta[ , "${meta_seq_names}"])
            # meta[ , "${meta_seq_names}"] <- meta[ , "${meta_name_replacement}"]
            # write.table(meta, file = "./metadata2.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\\t")
            # end modification of the metadata file for the correct use of ggtree::"%<+%" in germ_tree_vizu.R that uses the column name meta_seq_names for that
        ' |& tee -a seq_name_replacement.log
    fi
    """
}
