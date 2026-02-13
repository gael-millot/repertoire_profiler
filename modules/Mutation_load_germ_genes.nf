
process Mutation_load_germ_genes {
    label 'immcantation'
    //publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{*.log}", overwrite: false
    //publishDir path: "${out_path}", mode: 'copy', pattern: "{*.tsv}", overwrite: false
    cache 'true'

    input:
    path translate_germ_ch // parallelization expected (by clonal groups)
    path meta_file
    val meta_legend
    val clone_mut_obs_seq
    val clone_mut_germ_seq
    val clone_mut_regionDefinition


    output:
    path "*_shm-pass.tsv", emit: mutation_load_ch
    path "*.log", emit: mutation_load_log_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    FILENAME=\$(basename -- ${translate_germ_ch}) # recover a file name without path
    cp -Lr ${translate_germ_ch} "./TEMPO.tsv" # to have the hard file, not the symlink, because modifications will be performed inside
    chmod 777 TEMPO.tsv
    rm \$FILENAME # remove the initial file to avoid to send it into the channel
    cp -rp TEMPO.tsv "\$FILENAME" # -p for preserve permissions
    rm TEMPO.tsv
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a Mutation_load_germ_genes.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a Mutation_load_germ_genes.log
    Rscript -e '
        options(warning.length = 8170)
        args <- commandArgs(trailingOnly = TRUE)  # recover arguments written after the call of the Rscript
        tempo.arg.names <- c("file_name") # objects names exactly in the same order as in the bash code and recovered in args
        if(length(args) != length(tempo.arg.names)){
          tempo.cat <- paste0("\\n\\n========\\n\\nINTERNAL ERROR IN THE NEXTFLOW EXECUTION OF THE Mutation_load_germ_genes PROCESS\\n THE NUMBER OF ELEMENTS IN args (", length(args),") IS DIFFERENT FROM THE NUMBER OF ELEMENTS IN tempo.arg.names (", length(tempo.arg.names),")\\nargs:", paste0(args, collapse = ","), "\\ntempo.arg.names:", paste0(tempo.arg.names, collapse = ","), "\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n")
          stop(tempo.cat)
        }
        for(i2 in 1:length(tempo.arg.names)){
          assign(tempo.arg.names[i2], args[i2])
        }
        VDJ_db <- read.table(file_name, sep = "\\t", header = TRUE)
        if( ! "${clone_mut_obs_seq}" %in% names(VDJ_db)){
            stop(paste0("\\n\\n============\\n\\nERROR IN THE Mutation_load_germ_genes PROCESS OF NEXTFLOW\\nTHE clone_mut_obs_seq PARAMETER OF THE nextflow.config FILE MUST BE A COLUMN NAME OF THE clone_assigned_seq.tsv FILE.\\nPLEASE, CHOOSE AMONG: ", paste0(names(VDJ_db), collapse = " "), "\\n\\n============\\n\\n"), call. = FALSE)
        }
        if( ! "${clone_mut_germ_seq}" %in% names(VDJ_db)){
            stop(paste0("\\n\\n============\\n\\nERROR IN THE Mutation_load_germ_genes PROCESS OF NEXTFLOW\\nTHE clone_mut_germ_seq PARAMETER OF THE nextflow.config FILE MUST BE A COLUMN NAME OF THE clone_assigned_seq.tsv FILE.\\nPLEASE, CHOOSE AMONG: ", paste0(names(VDJ_db), collapse = " "), "\\n\\n============\\n\\n"), call. = FALSE)
        }
        # Calculate R and S mutation counts
        regionDefinition <- if("${clone_mut_regionDefinition}" == "NULL"){NULL}else{eval(parse(text = "shazam::${clone_mut_regionDefinition}"))}
        if( ! (( ! is.null(regionDefinition)) & grepl(x = "${clone_mut_obs_seq}", pattern = "alignment") & grepl(x = "${clone_mut_germ_seq}", pattern = "alignment"))){
            stop(paste0("\\n\\n============\\n\\nERROR IN THE Mutation_load_germ_genes PROCESS OF NEXTFLOW\\nTHE clone_mut_obs_seq AND clone_mut_germ_seq PARAMETERS OF THE nextflow.config FILE MUST BE ALIGNMENT SEQUENCES (COLUMN NAMES OF THE clone_assigned_seq.tsv FILE CONTAINING \\"alignment\\")\\nIF THE clone_mut_regionDefinition PARAMETER IS NOT \\"NULL\\".\\n\\n============\\n\\n"), call. = FALSE)
        }
        VDJ_db <- shazam::observedMutations(VDJ_db, sequenceColumn="${clone_mut_obs_seq}",
                                    germlineColumn="${clone_mut_germ_seq}",
                                    regionDefinition=regionDefinition,
                                    frequency=FALSE, 
                                    nproc=1)
        # Calculate combined R and S mutation counts
        VDJ_db <- shazam::observedMutations(VDJ_db, sequenceColumn="${clone_mut_obs_seq}",
                                    germlineColumn="${clone_mut_germ_seq}",
                                    regionDefinition=regionDefinition,
                                    frequency=FALSE, 
                                    combine=TRUE,
                                    nproc=1)

        # Calculate R and S mutation frequencies
        VDJ_db <- shazam::observedMutations(VDJ_db, sequenceColumn="${clone_mut_obs_seq}",
                                    germlineColumn="${clone_mut_germ_seq}",
                                    regionDefinition=regionDefinition,
                                    frequency=TRUE, 
                                    nproc=1)

        # Calculate combined R and S mutation frequencies
        VDJ_db <- shazam::observedMutations(VDJ_db, sequenceColumn="${clone_mut_obs_seq}",
                                    germlineColumn="${clone_mut_germ_seq}",
                                    regionDefinition=regionDefinition,
                                    frequency=TRUE, 
                                    combine=TRUE,
                                    nproc=1)

        VDJ_db <- dplyr::arrange(VDJ_db, clone_id)
        # write.table(VDJ_db, file = paste0("./", tools::file_path_sans_ext(file_name), "_shm-pass.tsv"), row.names = FALSE, quote = FALSE, sep = "\\t")
        # Reposition the column starting with "initial_" to the second position
        initial_col <- grep("^initial_", names(VDJ_db), value = TRUE)
        if (length(initial_col) > 0) {
            tempo <- VDJ_db[[initial_col]]
            VDJ_db <- VDJ_db[ , !(names(VDJ_db) %in% initial_col)]
            VDJ_db <- data.frame(VDJ_db[1], tempo, VDJ_db[2:length(VDJ_db)])
            names(VDJ_db)[2] <- initial_col
        }
        # Check if a column equal to meta_legend in lowercase exists in 'data' 
        if("${meta_file}" != "NULL" & "${meta_legend}" != "NULL" ){
            tempo_log <- names(VDJ_db) == tolower("${meta_legend}")
            if(any(tempo_log, na.rm = TRUE)){
                names(VDJ_db)[tempo_log] <- "${meta_legend}"
            }
        }
        # end Check if a column equal to meta_legend in lowercase exists in 'data'
        # add germline vdj genes
        required_inputs <- c("germline_v_call", "germline_d_call", "germline_j_call")
        if( !all(required_inputs %in% names(VDJ_db)) ){
            stop(paste0("\\n\\n========\\n\\nERROR IN THE NEXTFLOW EXECUTION OF THE germline_genes PROCESS\\nMISSING germline_v_call, germline_d_call OR germline_j_call FOLLOWING THE GermlineGenes PROCESS (OUTPUT germ_tree_ch)\\n\\n========\\n\\n"), call. = FALSE)
        }
        tempo_v <- strsplit(VDJ_db\$germline_v_call, ",")
        sub_v <- sapply(X = tempo_v, FUN = function(x){y <- sub(pattern = "\\\\*.*", replacement = "", x = x) ; paste0(unique(y), collapse = ",")})
        sub_d <- sapply(VDJ_db\$germline_d_call, function(x) {
            if (is.na(x) || x == "") { # The germline_d_call column can have only empty values if the sequences are light chain
                return(NA)
            } else {
                y <- sub(pattern = "\\\\*.*", replacement = "", x = unlist(strsplit(x, ",")))
                return(paste0(unique(y), collapse = ","))
            }
        })
        tempo_j <- strsplit(VDJ_db\$germline_j_call, ",")
        sub_j <- sapply(X = tempo_j, FUN = function(x){y <- sub(pattern = "\\\\*.*", replacement = "", x = x) ; paste0(unique(y), collapse = ",")})
        VDJ_db <- data.frame(VDJ_db, clonal_germline_v_gene = sub_v, clonal_germline_d_gene = sub_d, clonal_germline_j_gene = sub_j)
        # end add germline vdj genes
        write.table(VDJ_db, file = paste0("./", tools::file_path_sans_ext(file_name), "_shm-pass.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\\t")
    ' "\${FILENAME}" |& tee -a Mutation_load_germ_genes.log
    # cp ./tempo_shm-pass.tsv \${FILENAME}_shm-pass.tsv
    # rm tempo_shm-pass.tsv
    """
}
