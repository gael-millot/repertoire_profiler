
process Closest_germline {
    label 'immcantation'
    //publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{*.log}", overwrite: false
    //publishDir path: "${out_path}", mode: 'copy', pattern: "{*.tsv}", overwrite: false
    cache 'true'

    input:
    path clone_split_ch // parallelization expected (by clonal groups)
    val igblast_organism
    val igblast_variable_ref_files
    val clone_germline_kind
    path meta_file
    val meta_legend

    output:
    tuple path("*_germ-pass.tsv"), path("clonal_germline_seq.fasta"), emit: closest_ch, optional: true
    path "failed_clonal_germline_seq.tsv", emit: failed_clonal_germline_ch
    path "*.log", emit: closest_log_ch

    script:
    """
    #!/bin/bash -ue
    # set -o pipefail # inactivated because ParseDb.py returns error that are handled
    FILENAME=\$(basename -- ${clone_split_ch}) # recover a file name without path
    FILE=\${FILENAME%.*} # file name without extension
    cp -Lr ${clone_split_ch} "./TEMPO.tsv" # to have the hard file, not the symlink, because modifications will be performed inside
    chmod 777 TEMPO.tsv
    rm \$FILENAME # remove the initial file to avoid to send it into the channel
    cp -rp TEMPO.tsv "\$FILENAME" # -p for preserve permissions
    rm TEMPO.tsv
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a Closest_germline.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a Closest_germline.log
    # variables

    REPO_PATH="/usr/local/share/germlines/imgt/${igblast_organism}/vdj" # path where the imgt_human_IGHV.fasta, imgt_human_IGHD.fasta and imgt_human_IGHJ.fasta files are in the docker container
    VDJ_FILES_INI="${igblast_variable_ref_files}"
    VDJ_FILES_INI="\${VDJ_FILES_INI// NULL / }" # remove the string NULL if exists
    VDJ_FILES=\$(awk -v var1="\${VDJ_FILES_INI}" -v var2="\${REPO_PATH}" 'BEGIN{ORS=" " ; split(var1, array1, " ") ; for (key in array1) {print var2"/"array1[key]}}')
    # end variables
    if [[ "${clone_germline_kind}" == "full" ]] ; then
        Rscript -e '
            args = commandArgs(trailingOnly=TRUE)
            db <- read.table(args[1], sep = "\\t", header = TRUE)
            write.table(db[, c("sequence_id", "germline_alignment", "clone_id")], file = "tempo.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\\t") # germline_alignment column data
        ' \$FILENAME |& tee -a Closest_germline.log
    fi
    CreateGermlines.py -d \$FILENAME -g "${clone_germline_kind}" --cloned -r \${VDJ_FILES} |& tee -a Closest_germline.log
    # -r: List of folders and/or fasta files (with .fasta, .fna or .fa extension) with germline sequences. When using the default Change-O sequence and coordinate fields, these reference sequences must contain IMGT-numbering spacers (gaps) in the V segment.
    set -o pipefail
    if [[ -s \${FILE}_germ-pass.tsv ]] ; then # -s means CreateGermlines.py did not fail. Warning wildcard *_germ-pass.tsv does not work.
        if [[ "${clone_germline_kind}" == "full" ]] ; then
            Rscript -e '
                args = commandArgs(trailingOnly=TRUE)
                db1 <- read.table(args[1], sep = "\\t", header = TRUE)
                db2 <- read.table(args[2], sep = "\\t", header = TRUE)
                if( ! all(sort(db1\$sequence_id) == sort(db2\$sequence_id))){
                    stop(paste0("\\n\\n================\\n\\nERROR IN Closest_germline PROCESS.\\ndb1 and db2 DO NOT HAVE THE SAME NAMES IN THE sequence_id COLUMN.\\n\\n================\\n\\n"), call. = FALSE)
                }else{
                    db2 <- db2[match(db2\$sequence_id, db1\$sequence_id), ]
                }
                coord <- which(names(db1) == "germline_v_call")
                tempo <- data.frame(db1[1:(coord - 1)], germline_alignment_full = db1\$germline_alignment, db1[coord:length(db1)]) # add the initial germline_alignment column data in the modified germline_alignment column
                tempo\$germline_alignment <- db2\$germline_alignment
                write.table(tempo, file = args[1], row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\\t")
            ' *_germ-pass.tsv tempo.tsv |& tee -a Closest_germline.log
        fi
        Rscript -e '
            options(warning.length = 8170)
            args = commandArgs(trailingOnly=TRUE)
            db <- read.table(args[1], sep = "\\t", header = TRUE)
            productive <- read.table(args[2], sep = "\\t", header = TRUE)
            # put back meta_legend column name as in productive_seq.tsv, beacause in lowercase
            if("${meta_file}" != "NULL" & "${meta_legend}" != "NULL" ){
                tempo_log1 <- names(db) == tolower("${meta_legend}")
                if(any(tempo_log1, na.rm = TRUE)){
                    names(db)[tempo_log1] <- "${meta_legend}"
                }
                tempo_log2 <- names(productive) == tolower("${meta_legend}")
                if(any(tempo_log2, na.rm = TRUE)){
                    names(productive)[tempo_log2] <- "${meta_legend}"
                }
            }
            # end put back meta_legend column names as in productive_seq.tsv, beacause in lowercase
            # reorder as in productive_seq.tsv
            if( ! all(names(productive) %in% names(db))){
                stop(paste0("\\n\\n================\\n\\nERROR IN Closest_germline PROCESS.\\nNAMES OF productive_seq.tsv SHOULD ALL BE IN THE OUTPUT .tsv FILE OF CreateGermlines.py.\\nNAMES OF productive_seq.tsv:\\n", paste0(sort(names(productive)), collapse = " "), "\\nNAMES OF THE OUTPUT:\\n", paste0(sort(names(db)), collapse = " "), "\\n\\n================\\n\\n"), call. = FALSE)
            }else{
                tempo_log <- names(db) %in% names(productive)
                tempo_db1 <- db[tempo_log]
                tempo_db2 <- db[ ! tempo_log]
                tempo_db1 <- tempo_db1[ , match(names(productive), names(tempo_db1))] # reorder as in productive_seq.tsv
                db <- data.frame(tempo_db1, tempo_db2)
            }
            # end reorder as in productive_seq.tsv
            if("${clone_germline_kind}" == "dmask"){
                germline_col_name <- "germline_alignment_d_mask"
            }else if("${clone_germline_kind}" == "vonly"){
                germline_col_name <- "germline_alignment_v_region"
            }else if("${clone_germline_kind}" == "full"){
                germline_col_name <- "germline_alignment_full"
            }
            if(length(unique(db[, germline_col_name])) != 1){
                stop(paste0("\\n\\n================\\n\\nERROR IN Closest_germline PROCESS.\\nSEQUENCES ARE NOT IDENTICAL IN THE ", germline_col_name, " COLUMN:\\n", paste0(db[, germline_col_name], collapse = "\\n"), "\\n\\n================\\n\\n"), call. = FALSE)
            }
            # add clonal_germline_seq columns
            db <- data.frame(db, TEMPO_NAME_1 = db[ , germline_col_name])
            names(db)[names(db) == "TEMPO_NAME_1"] <- "clonal_germline_sequence_with_gaps"
            db <- data.frame(db, TEMPO_NAME_2 = gsub(x = db[ , germline_col_name], pattern = "\\\\.", replacement = ""))
            names(db)[names(db) == "TEMPO_NAME_2"] <- "clonal_germline_sequence_no_gaps"
            # end add clonal_germline_seq columns
            # make fasta file
            fasta <- file("clonal_germline_seq.fasta", "w")
            writeLines(paste0(">clone_id_", db\$clone_id[1]), fasta) # first line taken because clonal germline seq are identical
            writeLines(db[1, "clonal_germline_sequence_no_gaps"], fasta)
            close(fasta)
            # end make fasta file
            # move clone_id
            clone_id <- db\$clone_id
            db <- db[ , ! names(db) %in% "clone_id"] # remove clone_id column
            pos <- which(names(db) == "c_gene")
            db <- data.frame(db[1:pos], clone_id = clone_id, db[(pos + 1):length(names(db))]) # reposition clone_id
            # end move clone_id
            write.table(db, file = args[1], row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\\t")
        ' *_germ-pass.tsv \$FILENAME |& tee -a Closest_germline.log
        echo -e "\\n\\nNOTE: EMPTY failed_clonal_germline_seq.tsv FILE RETURNED FOLLOWING THE Closest_germline PROCESS\\n\\n" |& tee -a Closest_germline.log
        head -1 \${FILE}_germ-pass.tsv | cat > failed_clonal_germline_seq.tsv # keep the header
    else
        cp \$FILENAME failed_clonal_germline_seq.tsv
        Rscript -e '
            args = commandArgs(trailingOnly=TRUE)
            db <- read.table(args[1], sep = "\\t", header = TRUE)
            productive <- read.table(args[2], sep = "\\t", header = TRUE)
            # put back meta_legend column name as in productive_seq.tsv, beacause in lowercase
            if("${meta_file}" != "NULL" & "${meta_legend}" != "NULL" ){
                tempo_log1 <- names(db) == tolower("${meta_legend}")
                if(any(tempo_log1, na.rm = TRUE)){
                    names(db)[tempo_log1] <- "${meta_legend}"
                }
                tempo_log2 <- names(productive) == tolower("${meta_legend}")
                if(any(tempo_log2, na.rm = TRUE)){
                    names(productive)[tempo_log2] <- "${meta_legend}"
                }
            }
            # end put back meta_legend column names as in productive_seq.tsv, beacause in lowercase
            # reorder as in productive_seq.tsv
            if( ! all(names(productive) %in% names(db))){
                stop(paste0("\\n\\n================\\n\\nERROR IN Closest_germline PROCESS.\\nNAMES OF productive_seq.tsv SHOULD ALL BE IN THE OUTPUT .tsv FILE OF CreateGermlines.py.\\nNAMES OF productive_seq.tsv:\\n", sort(names(productive)), "\\nNAMES OF THE OUTPUT:\\n", sort(names(db)), "\\n\\n================\\n\\n"), call. = FALSE)
            }else{
                tempo_log <- names(db) %in% names(productive)
                tempo_db1 <- db[tempo_log]
                tempo_db2 <- db[ ! tempo_log]
                tempo_db1 <- tempo_db1[ , match(names(productive), names(tempo_db1))] # reorder as in productive_seq.tsv
                db <- data.frame(tempo_db1, tempo_db2)
            }
            # end reorder as in productive_seq.tsv
            write.table(db, file = args[1], row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\\t")
        ' failed_clonal_germline_seq.tsv \$FILENAME |& tee -a Closest_germline.log
    fi
    """
}

