nextflow.enable.dsl=2
/*
#########################################################################
##                                                                     ##
##     main.nf of repertoire profiler                                  ##
##                                                                     ##
##     Gael A. Millot                                                  ##
##     Chloe Taurel                                                    ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Computational Biology Department                                ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################
*/

//////// Processes


process workflowParam { // create a file with the workflow parameters in out_path
    label 'bash'
    publishDir "${out_path}/reports", mode: 'copy', overwrite: false
    cache 'false'

    input:
    val modules

    output:
    path "Run_info.txt"

    script:
    """
    echo "Project (empty means no .git folder where the main.nf file is present): " \$(git -C ${projectDir} remote -v | head -n 1) > Run_info.txt # works only if the main script run is located in a directory that has a .git folder, i.e., that is connected to a distant repo
    echo "Git info (empty means no .git folder where the main.nf file is present): " \$(git -C ${projectDir} describe --abbrev=10 --dirty --always --tags) >> Run_info.txt # idem. Provide the small commit number of the script and nextflow.config used in the execution
    echo "Cmd line: ${workflow.commandLine}" >> Run_info.txt
    echo "execution mode": ${system_exec} >> Run_info.txt
    modules=$modules # this is just to deal with variable interpretation during the creation of the .command.sh file by nextflow. See also \$modules below
    if [[ ! -z \$modules ]] ; then
        echo "loaded modules (according to specification by the user thanks to the --modules argument of main.nf): ${modules}" >> Run_info.txt
    fi
    echo "Manifest's pipeline version: ${workflow.manifest.version}" >> Run_info.txt
    echo "result path: ${out_path}" >> Run_info.txt
    echo "nextflow version: ${nextflow.version}" >> Run_info.txt
    echo -e "\\n\\nIMPLICIT VARIABLES:\\n\\nlaunchDir (directory where the workflow is run): ${launchDir}\\nprojectDir (directory where the main.nf script is located): ${projectDir}\\nworkDir (directory where tasks temporary files are created): ${workDir}" >> Run_info.txt
    echo -e "\\n\\nUSER VARIABLES:\\n\\nout_path: ${out_path}\\nsample_path: ${sample_path}" >> Run_info.txt
    """
}
//${projectDir} nextflow variable
//${workflow.commandLine} nextflow variable
//${workflow.manifest.version} nextflow variable
//Note that variables like ${out_path} are interpreted in the script block



process Unzip {
    label 'unzip'
    cache 'true'

    input:
    path zip
    val sample_path

    output:
    path "*", emit: unzip_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail # return the exit code of the first nonzero command in the pipeline, not just the last one (important when using tee)
    FILENAME_INI=\$(basename -- "${sample_path}")
    FILE_EXTENSION="\${FILENAME_INI##*.}"
    FILENAME="\${FILENAME_INI%.*}"
    if [[ ! "\${FILE_EXTENSION}" =~ zip ]] ; then
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nTHE FILE EXTENSION MUST BE \\".zip\\" AND NOT \${FILENAME_INI}\\n\\n========\\n\\n"
        exit 1
    else
        unzip ${zip}
    fi
    rm \$FILENAME_INI
    for file in "\$FILENAME"/*.* ; do
        # Check if the file is a regular file and not a directory
        if [ -f "\$file" ] ; then
            # Check if the file does not match the extensions .fasta or .fa
            if [[ ! "\$file" =~ \\.(fasta|fa|fas|fna|txt|seq)\$ ]] ; then
                echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nALL THE UNZIPPED FILE EXTENSION MUST BE fasta, fa, fna, txt, seq OR faa.\\nHERE A FILE IS:\\n\$file.\\n\\n========\\n\\n"
                exit 1
            fi
        else
            echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nTHE UNZIPPED FOLDER MUST CONTAIN ONLY FILES WITH THE FOLLOWING EXTENSION: fasta, fa, fna OR faa\\n\\n========\\n\\n"
            exit 1
        fi
    done
    cp "\$FILENAME"/*.* .
    rm -r \$FILENAME
    """
}



// next process is for repertoire
// It currently only accepts IG reference files and 'mouse' or 'human' species
process igblast_data_check { // cannot be igblast_data_check outside of process because the files are present in the docker
    label 'immcantation'
    //publishDir path: "${out_path}", mode: 'copy', pattern: "{*.tsv}", overwrite: false

    cache 'true'

    input:
    val igblast_organism
    val igblast_variable_ref_files
    val igblast_constant_ref_files

    output:
    path "*.tsv", emit: igblast_data_check_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    REPO_PATH_VAR="/usr/local/share/germlines/imgt/${igblast_organism}/vdj" # path where the imgt .fasta reference seq files are in the docker container
    REPO_PATH_CONST="/usr/local/share/germlines/imgt/${igblast_organism}/constant" # path where the imgt .fasta reference seq files are in the docker container
    VDJ_FILES=\$(awk -v var1="${igblast_variable_ref_files}" -v var2="\${REPO_PATH_VAR}" 'BEGIN{ORS=" " ; split(var1, array1, " ") ; for (key in array1) {print var2"/"array1[key]}}') # assemble files with their path
    CONST_FILES=\$(awk -v var1="${igblast_constant_ref_files}" -v var2="\${REPO_PATH_CONST}" 'BEGIN{ORS=" " ; split(var1, array1, " ") ; for (key in array1) {print var2"/"array1[key]}}') # assemble files with their path
    for i1 in \$VDJ_FILES ; do
        if [[ ! -e \${i1} ]] ; then
            echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nFILE DOES NOT EXISTS:\\n\${i1}\\n\\nINDICATED PATH:\\n\${REPO_PATH_VAR}\\n\\nCONTAINS:\\n"
            ls -la -w 1 \${REPO_PATH_VAR}
            echo -e "\\n\\n========\\n\\n"
            exit 1
        else
            FILENAME=\$(basename -- "\${i1}") # recover a file name without path
            FILENAME="\${FILENAME%.*}" # remove extension
            grep -E '^>.*\$' \${i1} | cut -f2 -d'|' > \${FILENAME}.tsv # detect line starting by > and extract the 2nd field after cutting by |
        fi
    done
    for i1 in \$CONST_FILES ; do
        if [[ ! -e \${i1} ]] ; then
            echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nFILE DOES NOT EXISTS:\\n\${i1}\\n\\nINDICATED PATH:\\n\${REPO_PATH_CONST}\\n\\nCONTAINS:\\n"
            ls -la -w 1 \${REPO_PATH_CONST}
            echo -e "\\n\\n========\\n\\n"
            exit 1
        else
            FILENAME=\$(basename -- "\${i1}") # recover a file name without path
            FILENAME="\${FILENAME%.*}" # remove extension
            grep -E '^>.*\$' \${i1} | cut -f2 -d'|' > \${FILENAME}.tsv # detect line starting by > and extract the 2nd field after cutting by |
        fi
    done
    """
    // write ${} between "" to make a single argument when the variable is made of several values separated by a space. Otherwise, several arguments will be considered
}






// Parse the tsv files generated by the igblast process to select the productive and unproductive sequences
// Output:
// - productive_seq_init.tsv: tsv file with the productive sequences
// - failed_productive_seq.tsv: tsv file with the unproductive sequences
process parseDb_filtering {
    label 'immcantation'
    cache 'true'

    input:
    path db_pass_ch // parallelization expected

    output:
    path "productive_seq_init.tsv", emit: select_ch, optional: true
    path "failed_productive_seq.tsv", emit: unselect_ch
    path "ParseDb_filtering.log", emit: parseDb_filtering_log_ch

    script:
    """
    #!/bin/bash -ue
    # set -o pipefail # inactivated because ParseDb.py returns error that are handled
    FILENAME=\$(basename -- "${db_pass_ch}") # recover a file name without path
    FILE=\${FILENAME%.*} # file name without extension
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a ParseDb_filtering.log
    if [[ -s ${db_pass_ch} ]]; then # -s means "exists and non empty". Thus, return FALSE is the file does not exists or is empty
        # ParseDb.py select -d ${db_pass_ch} -f productive -u TRUE T |& tee -a ParseDb_filtering.log
        ParseDb.py split -d ${db_pass_ch} -f productive |& tee -a ParseDb_filtering.log  # Used to create 2 files with content that depends of if the FALSE F or TRUE T value is found in the productive field (select command only creates a select file if values specifies in -u flag are found, in this case TRUE or T). If only F in the input file, _productive-T.tsv is not created
        if [ -s \${FILE}_productive-T.tsv ]; then
            cp \${FILE}_productive-T.tsv productive_seq_init.tsv |& tee -a ParseDb_filtering.log # can be empty file (only header)
        elif [ -s \${FILE}_productive-TRUE.tsv ]; then
            cp \${FILE}_productive-TRUE.tsv productive_seq_init.tsv |& tee -a ParseDb_filtering.log # can be empty file (only header)
        fi
        if [ -s \${FILE}_productive-F.tsv ]; then # see above for -s
            cp \${FILE}_productive-F.tsv failed_productive_seq.tsv |& tee -a ParseDb_filtering.log
        elif [ -s \${FILE}_productive-FALSE.tsv ]; then # if not TRUE or T, the value in the "productive" field can be either F or FALSE
            cp \${FILE}_productive-FALSE.tsv failed_productive_seq.tsv |& tee -a ParseDb_filtering.log
        elif [ -s \${FILE}_productive-NA.tsv ]; then # if not TRUE or T, the value in the "productive" field can be either F or FALSE
            cp \${FILE}_productive-NA.tsv failed_productive_seq.tsv |& tee -a ParseDb_filtering.log
        else
            echo -e "\\n\\nEMPTY failed_productive_seq.tsv FILE RETURNED FOLLOWING THE parseDb_filtering PROCESS\\n\\n" |& tee -a ParseDb_filtering.log
            # echo -n "" | cat > failed_productive_seq.tsv
            head -1 ${db_pass_ch} | cat > failed_productive_seq.tsv # only header in the file
        fi
    else
        echo -e "\\n\\n========\\n\\nINTERNAL ERROR IN NEXTFLOW EXECUTION\\n\\nEMPTY FILE GENERATED BY THE igblast PROCESS\\nCHECK THE ParseDb_filtering.log FILE IN THE report FOLDER INSIDE THE OUTPUT FOLDER\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n" # I said empty because the existence of the file has been checked after the igblast process
        exit 1
    fi
    """
}


// Trim the 5' end of the sequence if required (to remove the leader peptide) and Translate the nucleotidic sequences into amino acid sequences
process TrimTranslate {
    label 'seqkit'
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "productive_nuc/trimmed/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "productive_nuc/removed/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "productive_nuc/query/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "productive_nuc/align/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "productive_nuc/align_with_gaps/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "productive_aa/trimmed/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "productive_aa/igblast/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "productive_aa/query/*.fasta", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "productive_aa/align/*.fasta", overwrite: false
    cache 'true'

    input:
    path select_ch // parallelization expected

    output:
    path "trimtranslate.tsv", emit: trimtranslate_ch // productive file with column sequence_alignment_aa added
    path "productive_nuc/trimmed/*.*"
    path "productive_nuc/removed/*.*"
    path "productive_nuc/query/*.*"
    path "productive_nuc/align/*.*"
    path "productive_nuc/align_with_gaps/*.*"
    path "productive_aa/trimmed/*.*"
    path "productive_aa/igblast/*.*"
    path "productive_aa/query/*.*"
    path "productive_aa/align/*.*"
    path "trimtranslate.log", emit: trimtranslate_log_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    trimtranslate.sh ${select_ch} # |& tee -a trimtranslate.log not used because in trimtranslate.sh
    """
}

// Rename the sequence names in the data file
// The sequence names are replaced by the values of the column meta_name_replacement (specified in .config) of the metadata file
process seq_name_replacement {
    label 'r_ig_clustering'
    cache 'true'

    input:
    path trimtranslate_ch // parallelization expected
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
    FILENAME=\$(basename -- ${trimtranslate_ch}) # recover a file name without path
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a seq_name_replacement.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a seq_name_replacement.log
    # check first that the data file does not have the second column name starting by "initial_". Otherwise, with create unproper behavior in donut
    if [[ "${meta_file}" == "NULL" ]] ; then
        rm NULL # remove the initial file to avoid to send it into the channel
        echo -n "" > NULL # new hard file that can be sent into the channel
        chmod 777 NULL
        Rscript -e '
            seq <- read.table("./${trimtranslate_ch}", sep = "\\t", header = TRUE)
            if(grepl(x = names(seq)[2], pattern = "^initial_")){
                stop(paste0("\\n\\n============\\n\\nERROR IN THE seq_name_replacement PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS \\"NULL\\", THEN THE SECOND COLUMN OF THE DATA IN THE sample_path PARAMETER CANNOT HAVE THE NAME OF THE SECOND COLUNM STARTING BY \\"initial_\\"\\n\\n============\\n\\n"), call. = FALSE)
            }
        ' |& tee -a seq_name_replacement.log
        IFS='_' read -r -a TEMPO <<< "\${FILENAME}" # string split into array
        cat ${trimtranslate_ch} > ./\${TEMPO[0]}_renamed_seq.tsv |& tee -a seq_name_replacement.log
    else
        # if [[ "${meta_file}" != "NULL" && "${meta_name_replacement}" != "NULL" ]] ; then # or [[ "${meta_file}" -ne "NULL" && "${meta_name_replacement}" -ne "NULL" ]], but not !=
        Rscript -e '
            meta <- read.table("./${meta_file}", sep = "\\t", header = TRUE)
            seq <- read.table("./${trimtranslate_ch}", sep = "\\t", header = TRUE)
            id <- seq[1, 1] # Extract the name of the colum one of seq
            if( ! "${meta_seq_names}" %in% names(meta)){
                stop(paste0("\\n\\n============\\n\\nERROR IN THE seq_name_replacement PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS NOT \\"NULL\\", THEN THE meta_seq_names PARAMETER MUST BE A COLUMN NAME OF THE METADATA FILE.\\n\\n============\\n\\n"), call. = FALSE)
            }
            if(names(seq)[1] != "sequence_id"){
                stop(paste0("\\n\\n============\\n\\nINTERNAL ERROR IN THE seq_name_replacement PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS NOT \\"NULL\\", THEN THE TABLE GENERATED USING THE sample_path PARAMETER MUST CONTAIN sequence_id AS FIRST COLUMN NAME.\\n\\n============\\n\\n"), call. = FALSE)
            }
            if("${meta_name_replacement}" != "NULL"){
                if( ! "${meta_name_replacement}" %in% names(meta)){
                    stop(paste0("\\n\\n============\\n\\nERROR IN THE seq_name_replacement PROCESS OF NEXTFLOW\\nIF NOT \\"NULL\\", THE meta_name_replacement PARAMETER MUST BE A COLUMN NAME OF THE meta_path PARAMETER (METADATA FILE): ", "${meta_name_replacement}", "\\n\\n============\\n\\n"), call. = FALSE)
                }
                seq <- data.frame(seq[1], tempo = seq[1], seq[2:length(seq)]) # the second column is created to keep the initial sequence names, before replacement
                for(i2 in 1:nrow(meta)){
                    if(sum(seq[ , 2] %in% meta[i2, "${meta_seq_names}"]) > 1){
                        stop(paste0("\\n\\n============\\n\\nERROR IN THE seq_name_replacement PROCESS OF NEXTFLOW\\nIN THE METADATA FILE, A SEQUENCE NAME CANNOT BELONG TO SEVERAL VALUES OF THE meta_name_replacement PARAMETER COLUMN NAME OF THE meta_path PARAMETER\\nTHE METAFILE IS: ${meta_file}\\nTHE COLUM NAME IS: ", "${meta_name_replacement}", "\\nTHE PROBLEMATIC REPLACEMENT NAME IN THE METAFILE IS: ", paste(meta[i2, "${meta_seq_names}"], collapse = " "), "\\n\\n============\\n\\n"), call. = FALSE)
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
                    stop(paste0("\\n\\n============\\n\\nERROR IN THE seq_name_replacement PROCESS OF NEXTFLOW\\nIF THE meta_legend PARAMETER IS NOT \\"NULL\\", THEN IT MUST BE A COLUMN NAME OF THE METADATA FILE.\\n\\n============\\n\\n"), call. = FALSE)
                }
                seq <- data.frame(seq, TEMPO = as.character(NA))
                names(seq)[length(seq)] <- "${meta_legend}"
                for(i2 in 1:nrow(meta)){
                    # ifelse(test = "${meta_name_replacement}" == "NULL", yes = 1, no = 2) because 1st column to use or 2nd column
                    if(sum(seq[ , ifelse(test = "${meta_name_replacement}" == "NULL", yes = 1, no = 2)] %in% meta[i2, "${meta_seq_names}"]) > 1 | sum(meta[i2, "${meta_seq_names}"] %in% seq[ , ifelse(test = "${meta_name_replacement}" == "NULL", yes = 1, no = 2)]) > 1){
                        stop(paste0("\\n\\n============\\n\\nERROR IN THE seq_name_replacement PROCESS OF NEXTFLOW\\nIN THE METADATA FILE, A SEQUENCE NAME CANNOT MATCH SEVERAL NAMES IN THE TABLE PROVIDED BY THE sample_path PARAMETER.\\nTHE METAFILE IS: ${meta_file}\\nTHE COLUM NAME IS: ", "${meta_name_replacement}", "\\nTHE PROBLEMATIC REPLACEMENT NAME IN THE METAFILE IS: ", paste(meta[i2, "${meta_seq_names}"], collapse = " "), "\\n\\n============\\n\\n"), call. = FALSE)
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


// Add the dist_to_nearest values and gene columns in the seq_name_replacement file, thereby creating the productive_seq.tsv file
// The gene columns are created with the allele columns, minus the allele part
// + isotype class (4 first characters of the c_gene column)
process data_assembly {
    label 'immcantation'
    publishDir path: "${out_path}/tsv", mode: 'copy', pattern: "{productive_seq.tsv}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{data_assembly.log}", overwrite: false
    cache 'true'

    input:
    path seq_name_replacement_ch2 // no parallelization
    path distToNearest_ch

    output:
    path "productive_seq.tsv", emit: productive_ch
    path "data_assembly.log"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    Rscript -e '
        db <- read.table("${seq_name_replacement_ch2}", header = TRUE, sep = "\\t")
        dtn <- read.table("${distToNearest_ch}", header = TRUE, sep = "\\t")


        # replace TRUE and FALSE of the read.table() conversion by the initial T and F
        # tempo.log <- sapply(db, FUN = function(x){class(x) == "logical"})
        # db[tempo.log] <- lapply(db[tempo.log], FUN = function(x){substr(as.character(x), 1, 1)})
        # end replace TRUE and FALSE of the read.table() conversion by the initial T and F


        if( ! any(c("sequence_id", "initial_sequence_id") %in% names(db))){
            tempo.cat <- paste0("\\n\\n========\\n\\nINTERNAL ERROR IN THE NEXTFLOW EXECUTION OF THE data_assembly PROCESS\\ndb SHOULD HAVE \\"sequence_id\\" AND ALSO POTENTIALLY \\"initial_sequence_id\\" AS COLUMN NAME. HERE:\\nNAMES: ", paste(names(db), collapse = " "), "\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n")
            stop(tempo.cat)
        }else{
            tempo.col.name <- ifelse("initial_sequence_id" %in% names(db), "initial_sequence_id", "sequence_id")
        }
        if(all(c("sequence_id", "initial_sequence_id") %in% names(db))){
            if(all(db\$sequence_id == db\$initial_sequence_id)){
                tempo.cat <- paste0("\\n\\n========\\n\\nERROR IN THE data_assembly PROCESS OF NEXTFLOW\\nTHE meta_path AND meta_name_replacement PARAMETERS ARE NOT \\"NULL\\" BUT NO SEQUENCE NAMES HAVE BEEN REPLACED WHEN USING THE meta_name_replacement COLUMN\\n\\n========\\n\\n")
                stop(tempo.cat)
            }
        }
        if(any(is.na(match(db[, tempo.col.name], dtn\$sequence_id)))){
            tempo.cat <- paste0("\\n\\n========\\n\\nINTERNAL ERROR IN THE NEXTFLOW EXECUTION OF THE data_assembly PROCESS\\nNO NA SHOULD APPEAR AT THAT STAGE WITH match()\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n")
            stop(tempo.cat)
        }
        if(sum(db[, tempo.col.name] %in% dtn\$sequence_id, na.rm = TRUE) != nrow(db)){
            tempo.cat <- paste0("\\n\\n========\\n\\nINTERNAL ERROR IN THE NEXTFLOW EXECUTION OF THE data_assembly PROCESS\\nsum(db[, tempo.col.name] %in% dtn\$sequence_id, na.rm = TRUE) SHOULD BE EQUAL TO nrow(db)\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n")
            stop(tempo.cat)
        }
        if(nrow(dtn) < nrow(db)){
            tempo.cat <- paste0("\\n\\n========\\n\\nINTERNAL ERROR IN THE NEXTFLOW EXECUTION OF THE data_assembly PROCESS\\ndtn CANNOT HAVE LESS ROWS THAN db. HERE:\\ndb: ", nrow(db), "\\ndtn: ", nrow(dtn), "\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n")
            stop(tempo.cat)
        }
        # selection of rows and ordering of dtn
        tempo.dtn <- dtn[dtn\$sequence_id %in% db[, tempo.col.name], ] # nb rows of dtn reduced to the one of db
        if(nrow(tempo.dtn) != nrow(db)){
            tempo.cat <- paste0("\\n\\n========\\n\\nINTERNAL ERROR IN THE NEXTFLOW EXECUTION OF THE data_assembly PROCESS\\ntempo.dtn SHOULD HAVE THE SAME NUMBER OF ROWS AS db. HERE:\\ndb: ", nrow(db), "\\ntempo.dtn: ", nrow(tempo.dtn), "\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n")
            stop(tempo.cat)
        }
        tempo.dtn <- tempo.dtn[match(db[, tempo.col.name], tempo.dtn\$sequence_id), ]
        if( ! all(db[, tempo.col.name] == tempo.dtn\$sequence_id)){
            tempo.cat <- paste0("\\n\\n========\\n\\nINTERNAL ERROR IN THE NEXTFLOW EXECUTION OF THE data_assembly PROCESS\\n", tempo.col.name, " COLUMNS OF db AND sequence_id COLUMN OF dtn SHOULD BE IDENTICAL. HERE THEY ARE\\ndb\$", tempo.col.name, ":\\n", paste0(db[, tempo.col.name], collapse = ","), "\\ndtn\$sequence_id:\\n", paste0(dtn\$sequence_id, collapse = ","), "\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n")
            stop(tempo.cat)
        }
        db3 <- data.frame(db, dist_nearest = tempo.dtn\$dist_nearest)
        if( ! "c_call" %in% names(db)){
            tempo.cat <- paste0("\\n\\n========\\n\\nINTERNAL ERROR IN THE NEXTFLOW EXECUTION OF THE data_assembly PROCESS\\ndb SHOULD HAVE \\"c_call\\" AS COLUMN NAME. HERE:\\nNAMES: ", paste(names(db), collapse = " "), "\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n")
            stop(tempo.cat)
        }

        #### remove allele info

        tempo_v <- strsplit(db3\$v_call, ",")
        sub_v <- sapply(X = tempo_v, FUN = function(x){y <- sub(pattern = "\\\\*.*", replacement = "", x = x) ; paste0(unique(y), collapse = ",")})
        tempo_j <- strsplit(db3\$j_call, ",")
        sub_j <- sapply(X = tempo_j, FUN = function(x){y <- sub(pattern = "\\\\*.*", replacement = "", x = x) ; paste0(unique(y), collapse = ",")})
        # The c_call column must be handled with more care because it can have empty values
        subclass <- sapply(db3\$c_call, function(x) {
            if (is.na(x) || x == "") { # The c_call column can have empty values if the sequences are too short
                return(NA)  # Si NA ou "", on affecte NA
            } else {
                # Traitement normal: suppression des suffixes après *
                y <- sub(pattern = "\\\\*.*", replacement = "", x = unlist(strsplit(x, ",")))
                return(paste0(unique(y), collapse = ","))
            }
        })
        class <- sapply(X = subclass, FUN = function(x) {
            if (is.na(x)) {
                return(NA)  # Si x est NA, on garde NA
            } else {
                y <- sub(pattern = "\\\\*.*", replacement = "", x = x) 
                y <- substr(y, 1, 4)
                return(paste0(unique(y), collapse = ","))
            }
        })
        db4 <- data.frame(db3, v_gene = sub_v, j_gene = sub_j, isotype_class = class, c_gene = subclass)

        #### end remove allele info

        write.table(db4, file = paste0("./productive_seq.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\\t")
    ' |& tee -a data_assembly.log
    """
}


process metadata_check { // cannot be in germ_tree_vizu because I have to use the clone_assigned_seq.tsv file for the check
    label 'immcantation'
    cache 'true'

    input:
    path data_assembly_ch // no parallelization
    path meta_file
    val meta_seq_names
    val meta_name_replacement
    val meta_legend

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    if [[ "${meta_file}" != "NULL" ]] ; then
        Rscript -e '
            df <- read.table("${data_assembly_ch}", header = TRUE, sep = "\\t")
            meta <- read.table("./${meta_file}", sep = "\\t", header = TRUE)
            meta_seq_names <- "${meta_seq_names}" # conversion here because if NULL, block the code
            meta_name_replacement <- "${meta_name_replacement}" # conversion here because if NULL, block the code
            meta_legend <- "${meta_legend}" # conversion here because if NULL, block the code
            if( ! meta_seq_names %in% names(meta)){
                stop(paste0("\\n\\n============\\n\\nERROR IN THE metadata_check PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS NOT \\"NULL\\", THEN THE meta_seq_names PARAMETER MUST BE A COLUMN NAME OF THE METADATA FILE.\\n\\n============\\n\\n"), call. = FALSE)
            }
            if(meta_name_replacement == "NULL" & meta_legend == "NULL"){
                stop(paste0("\\n\\n============\\n\\nERROR IN THE metadata_check PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS NOT \\"NULL\\", THEN THE meta_name_replacement AND meta_legend PARAMETERS CANNOT BE BOTH \\"NULL\\".\\n\\n============\\n\\n"), call. = FALSE)
            }
            if(names(df)[1] != "sequence_id"){
                stop(paste0("\\n\\n============\\n\\nINTERNAL ERROR IN THE metadata_check PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS NOT \\"NULL\\", THEN THE TABLE GENERATED USING THE sample_path PARAMETER MUST CONTAIN sequence_id AS FIRST COLUMN NAME.\\n\\n============\\n\\n"), call. = FALSE)
            }
            if(meta_name_replacement != "NULL"){
                if( ! meta_name_replacement %in% names(meta)){
                    stop(paste0("\\n\\n============\\n\\nERROR IN THE metadata_check PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS NOT \\"NULL\\" AND IF THE meta_name_replacement PARAMETER IS NOT \\"NULL\\", THEN meta_name_replacement MUST BE A COLUMN NAME OF THE METADATA FILE.\\n\\n============\\n\\n"), call. = FALSE)
                }
                if(names(df)[2] != "initial_sequence_id"){
                    stop(paste0("\\n\\n============\\n\\nINTERNAL ERROR IN THE metadata_check PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS NOT \\"NULL\\" AND IF THE meta_name_replacement PARAMETER IS NOT \\"NULL\\", THEN THE TABLE GENERATED USING THE sample_path PARAMETER MUST CONTAIN initial_sequence_id AS FIRST COLUMN NAME.\\n\\n============\\n\\n"), call. = FALSE)
                }
                if( ! any(meta[ , meta_name_replacement] %in% df[ , 1])){
                   stop(paste0("\\n\\n============\\n\\nERROR IN THE metadata_check PROCESS OF NEXTFLOW\\nTHE meta_file AND meta_name_replacement PARAMETERS OF THE nextflow.config FILE ARE NOT NULL BUT NO NAME REPLACEMENT PERFORMED\\nPROBABLY THAT THE ${meta_seq_names} COLUMN OF THE FILE INDICATED IN THE meta_path PARAMETER IS NOT MADE OF NAMES OF FASTA FILES INDICATED IN THE sample_path PARAMETER\\nFIRST ELEMENTS OF THE ${meta_seq_names} COLUMN (meta_seq_names PARAMETER) OF THE META DATA FILE ARE:\\n", paste(head(meta[ , meta_seq_names], 20), collapse = "\\n"), "\\nFIRST FASTA FILES NAMES ARE:\\n", paste(head(df[ , 1], 20), collapse = "\\n"), "\\n\\n\\n============\\n\\n"), call. = FALSE)
                }
            }
            if(meta_legend != "NULL"){
                if( ! meta_legend %in% names(meta)){
                    stop(paste0("\\n\\n============\\n\\nERROR IN THE metadata_check PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS NOT \\"NULL\\" AND IF THE meta_legend PARAMETER IS NOT \\"NULL\\", THEN THE meta_legend PARAMETER MUST BE A COLUMN NAME OF THE METADATA FILE.\\n\\n============\\n\\n"), call. = FALSE)
                }
                if( ! meta_legend %in% names(df)){
                    stop(paste0("\\n\\n============\\n\\nINTERNAL ERROR IN THE metadata_check PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS NOT \\"NULL\\" AND IF THE meta_legend PARAMETER IS NOT \\"NULL\\", THEN meta_legend MUST BE A COLUMN NAME OF clone_assigned_seq.tsv.\\n\\n============\\n\\n"), call. = FALSE)
                }
            }
        '
    fi
    """
}

process repertoire {
    label 'r_ig_clustering'
    publishDir path: "${out_path}/pdf", mode: 'copy', pattern: "{*_repertoire.pdf}", overwrite: false
    publishDir path: "${out_path}/figures/png", mode: 'copy', pattern: "{*.png}", overwrite: false
    publishDir path: "${out_path}/figures/svg", mode: 'copy', pattern: "{*.svg}", overwrite: false
    publishDir path: "${out_path}/repertoires", mode: 'copy', pattern: "{rep_*.tsv}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{repertoire.log}", overwrite: false
    cache 'true'

    input:
    path seq_name_replacement_ch2 // no parallelization
    path igblast_data_check_ch
    path cute_file

    output:
    path "*_repertoire.pdf"
    path "*.svg"
    path "*.png", emit: repertoire_png_ch
    path "rep_*.tsv"
    path "repertoire.log"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    repertoire.R \
"${seq_name_replacement_ch2}" \
"${igblast_data_check_ch}" \
"${cute_file}" \
"repertoire.log"
    """
}

// Compares cdr3 différences from sequences with same cdr3 length, same v and same j
// Number of substitutions normalized by the length of the cdr3
// NA means cdr3 sequences are exactly the same for the group
process distToNearest {
    label 'immcantation'
    //publishDir path: "${out_path}", mode: 'copy', pattern: "{nearest_distance.tsv}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{distToNearest.log}", overwrite: false
    cache 'true'

    input:
    path trimtranslate_ch2 // no parallelization
    val clone_model
    val clone_normalize

    output:
    path "nearest_distance.tsv", emit: distToNearest_ch
    path "distToNearest.log"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    Rscript -e '
         # WEIRD stuf: if db alone is returned, and if distToNearest_ch is used for the clone_assignment process and followings, everything is fine. But if db3 is returned with db3 <- data.frame(db, dist_nearest = db2\$dist_nearest) or db3 <- data.frame(db, caca = db2\$dist_nearest) or data.frame(db, caca = db\$sequence_id) or db3 <- data.frame(db, caca = as.numeric(db2\$dist_nearest)) or db3 <- data.frame(db[1:3], caca = db\$sequence_id, db[4:length(db)]), the get_germ_tree process cannot make trees, while the productive.tsv seem identical at the end, between the use of db or db3, except that the clone_id order is not the same
        db <- read.table("${trimtranslate_ch2}", header = TRUE, sep = "\\t")
        db2 <- shazam::distToNearest(db, sequenceColumn = "junction", locusColumn = "locus", model = "${clone_model}", normalize = "${clone_normalize}", nproc = 1)
        write.table(db2, file = paste0("./nearest_distance.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\\t")
    ' |& tee -a distToNearest.log
    """
}


process distance_hist {
    label 'immcantation'
    publishDir path: "${out_path}/figures/png", mode: 'copy', pattern: "{*.png}", overwrite: false
    publishDir path: "${out_path}/figures/svg", mode: 'copy', pattern: "{*.svg}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{distance_hist.log}", overwrite: false
    cache 'true'

    input:
    path distToNearest_ch // no parallelization
    path cute_file
    val clone_model
    val clone_normalize
    val clone_distance

    output:
    path "seq*.pdf", emit: histogram_pdf_ch
    path "*.png", emit: distance_hist_ch // png plot (but sometimes empty) sustematically returned
    path "*.svg"
    path "distance_hist.log"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    histogram.R \
"${distToNearest_ch}" \
"${clone_model}" \
"${clone_normalize}" \
"${clone_distance}" \
"${cute_file}" \
"distance_hist.log"
    """
}

process histogram_assembly {
    label 'r_ig_clustering'
    publishDir path: "${out_path}/pdf", mode: 'copy', pattern: "{seq_distance.pdf}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{histogram_assembly.log}", overwrite: false
    cache 'true'

    input:
    path histogram_pdf_ch

    output:
    path "seq_distance.pdf"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    Rscript -e '
    # assignation to prevent a returned element
        tempo <- qpdf::pdf_combine(input = list.files(path = ".", pattern = "^seq.*.pdf\$"), output = "./seq_distance.pdf")
    ' |& tee -a histogram_assembly.log
    """
}


process clone_assignment {
    label 'immcantation'
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{*.log}", overwrite: false
    publishDir path: "${out_path}/tsv", mode: 'copy', pattern: "failed_clone_assigned_seq.tsv", overwrite: false
    cache 'true'

    input:
    path productive_ch // no parallelization
    val clone_model
    val clone_normalize
    val clone_distance
    val clone_strategy
    path meta_file
    val meta_legend

    output:
    path "*_clone-pass.tsv", emit: clone_ch
    path "failed_clone_assigned_seq.tsv", emit: failed_clone_ch
    path "*.log"

    script:
    """
    #!/bin/bash -ue
    # set -o pipefail # inactivated because DefineClones.py returns error that are handled
    # if [[ -s ${productive_ch} ]]; then # see above for -s # no need anymore
    FILENAME=\$(basename -- ${productive_ch}) # recover a file name without path
    FILE=\${FILENAME%.*} # file name without extension
    DefineClones.py -d ${productive_ch} --act ${clone_strategy} --model ${clone_model} --norm ${clone_normalize} --dist ${clone_distance} --fail |& tee -a clone_assignment.log
    shopt -s nullglob
    files=(\${FILE}_clone-{pass,fail}.tsv) # expend the two possibilities but assign only the existing ones
    cp -Lr "\${files[0]}" "tempo.tsv" # copy whatever failed or succeeded
    Rscript -e '
        options(warning.length = 8170)
        args = commandArgs(trailingOnly=TRUE)
        db <- read.table(args[1], sep = "\\t", header = TRUE)
        productive <- read.table(args[2], sep = "\\t", header = TRUE)
        # put back meta_legend column name as in productive_seq.tsv, beacause in lowercase
        if("${meta_file}" != "NULL" & "${meta_legend}" != "NULL" ){
            tempo_log <- names(db) == tolower("${meta_legend}")
            if(any(tempo_log, na.rm = TRUE)){
                names(db)[tempo_log] <- "${meta_legend}"
            }
        }
        # end put back meta_legend column names as in productive_seq.tsv, beacause in lowercase
        # reorder as in productive_seq.tsv
        if( ! all(names(productive) %in% names(db))){
            stop(paste0("\\n\\n================\\n\\nERROR IN clone_assignment PROCESS.\\nNAMES OF productive_seq.tsv SHOULD ALL BE IN THE OUTPUT .tsv FILE OF DefineClones.py.\\nNAMES OF productive_seq.tsv:\\n", paste0(sort(names(productive)), collapse = " "), "\\nNAMES OF THE OUTPUT:\\n", paste0(sort(names(db)), collapse = " "), "\\n\\n================\\n\\n"), call. = FALSE)
        }else{
            tempo_log <- names(db) %in% names(productive)
            tempo_db1 <- db[tempo_log]
            tempo_db2 <- db[ ! tempo_log]
            tempo_db1 <- tempo_db1[ , match(names(productive), names(tempo_db1))] # reorder as in productive_seq.tsv
            db2 <- data.frame(tempo_db1, tempo_db2)
            write.table(db2, file = "tempo2.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\\t")
        }
        # end reorder as in productive_seq.tsv
    ' tempo.tsv ${productive_ch}
    if [ -s \${FILE}_clone-fail.tsv ]; then # see above for -s
        cp tempo2.tsv failed_clone_assigned_seq.tsv |& tee -a clone_assignment.log
    elif [ -s \${FILE}_clone-pass.tsv ] ; then # see above for -s
        cp -f tempo2.tsv \${FILE}_clone-pass.tsv |& tee -a clone_assignment.log # force overwriting
        set -o pipefail
        echo -e "\\n\\nNOTE: EMPTY failed_clone_assigned_seq.tsv FILE RETURNED FOLLOWING THE clone_assignment PROCESS\\n\\n" |& tee -a clone_assignment.log
        head -1 \${FILE}_clone-pass.tsv | cat > failed_clone_assigned_seq.tsv # header kept
    else
        set -o pipefail
        echo -e "\\n\\n========\\n\\nINTERNAL ERROR IN NEXTFLOW EXECUTION\\n\\nOUTPUT FILE OF DefineClones.py IN THE clone_assignment PROCESS SHOULT BE EITHER *_clone-pass.tsv OR *_clone-fail.tsv.\\nCHECK THE clone_assignment.log IN THE report FOLDER INSIDE THE OUTPUT FOLDER\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n"
        exit 1
    fi
    """
}

process split_by_clones { // split the file into multiple files according to the clone_id column
    label 'immcantation'
    cache 'true'

    input:
    path clone_ch // no parallelization

    output:
    path "*_productive_seq_clone-pass.tsv", emit: clone_split_ch // multiple files -> parall expected

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    FILENAME=\$(basename -- ${clone_ch}) # recover a file name without path
    cp -Lr ${clone_ch} "./TEMPO.tsv" # to have the hard file, not the symlink, because modifications will be performed inside
    chmod 777 TEMPO.tsv
    rm \$FILENAME # remove the initial file to avoid to send it into the channel
    cp -rp TEMPO.tsv "\$FILENAME" # -p for preserve permissions
    rm TEMPO.tsv

    # Determine the position of the clone_id column
    clone_id_col=\$(Rscript -e '
        args <- commandArgs(trailingOnly = TRUE)
        filename <- args[1]

        # Read the header of the file
        header <- read.table(filename, sep = "\\t", header = TRUE, nrows = 1)

        # Find the position of the clone_id column
        clone_id_col <- which(names(header) == "clone_id")

        if (length(clone_id_col) == 0) {
          cat("ERROR: clone_id column not found in the input file\\n")
          quit(status = 1)
        } else {
          cat(clone_id_col, "\\n")
        }
    ' \$FILENAME)

    if [ -z "\$clone_id_col" ]; then
        echo "ERROR IN THE SPLIT_BY_CLONES PROCESS: clone_id column not found in the input file"
        exit 1
    fi

    awk -v var1=\$FILENAME -v col=\$clone_id_col -F "\\t" '{
        print \$col
        if(NR == 1){
            header=\$0 ; next
        }else{
            clone_id=\$col
            if(system( "[ -f " \$col"_"var1 " ] " ) > 0){ # test if a file exists
                print header > \$col"_"var1
            }
            print \$0 > \$col"_"var1
        }
    }' \$FILENAME
    """
}



// This process adds the germline sequences for each v,d,j gene specified in the germline_(v|d|j)_call columns present in the tsv file
// Inputs:
//      - closest_ch: tsv file containing at least germline_(v|d|j)_call column (otherwise error), sequences already regrouped in clonal groups
//      - igblast_organism: value specified in nextflow.config
//      - igblast_variable_ref_files: also specified in nextflow.config, contains germline ref sequences of IMGT database
//                                     allele names already contained in the tsv will be used to search in thos files for the corresponding sequence
// Outputs:
//      - add_germ_ch: same tsv file as closest_ch, with the addition of the new germline_._seq columns
process AddGermlineSequences{
    label 'immcantation'
    cache 'true'

    input:
    path germline_coords_ch // parallelization expected (by clonal groups)
    val igblast_organism
    val igblast_variable_ref_files

    output:
    path "*_germ-seq.tsv", emit: add_germ_ch
    path "*.log", emit: add_germ_log_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    FILENAME=\$(basename -- ${germline_coords_ch}) # recover a file name without path
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a AddGermlineSequences.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a AddGermlineSequences.log

    REPO_PATH="/usr/local/share/germlines/imgt/${igblast_organism}/vdj" # path where the imgt_human_IGHV.fasta, imgt_human_IGHD.fasta and imgt_human_IGHJ.fasta files are in the docker container
    VDJ_FILES=\$(awk -v var1="${igblast_variable_ref_files}" -v var2="\${REPO_PATH}" 'BEGIN{ORS=" " ; split(var1, array1, " ") ; for (key in array1) {print var2"/"array1[key]}}')
    GermlineSequences.py -i \$FILENAME -r \${VDJ_FILES} |& tee -a GermlineSequences.log
    """
}




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
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a TranslateGermline.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a TranslateGermline.log
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
            cat(paste0("\\nWARNING: THE clonal_germline_sequence_no_gaps COLUMN CONTAINS ", length_no_gaps, " CHARACTERS WHEN GAPS ARE REMOVED, WHICH IS NOT A MULTIPLE OF 3. \\n"), file = "TranslateGermline.log", append = TRUE)
        }
        germ_dna <- Biostrings::DNAString(germ_nuc[1])
        # Catch a warning in the log file if raised
        withCallingHandlers(
            expr = {
                germ_aa <- Biostrings::translate(germ_dna, if.fuzzy.codon="X")
            },
            warning = function(w) {
                cat("WARNING OF Biostrings::translate FUNCTION: ", conditionMessage(w), "\\n", file = "TranslateGermline.log", append = TRUE)
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
    '|& tee -a TranslateGermline.log
    """

}


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


process Clone_id_count {
    label 'immcantation'
    publishDir path: "${out_path}/tsv", mode: 'copy', pattern: "clone_id_count.tsv", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "clone_id_count.log", overwrite: false
    cache 'true'

    input:
    path clone_assigned_seq_ch // no parallelization


    output:
    path "clone_id_count.tsv"
    path "clone_id_count.log"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    Rscript -e '
        obs <- read.table("./clone_assigned_seq.tsv", sep = "\\t", header = TRUE)
        clone_id_count <- as.data.frame(table(obs\$clone_id))
        names(clone_id_count) <- c("clone_id", "count")
        write.table(clone_id_count, file = "clone_id_count.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\\t")
    ' |& tee -a clone_id_count.log
    """
}


// Makes a fasta file with several sequences based on the sequences present in the tsv input in BOTH nucleotidic format and amino-acidic format
// The tsv input is expected to already only contain sequences of one same clonal group. For each clonal group, both nuc and aa fastas will be emitted paired up.
process Tsv2fasta {
    label 'r_ig_clustering'
    cache 'true'
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "{*_nuc/*.fasta}", overwrite: false
    publishDir path: "${out_path}/fasta", mode: 'copy', pattern: "{*_aa/*.fasta}", overwrite: false

    input:
    tuple path(all_files_ch), val(seq_kind) // parallelization expected (by clonal groups over align_clone_nb sequences)
    val align_seq
    val clone_germline_kind
    val align_clone_nb
    path cute_file

    output:
    tuple path("*_nuc/*.fasta"), path("*_aa/*.fasta"), val(seq_kind), emit: fasta_align_ch
    path "Tsv2fasta.log", emit: tsv2fasta_log_ch
    path "warning.txt", emit: warning_ch, optional: true

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    FILENAME=\$(basename -- ${all_files_ch}) # recover a file name without path
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a Tsv2fasta.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a Tsv2fasta.log
    Tsv2fasta.R \
    "${all_files_ch}" \
    "sequence_id" \
    "${align_seq}" \
    "${clone_germline_kind}" \
    "${align_clone_nb}" \
    "${cute_file}" \
    "${seq_kind}" \
    "Tsv2fasta.log"
    """
}


// Align the amino-acidic sequences that are already in fasta files (grouped by clonal groups)
process Abalign_align_aa {
    label 'abalign'
    
    input:
    tuple path(fasta_nuc), path(fasta_aa), val(seq_kind) // parallelization expected (by clonal groups over align_clone_nb sequences)
    val igblast_organism
    val igblast_heavy_chain
    val align_abalign_options
    
    output:
    tuple path(fasta_nuc), path(fasta_aa), path("*_align_aa.fasta"), val(seq_kind), emit: aligned_aa_ch
    path "Abalign_align_aa.log", emit: abalign_align_aa_log_ch
    
    script:
    if( ! (igblast_heavy_chain == "TRUE" || igblast_heavy_chain == "FALSE") ){
        error "\n\n========\n\nERROR IN Abalign_align_aa PROCESS\n\nINVALID heavy_chain PARAMETER:\n${heavy_chain}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    parms="-al"
    if(igblast_heavy_chain == "TRUE"){parms="-ah"}
    // Choose the species parameter for abalign
    switch (igblast_organism) {
        case "mouse":
            species = "MU"
            break
        case "human":
            species = "HS"
            break
        case "rabbit":
            species = "OC"
            break
        case "rat":
            species = "MM"
            break
        case "rhesus_monkey":
            species = "RM"
            break
        default:
            error "\n\n========\n\nERROR IN Abalign_align_aa PROCESS\n\nINVALID igblast_organism PARAMETER:\n${igblast_organism}\nMUST BE EITHER \"mouse\" OR \"human\" OR \"rabbit\" OR \"rat\" OR \"rhesus_monkey\"\n\n========\n\n"
    }
    """
    #!/bin/bash -ue
    set -o pipefail
    FILENAME=\$(basename -- ${fasta_aa}) # recover a file name without path
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a Abalign_align_aa.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a Abalign_align_aa.log
    /bin/Abalign_V2_Linux_Term/Abalign ${align_abalign_options} -i ${fasta_aa} ${parms} ${fasta_aa.baseName}_align_aa.fasta -sp ${species}. |& tee -a Abalign_align_aa.log || true
    # -g   (IMGT Numbering scheme, default)
    # -lc length.txt -lg 1 # single length computed spanning the indicated regions. For instance, -lc length.txt -lg 1,2,3,4,5,6,7 returns a single length
    """
}


// Abalign puts fasta headers in all caps. next script is meant to put those headers back to how they originally were
process Abalign_rename {
    label 'r_ext'
    publishDir path: "${out_path}/alignments/aa", mode: 'copy', pattern: "{*_aligned_aa.fasta}", overwrite: false
    
    input:
    tuple path(fasta_nuc), path(fasta_aa), path(fasta_aa_align), val(seq_kind) // parallelization expected (by clonal groups over align_clone_nb sequences)
    
    output:
    tuple path(fasta_nuc), path("*_aligned_aa.fasta"), val(seq_kind), emit: renamed_aligned_aa_ch
    path "*_failed_abalign_align.tsv", emit: failed_abalign_align_ch
    path "Abalign_rename.log", emit: abalign_rename_log_ch
    
    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    Rscript -e '
        options(warning.length = 8170)
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
            print(paste0("\\n\\nWARNING: ALIGNMENT FAILED FOR ", tempo_txt, "\\n\\n"))
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
        writeLines(df_align, con = "${fasta_aa.baseName}_aligned_aa_tempo.fasta")
    ' ${fasta_aa} ${fasta_aa.baseName}_align_aa.fasta |& tee -a Abalign_rename.log
    awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' ${fasta_aa.baseName}_aligned_aa_tempo.fasta > ${fasta_aa.baseName}_aligned_aa.fasta # remove \\n in seq
    """
}


// This process aligns the nucleotidic fasta files by transfering amino-acidic alignments onto the nucleotidic ones
process Abalign_align_nuc {
    label 'goalign'
    publishDir path: "${out_path}/alignments/nuc", mode: 'copy', pattern: "{*_aligned_nuc.fasta}", overwrite: false

    input:
    tuple path(fasta_nuc), path(aligned_aa), val(seq_kind) // parallelization expected (by clonal groups over align_clone_nb sequences)

    output:
    tuple path("*_aligned_nuc.fasta"), path(aligned_aa), val(seq_kind), emit: aligned_all_ch
    path "Abalign_align_nuc.log", emit: abalign_align_nuc_log_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    FILENAME=\$(basename -- ${fasta_nuc}) # recover a file name without path
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a Abalign_align_nuc.log

    goalign codonalign -i ${aligned_aa} -f ${fasta_nuc} -o ${fasta_nuc.baseName}_aligned_nuc_tempo.fasta |& tee -a Abalign_align_nuc.log
    awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' ${fasta_nuc.baseName}_aligned_nuc_tempo.fasta > ${fasta_nuc.baseName}_aligned_nuc.fasta # remove \\n in seq
    """
}

process  Mafft_align {
    label 'mafft'
    publishDir path: "${out_path}/alignments/aa", mode: 'copy', pattern: "{*_aligned_aa.fasta}", overwrite: false
    publishDir path: "${out_path}/alignments/nuc", mode: 'copy', pattern: "{*_aligned_nuc.fasta}", overwrite: false

    input:
    tuple path(fasta_nuc), path(fasta_aa), val(seq_kind) // parallelization expected (by clonal groups over align_clone_nb sequences)
    val align_mafft_all_options
    val align_mafft_clonal_options

    output:
    tuple path("*_aligned_nuc.fasta"), path("*_aligned_aa.fasta"), val(seq_kind), emit: aligned_all_ch
    path "Mafft_align.log", emit: mafft_align_log_ch

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
        PARAM=${align_mafft_all_options}
    elif [[ "${seq_kind}" == "CLONE" ]] ; then
        PARAM=${align_mafft_clone_options}
    fi

    mafft \$PARAM ${fasta_nuc} | awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' > ${fasta_nuc.baseName}_aligned_nuc_tempo.fasta
    awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' ${fasta_nuc.baseName}_aligned_nuc_tempo.fasta > ${fasta_nuc.baseName}_aligned_nuc.fasta # remove \\n in seq
    mafft \$PARAM ${fasta_aa} | awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' > ${fasta_aa.baseName}_aligned_aa_tempo.fasta
    awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' ${fasta_aa.baseName}_aligned_aa_tempo.fasta > ${fasta_aa.baseName}_aligned_aa.fasta # remove \\n in seq
    echo "" > Mafft_align.log
    """
}




process Tree {
    publishDir path: "${out_path}/phylo/aa", mode: 'copy', pattern: "{*_aligned_aa.fasta.treefile}", overwrite: false
    publishDir path: "${out_path}/phylo/nuc", mode: 'copy', pattern: "{*_aligned_nuc.fasta.treefile}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{*.log}", overwrite: false

    label 'iqtree'

    input:
    tuple path(fasta_nuc_alignments), path(fasta_aa_alignments)  // parallelization expected (by clonal groups over align_clone_nb sequences)
    path phylo_tree_model_file

    output:
    path "*_aligned_aa.fasta.treefile", emit: tree_aa_ch
    path "*_aligned_nuc.fasta.treefile", emit: tree_nuc_ch
    path "*.log", emit: tree_log_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    iqtree -nt ${task.cpus} -s ${fasta_aa_alignments} -m ${phylo_tree_model_file}+I+R6 --seed 123456789
    iqtree -nt ${task.cpus} -s ${fasta_nuc_alignments} -m ${phylo_tree_model_file}+GTR+I+R6 --seed 123456789
    """
}

process Meta2ITOL  {
    label 'tabletoitol'

    input:
    path meta_file // no parallelization
    val meta_seq_names

    output:
    path "iTOL*", emit: itol_out_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    if [ "${meta_file}" == "NULL" ]; then
        echo "No metadata provided; skipping iTOL annotation" > skip.iTOL.txt
        touch iTOL.empty # create the empty file iTOL.empty
        exit 0
    else
        table2itol.R -i ${meta_seq_names} ${meta_file}
        # create three files: iTOL_colorstrip-Isotype.txt, iTOL_colorstrip-Name.txt, iTOL_gradient-KD.txt
    fi
    """
}


process ITOL{
    label 'gotree'

    input:
    path tree
    path itol_files
    val phylo_tree_itolkey

    output:
    path "*_aligned_aa.fasta_itol_url.txt", emit: itol_aa_ch, optional: true
    path "*_aligned_nuc.fasta_itol_url.txt", emit: itol_nuc_ch, optional: true

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    # Check if itol_files is made of iTOL.empty # -s means exists and non-empty
    if [ ! -s "iTOL.empty" ]; then 
        gotree upload itol --project gotree_uploads --user-id $phylo_tree_itolkey -i $tree > ${tree.baseName}_itol_url.txt 2>&1
    else
        # add metadata
        gotree upload itol --project gotree_uploads --user-id $phylo_tree_itolkey -i $tree $itol_files > ${tree.baseName}_itol_url.txt 2>&1 
    fi
    """
}


// Creates donut plots for the constant and variable region ; grouped by same allele or genes
// Inputs:
//      - kind: kind of donut plot to display. can be "all", "annotated" or "tree"
//      - data: tsv file containing the sequences to plot 
//      - col: columns in the data file to take into account when plotting the donut (can be c_call if the donut should be grouped by same constant region alleles for instance)
//      - donut_*: parameters for the display of the donut. These are entered in the nextflow.config file
//      - cute_file: file containing R functions used in the donut.R script
//      - igblast_variable_ref_files & igblast_constant_ref_files: only used to be certain that all studied loci are displayed in the donut's title (if we included IGL and IGK in the study but only one was found, both loci still need to be in the title)
process donut {
    label 'r_ext'
    //publishDir path: "${out_path}", mode: 'copy', pattern: "{*.tsv}", overwrite: false
    //publishDir path: "${out_path}", mode: 'copy', pattern: "{*.pdf}", overwrite: false
    publishDir path: "${out_path}/figures/png", mode: 'copy', pattern: "{*.png}", overwrite: false
    publishDir path: "${out_path}/figures/svg", mode: 'copy', pattern: "{*.svg}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{*.log}", overwrite: false
    cache 'true'

    input:
    tuple val(kind), path(data), val(col) // 12 parallelization expected
    val donut_palette
    val donut_hole_size
    val donut_hole_text
    val donut_hole_text_size
    val donut_border_color
    val donut_border_size
    val donut_annotation_distance
    val donut_annotation_size
    val donut_annotation_force
    val donut_annotation_force_pull
    val donut_legend_width
    val donut_legend_text_size
    val donut_legend_box_size
    val donut_legend_box_space
    val donut_legend_limit
    path cute_file
    path igblast_data_check_ch

    output:
    path "*.tsv", emit: donut_tsv_ch, optional: true
    path "*.pdf", emit: donut_pdf_ch, optional: true
    path "*.png", emit: donuts_png
    path "*.svg"
    path "*.log"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    FILENAME=\$(basename -- ${data}) # recover a file name without path
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\nKIND: ${kind}\\nCOL: ${col}\\n\\n################################\\n\\n" |& tee -a ${kind}_donut.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a ${kind}_donut.log
    donut.R \
"${data}" \
"${kind}" \
"${col}" \
"${donut_palette}" \
"${donut_hole_size}" \
"${donut_hole_text}" \
"${donut_hole_text_size}" \
"${donut_border_color}" \
"${donut_border_size}" \
"${donut_annotation_distance}" \
"${donut_annotation_size}" \
"${donut_annotation_force}" \
"${donut_annotation_force_pull}" \
"${donut_legend_width}" \
"${donut_legend_text_size}" \
"${donut_legend_box_size}" \
"${donut_legend_box_space}" \
"${donut_legend_limit}" \
"${cute_file}" \
"${igblast_data_check_ch}" \
"${kind}_donut.log"
    """
}

process donut_assembly {
    label 'r_ig_clustering'
    publishDir path: "${out_path}/pdf", mode: 'copy', pattern: "{donuts.pdf}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{donut_assembly.log}", overwrite: false
    cache 'true'

    input:
    path donut_pdf_ch2

    output:
    path "donuts.pdf", emit: donut_assembly_ch
    path "donut_assembly.log"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    Rscript -e '
        files <- list.files(path = ".", pattern = ".pdf\$")

        sorted_files <- files[order(
            !grepl("allele", files),  # Alleles will be displayed before genes
            grepl("gene", files),     
            !grepl("vj_", files)     # Amongst alleles and genes, vj will be displayed before c
        )]

        qpdf::pdf_combine(input = sorted_files, output = "./donuts.pdf")
    ' |& tee -a donut_assembly.log
    """
}

    /*
// alakazam::readChangeoDb() ; dowser::formatClones() ; dowser::getTrees()
process get_germ_tree {
    label 'immcantation_10cpu'
    publishDir path: "${out_path}/RData", mode: 'copy', pattern: "{*_get_germ_tree_cloneID.RData}", overwrite: false
    cache 'true'

    input:
    path mutation_load_ch // parallelization expected
    path meta_file // just to determine if metadata have been provided (TRUE means NULL) meta_file_ch not required here
    path cute_file
    val align_clone_nb
    val germ_tree_duplicate_seq
    val igphylm_exe_path // warning: here val and not path because we do not want the igphyml file to be imported in the work dir

    output:
    path "*_get_germ_tree_cloneID.RData", emit: rdata_germ_tree_ch, optional: true
    path "germ_tree_dismissed_seq.tsv", emit: no_germ_tree_ch
    path "seq_for_germ_tree.tsv", emit: germ_tree_ch
    path "germ_tree_dismissed_clone_id.tsv", emit: no_cloneID_ch
    path "germ_tree_clone_id.tsv", emit: cloneID_ch
    path "get_germ_tree.log", emit: get_germ_tree_log_ch
    //path "HLP10_germ_tree_parameters.tsv"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    if [[ ! -s ${mutation_load_ch} ]]; then
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nEMPTY ${mutation_load_ch} FILE AS INPUT OF THE mutation_load PROCESS\\nCHECK THE mutation_load.log IN THE report FOLDER INSIDE THE OUTPUT FOLDER\\n\\n========\\n\\n"
        exit 1
    fi
    FILENAME=\$(basename -- ${mutation_load_ch}) # recover a file name without path
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a get_germ_tree.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a get_germ_tree.log
    get_germ_tree.R \
"${mutation_load_ch}" \
"${meta_file}" \
"${align_clone_nb}" \
"${germ_tree_duplicate_seq}" \
"${igphylm_exe_path}" \
"${cute_file}" \
"get_germ_tree.log"
    """
}

process germ_tree_vizu {
    label 'r_ig_clustering'
    publishDir path: "${out_path}/pdf", mode: 'copy', pattern: "{germ_tree.pdf}", overwrite: false
    publishDir path: "${out_path}/pdf", mode: 'copy', pattern: "{germ_no_tree.pdf}", overwrite: false
    publishDir path: "${out_path}/figures/png", mode: 'copy', pattern: "{*.png}", overwrite: false
    publishDir path: "${out_path}/figures/svg", mode: 'copy', pattern: "{*.svg}", overwrite: false
    publishDir path: "${out_path}/RData", mode: 'copy', pattern: "{all_trees.RData}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{germ_tree_vizu.log}", overwrite: false
    cache 'true'

    input:
    path rdata_germ_tree_ch2 // no more parallelization
    val germ_tree_kind
    val align_clone_nb
    val germ_tree_duplicate_seq
    val germ_tree_leaf_color
    val germ_tree_leaf_shape
    val germ_tree_leaf_size
    val germ_tree_label_size
    val germ_tree_label_hjust
    val germ_tree_label_rigth
    val germ_tree_label_outside
    val germ_tree_right_margin
    val germ_tree_legend
    path data_assembly_ch
    path meta_file
    val meta_legend
    path cute_file

    output:
    path "*.RData", optional: true
    path "germ_tree.pdf"
    path "germ_no_tree.pdf"
    path "*.png", emit: germ_tree_vizu_ch // png plot (but sometimes empty) systematically returned
    path "*.svg"
    path "*germ_tree_dup_seq_not_displayed.tsv", emit: germ_tree_dup_seq_not_displayed_ch
    path "germ_tree_vizu.log"
    //path "HLP10_germ_tree_parameters.tsv"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    germ_tree_vizu.R \
"${germ_tree_kind}" \
"${align_clone_nb}" \
"${germ_tree_duplicate_seq}" \
"${germ_tree_leaf_color}" \
"${germ_tree_leaf_shape}" \
"${germ_tree_leaf_size}" \
"${germ_tree_label_size}" \
"${germ_tree_label_hjust}" \
"${germ_tree_label_rigth}" \
"${germ_tree_label_outside}" \
"${germ_tree_right_margin}" \
"${germ_tree_legend}" \
"${data_assembly_ch}" \
"${meta_file}" \
"${meta_legend}" \
"${cute_file}" \
"germ_tree_vizu.log"
    """
}
    */

// The render function only creates the html file with the rmardown
// Inputs:
//      template_rmd: rmardown template file used to create the html (file path is defined in config.nextflow)
//      nb_input: number of fasta sequences in initial input
//      nb_seq_aligned: number of sequences that igblast could analyse
//      nb_seq_not_aligned: number of sequences that igblast could not alayse
//      donuts_png: provide links to the donut images needed in the rmd inside the work folder
//      repertoire_png: provide links to the repertoire images needed in the rmd inside the work folder
//      repertoire_constant_ch: names of the constant gene repertoire files to be displayed
//      repertoire_vj_ch: names of the variable gene repertoire files to be displayed
//      itol_subscription: nextflow.config parameter to know if user has paid the subscription to itol automated visualization of trees, process ITOL is only executed if TRUE
//      heavy_chain: to know if the analyzed data is VL or VH, because "Amino acid sequences phylogeny" section in html report is only displayed for VH
// Outputs:
//      "report.html": finalized html report for a specific run
//      "print_report.log": will contain any error or warning messages produced by rmardown::render
process print_report{
    label 'r_ig_clustering'

    publishDir path: "${out_path}", mode: 'copy', pattern: "{report.html}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{alignments_viz.html}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{print_report.log}", overwrite: false
    cache 'false'

    input:
    file template_rmd
    file alignments_viz_rmd
    val nb_input
    val nb_seq_aligned
    val nb_seq_not_aligned
    val nb_productive
    val nb_unproductive
    val nb_clone_assigned
    val nb_failed_clone
    path donuts_png
    path reperoire_png
    val repertoire_constant_ch
    val repertoire_vj_ch
    val align_soft
    val itol_subscription

    output:
    file "report.html"
    file "alignments_viz.html"
    file "print_report.log"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    cp ${template_rmd} report_file.rmd
    cp ${alignments_viz_rmd} alignments_vizu.rmd
    cp -r "${out_path}/tsv" .
    cp -r "${out_path}/pdf" .
    cp -r "${projectDir}/bin/doc_images" .
    Rscript -e '

        # Find the constant and vj repertoires to be displayed in the html file (file names differ depending on light/heavy chain)
        constant_files <- as.character("${repertoire_constant_ch}")
        vj_files <- as.character("${repertoire_vj_ch}")
        cleaned_repertoire_constant <- gsub("^\\\\[\\\\[|\\\\]\\\\]\$", "", constant_files)
        cleaned_repertoire_vj <- gsub("^\\\\[\\\\[|\\\\]\\\\]\$", "", vj_files)
        constant_paths <- strsplit(cleaned_repertoire_constant, ",\\\\s*")[[1]]
        vj_paths <- strsplit(cleaned_repertoire_vj, ",\\\\s*")[[1]]
        constant_names <- basename(constant_paths)
        vj_names <- basename(vj_paths)
        # Verification that the resulting file names are as expected
        constant_rep <- constant_names[grepl("^IG.C_.*gene_non-zero\\\\.png\$", constant_names)]
        vj_rep <- vj_names[grepl("^rep_gene_IG.V_.*non-zero\\\\.png\$", vj_names)]
        if(length(constant_rep) == 0 || length(vj_rep) == 0){
            stop(paste0("\\n\\n========\\n\\nERROR IN print_report PROCESS\\n\\nTHE REPERTOIRE PNG FILES TO BE DISPLAYED WERE NOT FOUND\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n"), call. = FALSE)
        }

        rmarkdown::render(
        input = "report_file.rmd",
        output_file = "report.html",
        # list of the variables waiting to be replaced in the rmd file:
        params = list(nb_input = ${nb_input},
            nb_seq_aligned = ${nb_seq_aligned}, 
            nb_seq_not_aligned = ${nb_seq_not_aligned},
            nb_productive = ${nb_productive},
            nb_unproductive = ${nb_unproductive},
            nb_clone_assigned = ${nb_clone_assigned},
            nb_failed_clone = ${nb_failed_clone},
            constant_rep = constant_rep,
            vj_rep = vj_rep,
            align_soft = "${align_soft}",
            itol_subscription = ${itol_subscription}
        ),
        # output_dir = ".",
        # intermediates_dir = "./",
        # knit_root_dir = "./",
        run_pandoc = TRUE,
        quiet = TRUE,
        clean = TRUE
        )

        # dir.create("reports", showWarnings = FALSE, recursive = TRUE)
        rmarkdown::render(
        input = "alignments_vizu.rmd",
        output_file = "alignments_viz.html", 
        run_pandoc = TRUE,
        quiet = TRUE,
        clean = TRUE
        )
    ' |& tee -a print_report.log
    """
}


// Save the config file and the log file for a specific run
process backup {
    label 'bash'
    publishDir "${out_path}/reports", mode: 'copy', overwrite: false // since I am in mode copy, all the output files will be copied into the publishDir. See \\wsl$\Ubuntu-20.04\home\gael\work\aa\a0e9a739acae026fb205bc3fc21f9b
    cache 'false'

    input:
    path config_file
    path log_file

    output:
    path "${config_file}" // warning message if we use path config_file
    path "${log_file}" // warning message if we use path log_file
    path "Log_info.txt"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    echo -e "full .nextflow.log is in: ${launchDir}\nThe one in the result folder is not complete (miss the end)" > Log_info.txt
    """
}

//////// End Processes


//////// Modules

include { Igblast_query as Igblast_query } from './modules/igblast_query'
include { Igblast_germline_coords as Igblast_germline_coords  } from './modules/igblast_germline'
include { Closest_germline as Closest_germline  } from './modules/closest_germline'
include { Gff as GffNuc } from './modules/gff'
include { Gff as GffAa  } from './modules/gff'
include { PrintAlignment as PrintAlignmentNuc } from './modules/print_alignment'
include { PrintAlignment as PrintAlignmentAa  } from './modules/print_alignment'

//////// end Modules


//////// Workflow


workflow {

    print("\n\nINITIATION TIME: ${workflow.start}")

    //////// Options of nextflow run

    // --modules (it is just for the process workflowParam)
    params.modules = "" // if --module is used, this default value will be overridden
    // end --modules (it is just for the process workflowParam)

    //////// end Options of nextflow run


    //////// Variables

    modules = params.modules // remove the dot -> can be used in bash scripts
    config_file = workflow.configFiles[0] // better to use this than config_file = file("${projectDir}/nextflow.config") because the latter is not good if -c option of nextflow run is used
    log_file = file("${launchDir}/.nextflow.log")

    //////// end Variables


    //////// Checks
    //// check of the bin folder
    tested_files_bin = ["circos.R", "circos_data_prep.R", "cute_little_R_functions_v12.8.R", "defineGroups.pl", "donut.R", "fields_not_kept.txt", "germ_tree_vizu.R", "GermlineSequences.py", "get_germ_tree.R", "histogram.R", "parse_coordinates.R", "repertoire.R", "repertoire_profiler_template.rmd", "trimtranslate.sh", "Tsv2fasta.R", "Gff.R"]
    for(i1 in tested_files_bin){
        if( ! (file("${projectDir}/bin/${i1}").exists()) ){
            error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nTHE ${i1} FILE MUST BE PRESENT IN THE ./bin FOLDER, WHERE THE main.nf file IS PRESENT\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
        }
    }
    //// end check of the bin folder

    //// check of the modules folder
    tested_files_modules = ["print_alignment.nf", "igblast_query.nf", "igblast_germline.nf", "closest_germline.nf", "gff.nf"]
    for(i1 in tested_files_modules){
        if( ! (file("${projectDir}/modules/${i1}").exists()) ){
            error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nTHE ${i1} FILE MUST BE PRESENT IN THE ./modules FOLDER, WHERE THE main.nf file IS PRESENT\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
        }
    }
    //// end check of the modules folder

    //// check of config file parameters
    // Data
    if( ! (sample_path in String || sample_path in GString) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN nextflow.config FILE:\n${sample_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (file(sample_path).exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN nextflow.config FILE (DOES NOT EXIST): ${sample_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    if( ! (meta_path in String || meta_path in GString) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID meta_path PARAMETER IN nextflow.config FILE:\n${meta_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if(meta_path != "NULL"){
        if( ! (file(meta_path).exists()) ){
            error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID meta_path PARAMETER IN nextflow.config FILE (DOES NOT EXIST): ${meta_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
        }
    }
    if( ! (meta_seq_names in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID meta_seq_names PARAMETER IN nextflow.config FILE:\n${meta_seq_names}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (meta_name_replacement in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID meta_name_replacement PARAMETER IN nextflow.config FILE:\n${meta_name_replacement}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }

    // Ig annotation
    if( ! (igblast_organism in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_organism PARAMETER IN nextflow.config FILE:\n${igblast_organism}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (igblast_organism =~ /^(mouse|human|rabbit|rat|rhesus_monkey)$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_organism PARAMETER IN nextflow.config FILE:\n${igblast_organism}\nMUST BE EITHER \"mouse\", \"human\", \"rabbit\", \"rat\" OR \"rhesus_monkey\"\n\n========\n\n"
    }else if( ! (igblast_organism =~ /^(mouse|human)$/)){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_organism PARAMETER IN nextflow.config FILE: ${igblast_organism}\n\nTHE repertoire PROCESS CURRENTLY ONLY SUPPORTS mouse AND human SPECIES\nTHEREFORE igblast_organism MUST BE EITHER \"mouse\" OR \"human\"\n\n========\n\n"
    }
    if( ! (igblast_loci in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_loci PARAMETER IN nextflow.config FILE:\n${igblast_loci}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (igblast_loci == "ig" || igblast_loci == "tr") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_loci PARAMETER IN nextflow.config FILE:\n${igblast_loci}\nMUST BE EITHER \"ig\" OR \"tr\"\n\n========\n\n"
    }
    if( ! (igblast_heavy_chain in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_heavy_chain PARAMETER IN nextflow.config FILE:\n${igblast_heavy_chain}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (igblast_heavy_chain == "TRUE" || igblast_heavy_chain == "FALSE" || igblast_heavy_chain == "true" ||igblast_heavy_chain == "false") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_heavy_chain PARAMETER IN nextflow.config FILE:\n${igblast_heavy_chain}\nMUST BE EITHER \"TRUE\" OR \"FALSE\" (OR \"true\" OR \"false\") \n\n========\n\n"
    }
    if( ! (igblast_lambda_chain in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_lambda_chain PARAMETER IN nextflow.config FILE:\n${igblast_lambda_chain}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (igblast_lambda_chain == "TRUE" || igblast_lambda_chain == "FALSE" || igblast_lambda_chain == "true" ||igblast_lambda_chain == "false") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_lambda_chain PARAMETER IN nextflow.config FILE:\n${igblast_lambda_chain}\nMUST BE EITHER \"TRUE\" OR \"FALSE\" (OR \"true\" OR \"false\") \n\n========\n\n"
    }
    if( ! (igblast_kappa_chain in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_kappa_chain PARAMETER IN nextflow.config FILE:\n${igblast_kappa_chain}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (igblast_kappa_chain == "TRUE" || igblast_kappa_chain == "FALSE" || igblast_kappa_chain == "true" ||igblast_kappa_chain == "false") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_kappa_chain PARAMETER IN nextflow.config FILE:\n${igblast_kappa_chain}\nMUST BE EITHER \"TRUE\" OR \"FALSE\" (OR \"true\" OR \"false\") \n\n========\n\n"
    }
    // Checking of studied chain coherence
    heavy = igblast_heavy_chain == "TRUE" || igblast_heavy_chain == "true"
    lambda = igblast_lambda_chain == "TRUE" || igblast_lambda_chain == "true"
    kappa  = igblast_kappa_chain == "TRUE" || igblast_kappa_chain == "true"
    if (heavy && (lambda || kappa)) {
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID PARAMETERS IN nextflow.config FILE:\nHEAVY (igblast_heavy_chain) AND LIGHT (igblast_lambda_chain, igblast_kappa_chain) CHAIN LOCI CANNOT BOTH BE TRUE.\nHERE ARE THEIR CURRENT VALUES: \nigblast_heavy_chain: ${igblast_heavy_chain}\nigblast_lambda_chain: ${igblast_lambda_chain}\nigblast_kappa_chain: ${igblast_kappa_chain}\n\n========\n\n"
    }
    if (!heavy && !lambda && !kappa) {
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID PARAMETERS IN nextflow.config FILE:\nAT LEAST ONE OF igblast_heavy_chain, igblast_lambda_chain OR igblast_kappa_chain MUST BE TRUE\nHERE ARE THEIR CURRENT VALUES: \nigblast_heavy_chain: ${igblast_heavy_chain}\nigblast_lambda_chain: ${igblast_lambda_chain}\nigblast_kappa_chain: ${igblast_kappa_chain}\n\n========\n\n"
    }
    if ([heavy, lambda, kappa].count { it } > 2) {
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID PARAMETERS IN nextflow.config FILE:\nONLY ONE OF igblast_heavy_chain AND LIGHT CHAINS (igblast_lambda_chain, igblast_kappa_chain) CAN BE TRUE AT A TIME\nHERE ARE THEIR CURRENT VALUES: \nigblast_heavy_chain: ${igblast_heavy_chain}\nigblast_lambda_chain: ${igblast_lambda_chain}\nigblast_kappa_chain: ${igblast_kappa_chain}\n\n========\n\n"
    }

    // Clonal groups (clustering) and mutation load
    if( ! (clone_strategy in String) ){
                error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_strategy PARAMETER IN nextflow.config FILE:\n${clone_strategy}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (clone_strategy == "first" || clone_strategy == "set") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_strategy PARAMETER IN nextflow.config FILE:\n${clone_strategy}\nMUST BE EITHER \"first\" OR \"set\"\n\n========\n\n"
    }
    if( ! (clone_model in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_model PARAMETER IN nextflow.config FILE:\n${clone_model}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (clone_model == "ham" || clone_model == "aa" || clone_model == "hh_s1f" || clone_model == "hh_s5f" || clone_model == "mk_rs1nf" || clone_model == "mk_rs5nf" || clone_model == "m1n_compat" || clone_model == "hs1f_compat") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_model PARAMETER IN nextflow.config FILE:\n${clone_model}\nMUST BE EITHER \"ham\", \"aa\", \"hh_s1f\", \"hh_s5f\", \"mk_rs1nf\", \"mk_rs5nf\", \"m1n_compat\", \"hs1f_compat\"\n\n========\n\n"
    }
    if( ! (clone_normalize in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_normalize PARAMETER IN nextflow.config FILE:\n${clone_normalize}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (clone_normalize == "len" || clone_normalize == "none") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_normalize PARAMETER IN nextflow.config FILE:\n${clone_normalize}\nMUST BE EITHER \"len\" OR \"none\"\n\n========\n\n"
    }
    if( ! (clone_distance in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_distance PARAMETER IN nextflow.config FILE:\n${clone_distance}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (clone_distance =~ /^((1)|(0)|(0\.[0-9]*))$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_distance PARAMETER IN nextflow.config FILE:\n${clone_distance}\nMUST BE A POSITIVE PROPORTION VALUE\n\n========\n\n"
    }
    if( ! (clone_germline_kind in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_germline_kind PARAMETER IN nextflow.config FILE:\n${clone_germline_kind}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (clone_germline_kind =~ /^(dmask|full|vonly)$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_germline_kind PARAMETER IN nextflow.config FILE:\n${clone_germline_kind}\nMUST BE EITHER \"dmask\", \"full\", \"vonly\".\n\n========\n\n"
    }

    if( ! (clone_mut_obs_seq in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_mut_obs_seq PARAMETER IN nextflow.config FILE:\n${clone_mut_obs_seq}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (clone_mut_germ_seq in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_mut_germ_seq PARAMETER IN nextflow.config FILE:\n${clone_mut_germ_seq}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (clone_mut_regionDefinition in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_mut_regionDefinition PARAMETER IN nextflow.config FILE:\n${clone_mut_regionDefinition}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (clone_mut_regionDefinition =~ /^(NULL|IMGT_V|IMGT_V_BY_CODONS|IMGT_V_BY_REGIONS|IMGT_V_BY_SEGMENTS|IMGT_VDJ|IMGT_VDJ_BY_REGIONS)$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_mut_regionDefinition PARAMETER IN nextflow.config FILE:\n${clone_mut_regionDefinition}\nMUST BE EITHER \"NULL\", \"IMGT_V\", \"IMGT_V_BY_CODONS\", \"IMGT_V_BY_REGIONS\", \"IMGT_V_BY_SEGMENTS\", \"IMGT_VDJ\", \"IMGT_VDJ_BY_REGIONS\".\n\n========\n\n"
    }

    // Aligments
    if( ! (align_clone_nb in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID align_clone_nb PARAMETER IN nextflow.config FILE:\n${align_clone_nb}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ( ! (align_clone_nb =~/^\d+$/)) || align_clone_nb.toInteger() < 2 ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID align_clone_nb PARAMETER IN nextflow.config FILE:\n${align_clone_nb}\nMUST BE A POSITIVE INTEGER VALUE EQUAL OR GREATER TO 2\n\n========\n\n"
    }
    if( ! (align_soft in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID align_soft PARAMETER IN nextflow.config FILE:\n${align_soft}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (align_soft =~ /^(abalign|mafft)$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID align_soft PARAMETER IN nextflow.config FILE:\n${align_soft}\nMUST BE EITHER \"abalign\", OR \"mafft\".\n\n========\n\n"
    }
    if( ! (align_seq in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID align_seq PARAMETER IN nextflow.config FILE:\n${align_seq}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (align_seq =~ /^(query|igblast_full|trimmed|fwr1|fwr2|fwr3|fwr4|cdr1|cdr2|cdr3|junction|sequence_alignment|v_sequence_alignment|d_sequence_alignment|j_sequence_alignment|c_sequence_alignment|germline_alignment|v_germline_alignment|d_germline_alignment|j_germline_alignment|c_germline_alignment)$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID align_seq PARAMETER IN nextflow.config FILE:\n${align_seq}\nMUST BE EITHER \"query\", \"igblast_full\", \"trimmed\", \"fwr1\", \"fwr2\", \"fwr3\", \"fwr4\", \"cdr1\", \"cdr2\", \"cdr3\", \"junction\", \"sequence_alignment\", \"v_sequence_alignment\", \"d_sequence_alignment\", \"j_sequence_alignment\", \"c_sequence_alignment\", \"germline_alignment\", \"v_germline_alignment\", \"d_germline_alignment\", \"j_germline_alignment\", \"c_germline_alignment\".\n\n========\n\n"
    }
    if( ! (align_abalign_options in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID align_abalign_options PARAMETER IN nextflow.config FILE:\n${align_abalign_options}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (align_mafft_all_options in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID align_mafft_all_options PARAMETER IN nextflow.config FILE:\n${align_mafft_all_options}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (align_mafft_clonal_options in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID align_mafft_clonal_options PARAMETER IN nextflow.config FILE:\n${align_mafft_clonal_options}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (meta_legend in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID meta_legend PARAMETER IN nextflow.config FILE:\n${meta_legend}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (germ_tree_kind in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_kind PARAMETER IN nextflow.config FILE:\n${germ_tree_kind}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (germ_tree_duplicate_seq in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_duplicate_seq PARAMETER IN nextflow.config FILE:\n${germ_tree_duplicate_seq}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (germ_tree_duplicate_seq == "TRUE" || germ_tree_duplicate_seq == "FALSE") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_duplicate_seq PARAMETER IN nextflow.config FILE:\n${germ_tree_duplicate_seq}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if( ! (germ_tree_leaf_color in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_leaf_color PARAMETER IN nextflow.config FILE:\n${germ_tree_leaf_color}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (germ_tree_leaf_shape in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_leaf_shape PARAMETER IN nextflow.config FILE:\n${germ_tree_leaf_shape}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (germ_tree_leaf_shape =~  /^[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_leaf_shape PARAMETER IN nextflow.config FILE:\n${germ_tree_leaf_shape}\nMUST BE A POSITIVE INTEGER VALUE\n\n========\n\n"
    }
    if( ! (germ_tree_leaf_size in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_leaf_size PARAMETER IN nextflow.config FILE:\n${germ_tree_leaf_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (germ_tree_leaf_size =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_leaf_size PARAMETER IN nextflow.config FILE:\n${germ_tree_leaf_size}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (germ_tree_label_size in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_size PARAMETER IN nextflow.config FILE:\n${germ_tree_label_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (germ_tree_label_size =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_size PARAMETER IN nextflow.config FILE:\n${germ_tree_label_size}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (germ_tree_label_hjust in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_hjust PARAMETER IN nextflow.config FILE:\n${germ_tree_label_hjust}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (germ_tree_label_hjust =~  /^\-{0,1}[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_hjust PARAMETER IN nextflow.config FILE:\n${germ_tree_label_hjust}\nMUST BE A NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (germ_tree_label_rigth in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_rigth PARAMETER IN nextflow.config FILE:\n${germ_tree_label_rigth}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (germ_tree_label_rigth == "TRUE" || germ_tree_label_rigth == "FALSE") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_rigth PARAMETER IN nextflow.config FILE:\n${germ_tree_label_rigth}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if( ! (germ_tree_label_outside in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_outside PARAMETER IN nextflow.config FILE:\n${germ_tree_label_outside}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (germ_tree_label_outside == "TRUE" || germ_tree_label_outside == "FALSE") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_outside PARAMETER IN nextflow.config FILE:\n${germ_tree_label_outside}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if( ! (germ_tree_right_margin in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_right_margin PARAMETER IN nextflow.config FILE:\n${germ_tree_right_margin}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (germ_tree_right_margin =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_right_margin PARAMETER IN nextflow.config FILE:\n${germ_tree_right_margin}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (germ_tree_legend in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_legend PARAMETER IN nextflow.config FILE:\n${germ_tree_legend}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (germ_tree_legend == "TRUE" || germ_tree_legend == "FALSE") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_legend PARAMETER IN nextflow.config FILE:\n${germ_tree_legend}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }



    if( ! (donut_palette in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_palette PARAMETER IN nextflow.config FILE:\n${donut_palette}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (donut_hole_size in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_hole_size PARAMETER IN nextflow.config FILE:\n${donut_hole_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_hole_size =~  /^((1)|(0)|(0\.[0-9]*))$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_hole_size PARAMETER IN nextflow.config FILE:\n${donut_hole_size}\nMUST BE A POSITIVE PROPORTION VALUE\n\n========\n\n"
    }
    if( ! (donut_hole_text in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_hole_text PARAMETER IN nextflow.config FILE:\n${germ_tree_legend}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_hole_text == "TRUE" || donut_hole_text == "FALSE") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_hole_text PARAMETER IN nextflow.config FILE:\n${donut_hole_text}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if( ! (donut_hole_text_size in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_hole_text_size PARAMETER IN nextflow.config FILE:\n${donut_hole_text_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_hole_text_size =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_hole_text_size PARAMETER IN nextflow.config FILE:\n${donut_hole_text_size}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (donut_border_color in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_border_color PARAMETER IN nextflow.config FILE:\n${donut_border_color}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (donut_border_size in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_border_size PARAMETER IN nextflow.config FILE:\n${donut_border_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (donut_annotation_distance in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_distance PARAMETER IN nextflow.config FILE:\n${donut_annotation_distance}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_annotation_distance =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_distance PARAMETER IN nextflow.config FILE:\n${donut_annotation_distance}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (donut_annotation_size in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_size PARAMETER IN nextflow.config FILE:\n${donut_annotation_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_annotation_size =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_size PARAMETER IN nextflow.config FILE:\n${donut_annotation_size}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (donut_annotation_force in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_force PARAMETER IN nextflow.config FILE:\n${donut_annotation_force}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_annotation_force =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_force PARAMETER IN nextflow.config FILE:\n${donut_annotation_force}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (donut_annotation_force_pull in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_force_pull PARAMETER IN nextflow.config FILE:\n${donut_annotation_force_pull}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_annotation_force_pull =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_force_pull PARAMETER IN nextflow.config FILE:\n${donut_annotation_force_pull}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (donut_legend_width in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_width PARAMETER IN nextflow.config FILE:\n${donut_legend_width}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_legend_width =~  /^((1)|(0)|(0\.[0-9]*))$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_width PARAMETER IN nextflow.config FILE:\n${donut_legend_width}\nMUST BE A POSITIVE PROPORTION VALUE\n\n========\n\n"
    }
    if( ! (donut_legend_text_size in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_text_size PARAMETER IN nextflow.config FILE:\n${donut_legend_text_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_legend_text_size =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_text_size PARAMETER IN nextflow.config FILE:\n${donut_legend_text_size}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (donut_legend_box_size in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_box_size PARAMETER IN nextflow.config FILE:\n${donut_legend_box_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_legend_box_size =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_box_size PARAMETER IN nextflow.config FILE:\n${donut_legend_box_size}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (donut_legend_box_space in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_box_space PARAMETER IN nextflow.config FILE:\n${donut_legend_box_space}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_legend_box_space =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_box_space PARAMETER IN nextflow.config FILE:\n${donut_legend_box_space}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (donut_legend_limit in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_limit PARAMETER IN nextflow.config FILE:\n${donut_legend_limit}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_legend_width ==  "NULL") ){
        if( ! (donut_legend_width =~  /^((1)|(0)|(0\.[0-9]*))$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_limit PARAMETER IN nextflow.config FILE:\n${donut_legend_limit}\nMUST BE A POSITIVE PROPORTION VALUE IF NOT \"NULL\"\n\n========\n\n"
        }
    }
    if( ! (phylo_tree_model_path in String || phylo_tree_model_path in GString) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID phylo_tree_model_path PARAMETER IN nextflow.config FILE:\n${phylo_tree_model_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (file(phylo_tree_model_path).exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID phylo_tree_model_path PARAMETER IN nextflow.config FILE (DOES NOT EXIST): ${phylo_tree_model_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    if( ! (phylo_tree_itolkey in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID phylo_tree_itolkey PARAMETER IN nextflow.config FILE:\n${phylo_tree_itolkey}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (phylo_tree_itol_subscription in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID phylo_tree_itol_subscription PARAMETER IN nextflow.config FILE:\n${phylo_tree_itol_subscription}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (phylo_tree_itol_subscription == "TRUE" || phylo_tree_itol_subscription == "FALSE") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID phylo_tree_itol_subscription PARAMETER IN nextflow.config FILE:\n${phylo_tree_itol_subscription}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if( ! (cute_path in String || cute_path in GString) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID cute_path PARAMETER IN nextflow.config FILE:\n${cute_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (file(cute_path).exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID cute_path PARAMETER IN nextflow.config FILE (DOES NOT EXIST): ${cute_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    if( ! (igphylm_exe_path in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igphylm_exe_path PARAMETER IN nextflow.config FILE:\n${igphylm_exe_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    //// end check of config file parameters

    // below: those variable are already used in the config file. Thus, to late to check them. And not possible to check inside the config file
    // out_ini
    print("\n\nRESULT DIRECTORY: ${out_path}")
    //print("\n\nWARNING: PARAMETERS ALREADY INTERPRETED IN THE .config FILE:")
    //print("    system_exec: ${system_exec}")
    //print("    out_path: ${out_path_ini}")
    if("${system_exec}" == "slurm"){
        print("    queue: ${slurm_queue}")
        print("    qos: ${slurm_qos}")
    }
    if("${system_exec}" != "local"){
        print("    add_options: ${add_options}")
    }
    // if(igblast_variable_ref_files =~ /^.*IG(K|L)V.*$/){
    //     print("\n\nWARNING:\nLIGHT CHAIN DETECTED IN THE igblast_variable_ref_files parameter.\nBUT CLONAL GROUPING IS GENERALLY RESTRICTED TO HEAVY CHAIN SEQUENCES, AS THE DIVERSITY OF LIGHT CHAINS IS NOT SUFFICIENT TO DISTINGUISH CLONES WITH REASONABLE CERTAINTY")
    // }
    print("\n\nWARNING:\nTHE REPERTOIRE PROCESS CURRENTLY ONLY SUPPORTS ig REFERENCE FILES AND 'mouse' OR 'human' SPECIES")
    print("\n\nWARNING:\nTO MAKE THE REPERTOIRES AND DONUTS, THE SCRIPT CURRENTLY TAKES THE FIRST ANNOTATION OF THE IMGT ANNOTATION IF SEVERAL ARE PRESENTS IN THE v_call, j_call OR c_call COLUMN OF THE productive_seq.tsv FILE")
    print("\n\n")




    //////// end Checks


    //////// Variable modification

    // CONSTRUCTION OF THE igblast REFERENCE FILES PATHS
    // heavy, lambda and kappa were defined earlier in the "Checking of chain coherence section"

    igblast_variable_ref_files = ""
    igblast_constant_ref_files = ""

    def var_files = []
    def const_files = []

    if (heavy) {
        var_files += [
            "imgt_${igblast_organism}_IGHV.fasta",
            "imgt_${igblast_organism}_IGHD.fasta",
            "imgt_${igblast_organism}_IGHJ.fasta"
        ]
        const_files += "imgt_${igblast_organism}_IGHC.fasta"
    }
    if (lambda) {
        var_files += [
            "imgt_${igblast_organism}_IGLV.fasta",
            "imgt_${igblast_organism}_IGLJ.fasta"
        ]
        const_files += "imgt_${igblast_organism}_IGLC.fasta"
    }
    if (kappa) {
        var_files += [
            "imgt_${igblast_organism}_IGKV.fasta",
            "imgt_${igblast_organism}_IGKJ.fasta"
        ]
        const_files += "imgt_${igblast_organism}_IGKC.fasta"
    }

    igblast_variable_ref_files = var_files.join(' ').toString()
    igblast_constant_ref_files = const_files.join(' ').toString()


    //////// end Variable modification


    //////// Channels

    // fs_ch define below because can be a .zip file

    //////// end Channels


    //////// files import

    meta_file = file(meta_path) // in variable because a single file. If "NULL", will create a empty file, present in work folders, but that cannot be correctly linked. Thus, if the file has to be redirected into a channel inside a process, it will not work. Thus, in the first process using meta_file, I hard copy the NULL file if required (see below)
    cute_file = file(cute_path) // in variable because a single file
    phylo_tree_model_file  = file(phylo_tree_model_path)
    template_rmd = file(template_rmd_path)
    alignments_viz_rmd = file(alignments_viz_path)

    //////// end files import


    //////// Main

    if(sample_path =~ /.*\.zip$/){
        Unzip( // warning: it is a process defined above
            Channel.fromPath(sample_path),
            sample_path
        ) 
        fs_ch = Unzip.out.unzip_ch.flatten()
    }else{
        fs_ch = Channel.fromPath("${sample_path}/*.*", checkIfExists: false) // in channel because many files 
    }
    
    nb_input = fs_ch.count()

    workflowParam(
        modules
    )

    file("${out_path}/tsv").mkdirs()
    file("${out_path}/pdf").mkdirs()

    igblast_data_check(
        igblast_organism, 
        igblast_variable_ref_files,
        igblast_constant_ref_files
    )

    Igblast_query(
        fs_ch, 
        igblast_variable_ref_files, 
        igblast_organism, 
        igblast_loci
    )
    //tsv_ch1 = Igblast_query.out.db_pass_ch.map { tuple -> tuple[1] } // Only the igblast_db-pass.tsv files (without _igblast.tsv)
    tsv_ch1 = Igblast_query.out.db_pass_ch
    tsv_ch1.count().subscribe{ n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\n0 ANNOTATION SUCCEEDED BY THE Igblast_query PROCESS\n\nCHECK THAT THE igblast_organism, igblast_loci AND igblast_variable_ref_files ARE CORRECTLY SET IN THE nextflow.config FILE\n\n========\n\n"}}
    tsv_ch2 = tsv_ch1.collectFile(name: "igblast_seq.tsv", skip: 1, keepHeader: true) // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files
    tsv_ch2.subscribe{it -> it.copyTo("${out_path}/tsv")}
    db_unpass = Igblast_query.out.db_unpass_ch.collectFile(name: "failed_igblast_seq.tsv", skip: 1, keepHeader: true) // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files
    db_unpass.subscribe{it -> it.copyTo("${out_path}/tsv")}
    nb1 = tsv_ch2.countLines(keepHeader: true) 
    nb2 =  db_unpass.countLines(keepHeader: true)
    // nb1.view()
    //nb1.view()
    //fs_ch.count().view()
    fs_ch.count().combine(nb1).combine(nb2).subscribe{n,n1,n2 -> if(n != n1 + n2){error "\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nTHE NUMBER OF FILES IN THE igblast_aligned_seq.tsv (${n1}) AND igblast_unaligned_seq_name.tsv (${n2}) IS NOT EQUAL TO THE NUMBER OF SUBMITTED FASTA FILES (${n})\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n"}}
    Igblast_query.out.log_ch.collectFile(name: "Igblast_query_report.log").subscribe{it -> it.copyTo("${out_path}/reports")}


    parseDb_filtering(
        Igblast_query.out.db_pass_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE Igblast_query PROCESS\n\n========\n\n"}
    )
    // parseDb_filtering.out.unproductive_ch.count().subscribe{n -> if ( n == 0 ){print "\n\nWARNING: EMPTY failed_productive_seq.tsv FILE RETURNED FOLLOWING THE parseDb_filtering PROCESS\n\n"}else{it -> it.copyTo("${out_path}/failed_productive_seq.tsv")}} // see https://www.nextflow.io/docs/latest/script.html?highlight=copyto#copy-files
    parseDb_filtering.out.select_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nO PRODUCTIVE SEQUENCE FILES FOLLOWING THE parseDb_filtering PROCESS\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n"}}
    select_ch2 = parseDb_filtering.out.select_ch.collectFile(name: "productive_seq_init.tsv", skip: 1, keepHeader: true) // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files
    parseDb_filtering.out.unselect_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nO UNPRODUCTIVE SEQUENCE FILES FOLLOWING THE parseDb_filtering PROCESS\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n"}} // because an empty file must be present
    unselect_ch2 = parseDb_filtering.out.unselect_ch.collectFile(name: "failed_productive_seq.tsv", skip: 1, keepHeader: true) // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files
    unselect_ch2.subscribe{it -> it.copyTo("${out_path}/tsv")}
    nb1_b = select_ch2.countLines()
    nb2_b = unselect_ch2.countLines()
    nb1_b.subscribe { n -> if ( n == 1 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nO PRODUCTIVE SEQUENCE FOLLOWING THE parseDb_filtering PROCESS\nSEE THE EMPTY productive_seq.tsv FILE AND THE failed_productive_seq.tsv FILE IN THE OUTPUT FOLDER\n\n========\n\n"}}
    // nb1_b.map { "Nombre de lignes dans select_ch2 (productive_seq_init.tsv): $it (en comptant le header)" }.view()
    // nb2_b.map { "Nombre de lignes dans unselect_ch2 (failed_productive_seq.tsv): $it (en comptant le header)" }.view()
    // tsv_ch2.countLines().map { "Nombre de lignes dans tsv_ch2 (igblast_seq.tsv): $it (en comptant le header)" }.view()
    tsv_ch2.countLines().combine(nb1_b).combine(nb2_b).subscribe{n,n1,n2 -> if(n != n1 + n2 - 1){error "\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nTHE NUMBER OF FILES IN THE productive_seq.tsv (${n1}) AND failed_productive_seq.tsv (${n2} - 1) IS NOT EQUAL TO THE NUMBER OF FILES IN igblast_seq.tsv (${n})\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n"}} // n1 + n2 - 1 because header counted in both n1 and n2 while only one header in n
    parseDb_filtering.out.parseDb_filtering_log_ch.collectFile(name: "ParseDb_filtering.log").subscribe{it -> it.copyTo("${out_path}/reports")}


    TrimTranslate(
        parseDb_filtering.out.select_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE parseDb_filtering PROCESS\n\n========\n\n"}
    )
    TrimTranslate.out.trimtranslate_ch.count().subscribe { n -> if ( n == 1 ){error "\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nO PRODUCTIVE SEQUENCE FILES FOLLOWING THE TrimTranslate PROCESS\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n"}}
    trimtranslate_ch2 = TrimTranslate.out.trimtranslate_ch.collectFile(name: "trimtranslate.tsv", skip: 1, keepHeader: true) // productive file with column sequence_alignment_aa added  // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files
    TrimTranslate.out.trimtranslate_log_ch.collectFile(name: "trimtranslate.log").subscribe{it -> it.copyTo("${out_path}/reports")}



    distToNearest(
        trimtranslate_ch2.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE TrimTranslate PROCESS\n\n========\n\n"},
        clone_model,
        clone_normalize
    )
    //distToNearest.out.distToNearest_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE distToNearest PROCESS\n\n========\n\n"}}


    distance_hist(
        distToNearest.out.distToNearest_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE distToNearest PROCESS\n\n========\n\n"},
        cute_file, 
        clone_model,
        clone_normalize,
        clone_distance
    )
    //distance_hist.out.distance_hist_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE distance_hist PROCESS\n\n========\n\n"}}



    histogram_assembly(
        distance_hist.out.histogram_pdf_ch.collect().ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE distance_hist PROCESS\n\n========\n\n"}
    )


 
   seq_name_replacement(
        TrimTranslate.out.trimtranslate_ch,
        meta_file,
        meta_seq_names, 
        meta_name_replacement,
        meta_legend
    )
    seq_name_replacement.out.seq_name_replacement_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE seq_name_replacement PROCESS\n\n========\n\n"}}
    seq_name_replacement_ch2 = seq_name_replacement.out.seq_name_replacement_ch.collectFile(name: "replacement.tsv", skip: 1, keepHeader: true) // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files
    // seq_name_replacement_ch2.subscribe{it -> it.copyTo("${out_path}")}
    // tuple_seq_name_replacement = new Tuple("all", seq_name_replacement_ch2) # warning: this is not a channel but a variable now
    seq_name_replacement.out.seq_name_replacement_log_ch.collectFile(name: "seq_name_replacement.log").subscribe{it -> it.copyTo("${out_path}/reports")}


    data_assembly(
        seq_name_replacement_ch2.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE seq_name_replacement PROCESS\n\n========\n\n"}, 
        distToNearest.out.distToNearest_ch
    )
    data_assembly.out.productive_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE data_assembly PROCESS\n\n========\n\n"}}

    metadata_check(
        data_assembly.out.productive_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE data_assembly PROCESS\n\n========\n\n"},
        meta_file, 
        meta_seq_names, 
        meta_name_replacement,
        meta_legend
    )

    repertoire(
        data_assembly.out.productive_ch,
        igblast_data_check.out.igblast_data_check_ch,
        cute_file
    )

    clone_assignment(
        data_assembly.out.productive_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE data_assembly PROCESS\n\n========\n\n"}, 
        clone_model,
        clone_normalize,
        clone_distance,
        clone_strategy,
        meta_file,
        meta_legend
    )
    // clone_assignment.out.clone_ch.view()
    // nb_failed_clone = clone_assignment.out.failed_clone_ch.map { file -> file.countLines() == 0 ? 0: file.countLines() - 1} // either failed_clone_assigned_seq.tsv is empty (no lines) or has a header to remove to count the lines
    nb_failed_clone_1 = clone_assignment.out.failed_clone_ch.countLines() - 1 // Minus 1 because 1st line = column names // either failed_clone_assigned_seq.tsv is empty (no lines) or has a header to remove to count the lines

    split_by_clones(
        clone_assignment.out.clone_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE clone_assignment PROCESS\n\n========\n\n"}
    )
    tempo_test = split_by_clones.out.clone_split_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE split_by_clones PROCESS\n\n========\n\n"}


    Closest_germline(
        split_by_clones.out.clone_split_ch.flatten(), // flatten split the list into several objects (required for parallelization)
        igblast_organism, 
        igblast_variable_ref_files, 
        clone_germline_kind,
        meta_file,
        meta_legend
    )
    failed_clonal_germline_file = Closest_germline.out.failed_clonal_germline_ch.collectFile(name: "failed_clonal_germline.tsv", skip: 1, keepHeader: true)
    failed_clonal_germline_file.subscribe{it -> it.copyTo("${out_path}/tsv")} // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files
    nb_failed_clone_2 = failed_clonal_germline_file.countLines() - 1
    nb_failed_clone = nb_failed_clone_2.combine(nb_failed_clone_1).map { new_count, old_count -> new_count + old_count }.view() // to add the two kind of failure
    Closest_germline.out.closest_log_ch.collectFile(name: "Closest_germline.log").subscribe{it -> it.copyTo("${out_path}/reports")}


    Igblast_germline_coords(
        Closest_germline.out.closest_ch, // no ifEmpty{error} because can be empty for some of the parallelized processes
        igblast_organism, 
        igblast_loci
    )
    //tsv_ch1 = Igblast_query.out.db_pass_ch.map { tuple -> tuple[1] } // Only the igblast_db-pass.tsv files (without _igblast.tsv)
    Igblast_germline_coords.out.germline_coords_ch.count().subscribe{ n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\n0 ANNOTATION SUCCEEDED BY THE Igblast_germline_coords PROCESS\n\nCHECK THAT THE igblast_organism, igblast_loci AND igblast_variable_ref_files ARE CORRECTLY SET IN THE nextflow.config FILE\n\n========\n\n"}}
    Igblast_germline_coords.out.log_ch.collectFile(name: "Igblast_germline_coords_report.log").subscribe{it -> it.copyTo("${out_path}/reports")}


    AddGermlineSequences(
        Igblast_germline_coords.out.germline_coords_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE Igblast_germline_coords PROCESS\n\n========\n\n"},
        igblast_organism, 
        igblast_variable_ref_files
    )
    AddGermlineSequences.out.add_germ_log_ch.collectFile(name: "AddGermlineSequences.log").subscribe{it -> it.copyTo("${out_path}/reports")}



    TranslateGermline(
        AddGermlineSequences.out.add_germ_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE AddGermlineSequences PROCESS\n\n========\n\n"}
    )
    TranslateGermline.out.translate_germ_log_ch.collectFile(name: "TranslateGermline.log").subscribe{it -> it.copyTo("${out_path}/reports")}




    Mutation_load_germ_genes(
        TranslateGermline.out.translate_germ_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE TranslateGermline PROCESS\n\n========\n\n"},
        meta_file, 
        meta_legend,
        clone_mut_obs_seq,
        clone_mut_germ_seq,
        clone_mut_regionDefinition
    )
    clone_assigned_seq_ch = Mutation_load_germ_genes.out.mutation_load_ch.collectFile(name: "clone_assigned_seq.tsv", skip: 1, keepHeader: true) // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files
    nb_clone_assigned = clone_assigned_seq_ch.countLines() - 1 // Minus 1 because 1st line = column names
    clone_assigned_seq_ch.subscribe{it -> it.copyTo("${out_path}/tsv")}
    clone_assigned_seq_filtered_ch = Mutation_load_germ_genes.out.mutation_load_ch.filter{ file -> file.countLines() > align_clone_nb.toInteger() } // Only keep clonal groups that have a number of sequences superior to align_clone_nb (variable defined in nextflow.config)
    Mutation_load_germ_genes.out.mutation_load_log_ch.collectFile(name: "Mutation_load_germ_genes.log").subscribe{it -> it.copyTo("${out_path}/reports")} // 




    Clone_id_count(
        clone_assigned_seq_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE Mutation_load_germ_genes PROCESS\n\n========\n\n"}
    )





    /*

    get_germ_tree(
        Mutation_load_germ_genes.out.mutation_load_ch,
        meta_file, // first() because get_germ_tree process is a parallele one and because meta_file is single
        cute_file, 
        align_clone_nb,
        germ_tree_duplicate_seq,
        igphylm_exe_path
    )
    
    get_germ_tree.out.rdata_germ_tree_ch.count().subscribe { n -> if ( n == 0 ){
        print("\n\nWARNING: EMPTY OUTPUT FOLLOWING THE get_germ_tree PROCESS -> NO TREE RETURNED\n\n")
    }}
    rdata_germ_tree_ch2 = get_germ_tree.out.rdata_germ_tree_ch.collect()
    //rdata_germ_tree_ch2.view()

    get_germ_tree.out.no_germ_tree_ch.count().subscribe { n -> if ( n == 0 ){
        print("\n\nWARNING: ALL SEQUENCES IN TREES FOLLOWING THE get_germ_tree PROCESS -> EMPTY germ_tree_dismissed_seq.tsv FILE RETURNED\n\n")
    }}
    no_germ_tree_ch2 = get_germ_tree.out.no_germ_tree_ch.collectFile(name: "germ_tree_dismissed_seq.tsv", skip: 1, keepHeader: true)
    no_germ_tree_ch2.subscribe{it -> it.copyTo("${out_path}/tsv")}

    get_germ_tree.out.germ_tree_ch.count().subscribe { n -> if ( n == 0 ){
        print("\n\nWARNING: NO SEQUENCES IN TREES FOLLOWING THE get_germ_tree PROCESS -> EMPTY seq_for_germ_tree.tsv FILE RETURNED\n\n")
    }}
    germ_tree_ch2 = get_germ_tree.out.germ_tree_ch.collectFile(name: "seq_for_germ_tree.tsv", skip: 1, keepHeader: true)
    // germ_tree_ch2.subscribe{it -> it.copyTo("${out_path}/tsv")}
    germ_tree_ch3 = get_germ_tree.out.germ_tree_ch.flatten().filter{ file -> file.countLines() > align_clone_nb.toInteger() }

    get_germ_tree.out.no_cloneID_ch.count().subscribe { n -> if ( n == 0 ){
        print("\n\nWARNING: ALL SEQUENCES IN CLONAL GROUP FOLLOWING THE get_germ_tree PROCESS -> EMPTY germ_tree_dismissed_clone_id.tsv FILE RETURNED\n\n")
    }}
    no_cloneID_ch2 = get_germ_tree.out.no_cloneID_ch.collectFile(name: "germ_tree_dismissed_clone_id.tsv")
    no_cloneID_ch2.subscribe{it -> it.copyTo("${out_path}/tsv")}

    get_germ_tree.out.cloneID_ch.count().subscribe { n -> if ( n == 0 ){
        print("\n\nWARNING: NO CLONAL GROUP FOLLOWING THE get_germ_tree PROCESS -> EMPTY germ_tree_clone_id.tsv and germ_tree.pdf FILES RETURNED\n\n")
    }}
    cloneID_ch2 = get_germ_tree.out.cloneID_ch.collectFile(name: "germ_tree_clone_id.tsv")
    cloneID_ch2.subscribe{it -> it.copyTo("${out_path}/tsv")}

    get_germ_tree.out.get_germ_tree_log_ch.collectFile(name: "get_germ_tree.log").subscribe{it -> it.copyTo("${out_path}/reports")} // 



    germ_tree_vizu(
        rdata_germ_tree_ch2,
        germ_tree_kind,
        align_clone_nb,
        germ_tree_duplicate_seq,
        germ_tree_leaf_color,
        germ_tree_leaf_shape,
        germ_tree_leaf_size,
        germ_tree_label_size,
        germ_tree_label_hjust,
        germ_tree_label_rigth,
        germ_tree_label_outside,
        germ_tree_right_margin,
        germ_tree_legend,
        clone_assigned_seq_ch, // may be add .first()
        meta_file, 
        meta_legend,
        cute_file
    )
    germ_tree_vizu.out.germ_tree_vizu_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE germ_tree_vizu PROCESS\n\n========\n\n"}
    //germ_tree_vizu.out.germ_tree_vizu_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE germ_tree_vizu PROCESS\n\n========\n\n"}}
    germ_tree_vizu.out.germ_tree_dup_seq_not_displayed_ch.ifEmpty{
        print("\n\nWARNING: -> NO germ_tree_dup_seq_not_displayed.tsv FILE RETURNED\n\n")
    }
    germ_tree_dup_seq_not_displayed_ch2 = germ_tree_vizu.out.germ_tree_dup_seq_not_displayed_ch.flatten().collectFile(name: "germ_tree_dup_seq_not_displayed.tsv", skip: 1, keepHeader: true) // flatten split the list into several objects which is required by collectFile()
    germ_tree_dup_seq_not_displayed_ch2.subscribe{it -> it.copyTo("${out_path}/tsv")}

    */



    


    tempo1_ch = Channel.of("all", "annotated", "tree") // 1 channel with 3 values (not list)
    tempo2_ch = data_assembly.out.productive_ch.mix(data_assembly.out.productive_ch.mix(clone_assigned_seq_ch)) // 1 channel with 3 paths (do not use flatten() -> not list)
    tempo3_ch = tempo1_ch.merge(tempo2_ch) // 3 lists
    tempo4_ch = Channel.of("vj_allele", "c_allele", "vj_gene", "c_gene")
    tempo5_ch = tempo3_ch.combine(tempo4_ch) // 12 tuples

    donut(
        tempo5_ch,
        donut_palette,
        donut_hole_size,
        donut_hole_text,
        donut_hole_text_size,
        donut_border_color,
        donut_border_size,
        donut_annotation_distance,
        donut_annotation_size,
        donut_annotation_force,
        donut_annotation_force_pull,
        donut_legend_width,
        donut_legend_text_size,
        donut_legend_box_size,
        donut_legend_box_space,
        donut_legend_limit,
        cute_file,
        igblast_data_check.out.igblast_data_check_ch
    )
    donut.out.donut_pdf_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE donut PROCESS\n\n========\n\n"}
    donut.out.donut_pdf_ch.count().subscribe { n -> if ( n == 0 ){
        print("\n\nWARNING: EMPTY OUTPUT FOLLOWING THE donut PROCESS -> NO DONUT RETURNED\n\n")}}
    donut_pdf_ch2 = donut.out.donut_pdf_ch.collect()

    donut.out.donut_tsv_ch.count().subscribe { n -> if ( n == 0 ){
        print("\n\nWARNING: -> NO donut_stats.tsv FILE RETURNED FOLLOWING THE donut PROCESS\n\n")
    }}
    donut_tsv_ch2 = donut.out.donut_tsv_ch.collectFile(name: "donut_stats.tsv", skip: 1, keepHeader: true) // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files
    donut_tsv_ch2.subscribe{it -> it.copyTo("${out_path}/tsv")}



    donut_assembly(
        donut_pdf_ch2
    )
    //donut_assembly.out.donut_assembly_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE donut_assembly PROCESS\n\n========\n\n"}
    donut_assembly.out.donut_assembly_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE donut_assembly PROCESS\n\n========\n\n"}}

/*
    Reformat(
        aa_tsv_ch2,
        igblast_data_check.out.igblast_data_check_ch.collect() // To make sure the igblast_ref files are in the right format before executing following processes
    )
    fasta = Reformat.out.aa_fasta_ch


    DefineGroups(
        fasta,
        select_ch2
    )
    fastagroups = DefineGroups.out.groups_ch.flatten()
    align_input = fastagroups.combine(heavy_chain)
*/


    clone_labeled_ch = clone_assigned_seq_filtered_ch.map { file -> tuple(file, 'CLONE') }
    productive_labeled_ch = data_assembly.out.productive_ch.map { file -> tuple(file, 'ALL') }
    all_files_ch = clone_labeled_ch.mix(productive_labeled_ch)
    Tsv2fasta(
        all_files_ch,
        align_seq, 
        clone_germline_kind, 
        align_clone_nb, 
        cute_path
    )
    tsv2fasta_log_ch2 = Tsv2fasta.out.tsv2fasta_log_ch.collectFile(name: "Tsv2fasta.log")
    tsv2fasta_log_ch2.subscribe{it -> it.copyTo("${out_path}/reports")}
    fasta_align_warn = Tsv2fasta.out.warning_ch

    // Print warnings on the terminal:
    fasta_align_warn.filter { file(it). exists() }
                .map {file -> 
                    file.text  // contenu du fichier
                }
                .view()

    if(align_soft == "abalign" && (align_seq == "query" || align_seq == "igblast_full" || align_seq == "trimmed" || align_seq == "fwr1" || align_seq == "fwr2" || align_seq == "fwr3" || align_seq == "fwr4" || align_seq == "cdr1" || align_seq == "cdr2" || align_seq == "cdr3" || align_seq == "junction" || align_seq == "d_sequence_alignment" || align_seq == "j_sequence_alignment" || align_seq == "c_sequence_alignment" || align_seq == "d_germline_alignment" || align_seq == "j_germline_alignment" || align_seq == "c_germline_alignment")){
        align_soft = "mafft"
        print("\n\nWARNING: align_soft PARAMETER RESET TO \"mafft\" SINCE align_soft PARAMETER WAS SET TO \"abalign\" BUT THAT IT REQUIRES AT LEAST A V DOMAIN IN THE SEQUENCES,\nWHILE align_seq PARAMETER IS SET TO \"${align_seq}\"\n\n")
    }
    if(align_soft == "mafft"){
        Mafft_align(
            Tsv2fasta.out.fasta_align_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE Tsv2fasta PROCESS\n\n========\n\n"},
            align_mafft_all_options,
            align_mafft_clonal_options
        )
        align_aa_ch = Mafft_align.out.aligned_all_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE Mafft_align PROCESS\n\n========\n\n"}.map{ x, y, z -> [y, z] }
        align_nuc_ch = Mafft_align.out.aligned_all_ch.map{ x, y, z -> [x, z] }
        aligned_all_ch2 = Mafft_align.out.aligned_all_ch.map{ x, y, z -> [x, y] }
        Mafft_align.out.mafft_align_log_ch.collectFile(name: "Mafft_align.log").subscribe{it -> it.copyTo("${out_path}/reports")}
    }else if(align_soft == "abalign"){
        Abalign_align_aa(
            Tsv2fasta.out.fasta_align_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE Tsv2fasta PROCESS\n\n========\n\n"},
            igblast_organism,
            igblast_heavy_chain,
            align_abalign_options
        )
        Abalign_align_aa.out.abalign_align_aa_log_ch.collectFile(name: "Abalign_align_aa.log").subscribe{it -> it.copyTo("${out_path}/reports")}

        Abalign_rename(
            Abalign_align_aa.out.aligned_aa_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE Abalign_align_aa PROCESS\n\n========\n\n"}
        )
        Abalign_rename.out.failed_abalign_align_ch.collectFile(name: "failed_abalign_align.tsv", skip: 1, keepHeader: true).subscribe{it -> it.copyTo("${out_path}/tsv")}
        Abalign_rename.out.abalign_rename_log_ch.collectFile(name: "Abalign_rename.log").subscribe{it -> it.copyTo("${out_path}/reports")}

        // aligned_aa_only_ch = Abalign_align_aa.out.aligned_aa_ch.map { x, y, z -> y }
        Abalign_align_nuc(
            Abalign_rename.out.renamed_aligned_aa_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE Abalign_align_aa PROCESS\n\n========\n\n"}
        )
        align_nuc_ch = Abalign_align_nuc.out.aligned_all_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE Abalign_align_nuc PROCESS\n\n========\n\n"}.map{ x, y, z -> [x, z] }
        align_aa_ch = Abalign_align_nuc.out.aligned_all_ch.map{ x, y, z -> [y, z] }
        aligned_all_ch2 = Abalign_align_nuc.out.aligned_all_ch.map{ x, y, z -> [x, y] }
        Abalign_align_nuc.out.abalign_align_nuc_log_ch.collectFile(name: "Abalign_align_nuc.log").subscribe{it -> it.copyTo("${out_path}/reports")}
    }else{
        error "\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nINVALID align_soft PARAMETER IN nextflow.config FILE:\n${align_soft}\n\n========\n\n"
    }

    //Abalign_align_nuc.out.aligned_all_ch.map{ x, y, z -> [y, z] }.view()
    //Mutation_load_germ_genes.out.mutation_load_ch.view()
    branches_nuc = align_nuc_ch.branch {
        CLONE: it[1] == 'CLONE'
        ALL  : it[1] == 'ALL'
    }
    clone_for_gff_nuc_ch = branches_nuc.CLONE.combine(clone_assigned_seq_filtered_ch.first())
    all_for_gff_nuc_ch  = branches_nuc.ALL.combine(data_assembly.out.productive_ch.first())
    for_gff_nuc_ch = clone_for_gff_nuc_ch.mix(all_for_gff_nuc_ch)
    branches_aa = align_aa_ch.branch {
        CLONE: it[1] == 'CLONE'
        ALL  : it[1] == 'ALL'
    }
    clone_for_gff_aa_ch = branches_aa.CLONE.combine(clone_assigned_seq_filtered_ch.first())
    all_for_gff_aa_ch  = branches_aa.ALL.combine(data_assembly.out.productive_ch.first())
    for_gff_aa_ch = clone_for_gff_aa_ch.mix(all_for_gff_aa_ch)


    GffNuc( // module gff.nf
        for_gff_nuc_ch,
        align_seq, 
        align_clone_nb, 
        cute_path
    )
    GffNuc.out.gff_ch.flatten().ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE GffNuc PROCESS\n\n========\n\n"}.subscribe{it -> it.copyTo("${out_path}/alignments/nuc")} // flatten because several files. Otherwise, copyTo does not like it
    GffNuc.out.gff_log_ch.collectFile(name: "GffNuc.log").subscribe{it -> it.copyTo("${out_path}/reports")}


    //for_gff_aa_ch.view()
    GffAa( // module gff.nf
        for_gff_aa_ch,
        align_seq, 
        align_clone_nb, 
        cute_path
    )
    GffAa.out.gff_ch.flatten().ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE GffAa PROCESS\n\n========\n\n"}.subscribe{it -> it.copyTo("${out_path}/alignments/aa")}
    GffAa.out.gff_log_ch.collectFile(name: "GffAa.log").subscribe{it -> it.copyTo("${out_path}/reports")}
//GffAa.out.gff_ch.flatten().view()


    PrintAlignmentNuc( // module print_alignment.nf
        align_nuc_ch
    )
    PrintAlignmentNuc.out.alignment_html.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE PrintAlignmentNuc PROCESS\n\n========\n\n"}.subscribe{it -> it.copyTo("${out_path}/alignments/nuc")}

    PrintAlignmentAa( // module print_alignment.nf
        align_aa_ch
    )
    PrintAlignmentAa.out.alignment_html.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE PrintAlignmentAa PROCESS\n\n========\n\n"}.subscribe{it -> it.copyTo("${out_path}/alignments/aa")}

    Tree(
        aligned_all_ch2,
        phylo_tree_model_file
    )
    tree_aa = Tree.out.tree_aa_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE Tree PROCESS FOR NUC\n\n========\n\n"}
    tree_nuc = Tree.out.tree_nuc_ch
    Tree.out.tree_log_ch.collectFile(name: "tree.log").subscribe{it -> it.copyTo("${out_path}/reports")}

    Meta2ITOL (
        meta_file,
        meta_seq_names
    )

    // The ITOL process can only be executed if user has paid the subsription for automated visualization

    if(phylo_tree_itol_subscription == "TRUE"){

        tree = tree_aa.concat(tree_nuc)

        ITOL(
            tree,
            Meta2ITOL .out.itol_out_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE Meta2ITOL  PROCESS\n\n========\n\n"},
            phylo_tree_itolkey
        )

        itol_aa_ch = ITOL.out.itol_aa_ch
        itol_nuc_ch = ITOL.out.itol_nuc_ch

        itol_aa_ch.subscribe{it -> it.copyTo("${out_path}/phylo/aa")}
        itol_nuc_ch.subscribe{it -> it.copyTo("${out_path}/phylo/nuc")}

    }



    repertoire_constant_ch =repertoire.out.repertoire_png_ch
        .flatten()
        .filter { file -> 
            file.name =~ /^.*IG.C_.*gene_non-zero\.png$/
        }
    repertoire_vj_ch = repertoire.out.repertoire_png_ch
        .flatten()
        .filter { file ->
            file.name =~ /.*\/?(rep_gene_IG.V_.*non-zero\.png)$/
        }



    print_report(
        template_rmd,
        alignments_viz_rmd, 
        nb_input,
        nb1,
        nb2,
        nb1_b - 1, // Minus 1 because the 1st line = the column names
        nb2_b - 1, // Minus 1 because the 1st line = the column names
        nb_clone_assigned,
        nb_failed_clone,
        donut.out.donuts_png.collect(),
        repertoire.out.repertoire_png_ch.collect(),
        repertoire_constant_ch,
        repertoire_vj_ch,
        align_soft,
        phylo_tree_itol_subscription
    )



    backup(
        config_file, 
        log_file
    )

}

    //////// end Main


//////// end Processes