nextflow.enable.dsl=2
/*
#########################################################################
##                                                                     ##
##     main.nf                                                         ##
##                                                                     ##
##     repertoire_profiler                                             ##
##     Gael A. Millot                                                  ##
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




// next process is for repertoire
process igblast_data_check { // cannot be igblast_data_check outside of process because the files are present in the docker
    label 'immcantation'
    //publishDir path: "${out_path}", mode: 'copy', pattern: "{*.tsv}", overwrite: false
    cache 'true'

    input:
    val igblast_database_path
    val igblast_files

    output:
    path "*.tsv", emit: igblast_data_check_ch

    script:
    """
    #!/bin/bash -ue
    REPO_PATH="/usr/local/share/${igblast_database_path}" # path where the imgt_human_IGHV.fasta, imgt_human_IGHD.fasta and imgt_human_IGHJ.fasta files are in the docker container
    VDJ_FILES=\$(awk -v var1="${igblast_files}" -v var2="\${REPO_PATH}" 'BEGIN{ORS=" " ; split(var1, array1, " ") ; for (key in array1) {print var2"/"array1[key]}}') # assemble files with their path
    for i1 in \$VDJ_FILES ; do
        if [[ ! -e \${i1} ]] ; then
            echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nFILE DOES NOT EXISTS:\\n\${i1}\\n\\nINDICATED PATH:\\n\${REPO_PATH}\\n\\nCONTAINS:\\n"
            ls -la -w 1 \${REPO_PATH}
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



process igblast {
    label 'immcantation'
    //publishDir path: "${out_path}", mode: 'copy', overwrite: false
    cache 'true'

    input:
    path fs_ch // parallelization expected (for each fasta file)
    val igblast_database_path
    val igblast_files
    val igblast_organism
    val igblast_loci
    val igblast_aa

    output:
    path "*_igblast_db-pass.tsv", emit: tsv_ch1, optional: true
    path "igblast_unaligned_seq.tsv", emit: unaligned_seq_ch
    path "igblast_aligned_seq.tsv", emit: aligned_seq_ch
    path "*.log", emit: log_ch, optional: true

    script:
    """
    #!/bin/bash -ue

    # variables
    REPO_PATH="/usr/local/share/${igblast_database_path}" # path where the imgt_human_IGHV.fasta, imgt_human_IGHD.fasta and imgt_human_IGHJ.fasta files are in the docker container
    VDJ_FILES=\$(awk -v var1="${igblast_files}" -v var2="\${REPO_PATH}" 'BEGIN{ORS=" " ; split(var1, array1, " ") ; for (key in array1) {print var2"/"array1[key]}}')
    FILENAME=\$(basename -- "${fs_ch}") # recover a file name without path
    FILE=\${FILENAME%.*} # file name without extension
    FILE_EXTENSION="\${FILENAME##*.}" #  ## means "delete the longest regex starting at the beginning of the tested string". If nothing, delete nothing. Thus ##*. means delete the longest string finishing by a dot. Use # instead of ## for "delete the shortest regex starting at the beginning of the tested string"
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a igblast_report.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a igblast_report.log
    # end variables

    # checks
    if [[ ! "\${FILE_EXTENSION}" =~ fasta|fa|fna|txt|seq ]] ; then
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nINVALID FILE EXTENSION IN THE sample_path PARAMETER OF THE repertoire_profiler.config FILE:\\n${fs_ch}\\n\${FILENAME}\\nMUST BE fasta|fs_ch|txt|seq\\n\\n========\\n\\n"
        exit 1
    fi
    sed 's/\\r\$//' ${fs_ch} > tempo_file.fasta # remove carriage returns
    awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' tempo_file.fasta > \${FILE}.fa # remove \\n in the middle of the sequence # \${FILENAME}.fa is a trick to do not use ${fs_ch} and modify the initial file due to the link in the nextflow work folder
    TEMPO=\$(wc -l \${FILE}.fa | cut -f1 -d' ')
    if read -n 1 char <"\${FILE}.fa"; [[ \$char != ">" || \$TEMPO != 2 ]]; then
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nINVALID FASTA FILE IN THE sample_path OF THE repertoire_profiler.config FILE:\\n${fs_ch}\\nMUST BE A FASTA FILE (\'>\' AS FIRST CHARATER) MADE OF A SINGLE SEQUENCE\\n\\n========\\n\\n"
        exit 1
    fi
    # end checks

    # Alignment <-> annotate sequence using VDJ info
    # See https://changeo.readthedocs.io/en/stable/tools/AssignGenes.html for the details
    if [[ "${igblast_aa}" == "false" ]] ; then
        AssignGenes.py igblast -s \${FILE}.fa -b /usr/local/share/igblast --organism ${igblast_organism} --loci ${igblast_loci} --format blast |& tee -a igblast_report.log
    else
        # WARNING: does not work if the fasta file contains \\r (CR carriage return, instead or in addition of \\n, LF line feed) but ok here because removed above
        awk -v var1=\${FILENAME} '{lineKind=(NR-1)%2;}lineKind==0{record=\$0 ; next}lineKind==1{if(\$0~/^[NATGC]*\$/){print "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nFASTA FILE\\n"var1"\\nMUST BE AN AMINO ACID SEQUENCE IF THE igblast_aa PARAMETER IS SET TO true\\nHERE IT SEEMS ONLY MADE OF NUCLEOTIDES:\\n"\$0"\\n\\n========\\n\\n" ; exit 1}}' \${FILE}.fa
        AssignGenes.py igblast-aa -s \${FILE}.fa -b /usr/local/share/igblast --organism ${igblast_organism} --loci ${igblast_loci} |& tee -a igblast_report.log
    fi
    # convert to tsv
    # Also convert data from the web interface IMGT/HighV-QUEST
    if [[ "${igblast_aa}" == "false" ]] ; then
        MakeDb.py igblast -i ./\${FILE}_igblast.fmt7 -s ./\${FILE}.fa -r \${VDJ_FILES} --extended |& tee -a igblast_report.log
    else
        MakeDb.py igblast-aa -i ./\${FILE}_igblast.fmt7 -s ./\${FILE}.fa -r \${VDJ_FILES} --extended |& tee -a igblast_report.log
    fi
    # printing if no tsv file made
    if [[ ! -f ./\${FILE}_igblast_db-pass.tsv ]] ; then
        echo "\${FILENAME}" | cat > igblast_unaligned_seq.tsv
        echo -n "" | cat > igblast_aligned_seq.tsv
    else
        echo "\${FILENAME}" | cat > igblast_aligned_seq.tsv
        echo -n "" | cat > igblast_unaligned_seq.tsv
    fi
    """
    // write ${} between "" to make a single argument when the variable is made of several values separated by a space. Otherwise, several arguments will be considered
}


process parseDb_filtering {
    label 'immcantation'
    publishDir path: "${out_path}", mode: 'copy', pattern: "unproductive_seq.tsv", overwrite: false
    publishDir path: "${out_path}", mode: 'copy', pattern: "productive_seq.tsv", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{ParseDb_filtering.log}", overwrite: false
    cache 'true'

    input:
    path tsv_ch2 // no more parallelization
    val igblast_aa

    output:
    path "productive_seq.tsv", emit: select_ch
    path "unproductive_seq.tsv"
    path "ParseDb_filtering.log"

    script:
    """
    #!/bin/bash -ue
    if [[ "${igblast_aa}" == "true" ]]; then # if igblast_aa is true, then the productive column is empty because aa sequences means productive ones
        FILENAME=\$(basename -- "${tsv_ch2}") # recover a file name without path
        FILE=\${FILENAME%.*} # file name without extension
        ParseDb.py select -d ${tsv_ch2} -f v_call j_call -u ".+" --regex --logic any |& tee -a ParseDb_filtering.log #means look inside the -f v_call j_call fields of the input and return any lines that are non empty for at least one field(--logic any) # should be identical to cp ${tsv_ch2} "\${FILE}_parse-select.tsv" |& tee -a ParseDb_filtering.log
        echo -n "" > unproductive_seq.tsv |& tee -a ParseDb_filtering.log
    elif [[ -s ${tsv_ch2} ]]; then # -s means "exists and non empty". Thus, return FALSE is the file does not exists or is empty
        ParseDb.py select -d ${tsv_ch2} -f productive -u T |& tee -a ParseDb_filtering.log
        ParseDb.py split -d ${tsv_ch2} -f productive |& tee -a ParseDb_filtering.log
        if [ -f *_parse-select.tsv ]; then
            cp *_parse-select.tsv productive_seq.tsv |& tee -a ParseDb_filtering.log
        else
            echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nNO *_parse-select.tsv FILE GENERATED BY THE igblast PROCESS\\nCHECK THE ParseDb_filtering.log FILE IN THE report FOLDER\\n\\n========\\n\\n"
            exit 1
        fi
        if [ -s *_productive-F.tsv ]; then # see above for -s
            cp *_productive-F.tsv unproductive_seq.tsv |& tee -a ParseDb_filtering.log
        else
            echo -e "\n\nWARNING: EMPTY unproductive_seq.tsv FILE RETURNED FOLLOWING THE parseDb_filtering PROCESS\n\n" |& tee -a ParseDb_filtering.log
            echo -n "" | cat > unproductive_seq.tsv
        fi
    else
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nEMPTY FILE GENERATED BY THE igblast PROCESS\\nCHECK THE ParseDb_filtering.log FILE IN THE report FOLDER\\n\\n========\\n\\n" # I said empty because the existence of the file has been checked after the igblast process
        exit 1
    fi
    """
}



process translation {
    label 'seqkit'
    publishDir path: "${out_path}", mode: 'copy', pattern: "aa/*.fasta", overwrite: false
    publishDir path: "${out_path}", mode: 'copy', pattern: "aligned_seq/*.fasta", overwrite: false
    publishDir path: "${out_path}", mode: 'copy', pattern: "aa.tsv", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{translation.log}", overwrite: false
    cache 'true'

    input:
    path select_ch // no more parallelization
    val igblast_aa

    output:
    path "translation.tsv" , emit: translation_ch
    path "aa.tsv" , optional: true
    path "aligned_seq/*.*"
    path "aa/*.*" , optional: true
    path "translation.log", emit: translation_log_ch

    script:
    """
    #!/bin/bash -ue
    translation.sh ${select_ch} ${igblast_aa}
    """
}



process distToNearest {
    label 'immcantation'
    //publishDir path: "${out_path}", mode: 'copy', pattern: "{nearest_distance.tsv}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{distToNearest.log}", overwrite: false
    cache 'true'

    input:
    path translation_ch // no parallelization
    val igblast_aa
    val clone_model
    val clone_normalize

    output:
    path "nearest_distance.tsv", emit: distToNearest_ch
    path "distToNearest.log"

    script:
    """
    #!/bin/bash -ue
    Rscript -e '
         # WEIRD stuf: if db alone is returned, and if distToNearest_ch is used for the clone_assignment process and followings, everything is fine. But if db3 is returned with db3 <- data.frame(db, dist_nearest = db2\$dist_nearest) or db3 <- data.frame(db, caca = db2\$dist_nearest) or data.frame(db, caca = db\$sequence_id) or db3 <- data.frame(db, caca = as.numeric(db2\$dist_nearest)) or db3 <- data.frame(db[1:3], caca = db\$sequence_id, db[4:length(db)]), the get_tree process cannot make trees, while the productive.tsv seem identical at the end, between the use of db or db3, except that the clone_id order is not the same
        db <- read.table("${translation_ch}", header = TRUE, sep = "\\t")
        if("${clone_model}" != "aa" & "${igblast_aa}" == "true"){
          tempo.cat <- paste0("\\n\\n========\\n\\nERROR IN THE NEXTFLOW EXECUTION OF THE distToNearest PROCESS\\nclone_model PARAMETER SHOULD BE \\"aa\\" IF AA FASTA FILES ARE USED (igblast_aa PARAMETER SET TO \\"true\\"). HERE:\\nclone_model: ${clone_model}\\n\\n========\\n\\n")
          stop(tempo.cat)
        }
        db2 <- shazam::distToNearest(db, sequenceColumn = "junction", locusColumn = "locus", model = "${clone_model}", normalize = "${clone_normalize}", nproc = 1)
        write.table(db2, file = paste0("./nearest_distance.tsv"), row.names = FALSE, col.names = TRUE, sep = "\\t")
    ' |& tee -a distToNearest.log
    """
}


process distance_hist {
    label 'immcantation'
    publishDir path: "${out_path}/png", mode: 'copy', pattern: "{*.png}", overwrite: false
    publishDir path: "${out_path}/svg", mode: 'copy', pattern: "{*.svg}", overwrite: false
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
    label 'r_ext'
    publishDir path: "${out_path}", mode: 'copy', pattern: "{seq_distance.pdf}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{histogram_assembly.log}", overwrite: false
    cache 'true'

    input:
    path histogram_pdf_ch

    output:
    path "seq_distance.pdf"

    script:
    """
    #!/bin/bash -ue
    Rscript -e '
    # assignation to prevent a returned element
        tempo <- qpdf::pdf_combine(input = list.files(path = ".", pattern = "^seq.*.pdf\$"), output = "./seq_distance.pdf")
    ' |& tee -a histogram_assembly.log
    """
}


process clone_assignment {
    label 'immcantation'
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{*.log}", overwrite: false
    publishDir path: "${out_path}", mode: 'copy', pattern: "non_clone_assigned_sequence.tsv", overwrite: false
    cache 'true'

    input:
    path translation_ch // no parallelization
    val clone_model
    val clone_normalize
    val clone_distance

    output:
    path "*_clone-pass.tsv", emit: clone_ch
    path "non_clone_assigned_sequence.tsv", emit: failed_clone_ch
    path "*.log"

    script:
    """
    #!/bin/bash -ue
    if [[ -s ${translation_ch} ]]; then # see above for -s
        DefineClones.py -d ${translation_ch} --act set --model ${clone_model} --norm ${clone_normalize} --dist ${clone_distance} --fail |& tee -a clone_assignment.log
        if [ -s *_clone-fail.tsv ]; then # see above for -s
            cp *_clone-fail.tsv non_clone_assigned_sequence.tsv |& tee -a clone_assignment.log
        else
            echo -e "\n\nNOTE: EMPTY non_clone_assigned_sequence.tsv FILE RETURNED FOLLOWING THE clone_assignment PROCESS\n\n" |& tee -a clone_assignment.log
            echo -n "" | cat > non_clone_assigned_sequence.tsv
        fi
    else
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nEMPTY FILE GENERATED BY THE parseDb_filtering PROCESS\\nCHECK THE clone_assignment.log AND *_productive-F.tsv FILES IN THE report FOLDER\\n\\n========\\n\\n"
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
    path "*_translation_clone-pass.tsv", emit: clone_split_ch // multiple files -> parall expected

    script:
    """
    #!/bin/bash -ue
    FILENAME=\$(basename -- ${clone_ch}) # recover a file name without path
    cp -Lr ${clone_ch} "./TEMPO.tsv" # to have the hard file, not the symlink, because modifications will be performed inside
    chmod 777 TEMPO.tsv
    rm \$FILENAME # remove the initial file to avoid to send it into the channel
    cp -rp TEMPO.tsv "\$FILENAME" # -p for preserve permissions
    rm TEMPO.tsv
    awk -v var1=\$FILENAME -F "\\t" '{
        if(NR == 1){
            header=\$0 ; next
        }else{
            if(system( "[ -f " \$NF"_"var1 " ] " ) > 0){ # test if a file exists
                print header > \$NF"_"var1
            }
            print \$0 > \$NF"_"var1
        }
    }' \$FILENAME
    # \$NF is the last column of the file. The value in the last column is used as name
    #  print \$0 > \$NF"_"var1 print the line into the file named \$NF"_"var1
    # Warning: > append in awk, meaning that if the file already exists, a new line is added into this file
    """
}


process closest_germline {
    label 'immcantation'
    //publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{*.log}", overwrite: false
    //publishDir path: "${out_path}", mode: 'copy', pattern: "{*.tsv}", overwrite: false
    cache 'true'

    input:
    path clone_split_ch // parallelization expected
    val igblast_database_path
    val igblast_files

    output:
    path "*_germ-pass.tsv", emit: closest_ch
    path "*.log", emit: closest_log_ch

    script:
    """
    #!/bin/bash -ue
    FILENAME=\$(basename -- ${clone_split_ch}) # recover a file name without path
    cp -Lr ${clone_split_ch} "./TEMPO.tsv" # to have the hard file, not the symlink, because modifications will be performed inside
    chmod 777 TEMPO.tsv
    rm \$FILENAME # remove the initial file to avoid to send it into the channel
    cp -rp TEMPO.tsv "\$FILENAME" # -p for preserve permissions
    rm TEMPO.tsv
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a closest_germline.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a closest_germline.log
    # variables

    REPO_PATH="/usr/local/share/${igblast_database_path}" # path where the imgt_human_IGHV.fasta, imgt_human_IGHD.fasta and imgt_human_IGHJ.fasta files are in the docker container
    VDJ_FILES=\$(awk -v var1="${igblast_files}" -v var2="\${REPO_PATH}" 'BEGIN{ORS=" " ; split(var1, array1, " ") ; for (key in array1) {print var2"/"array1[key]}}')
    # end variables
    CreateGermlines.py -d \$FILENAME -g dmask --cloned -r \${VDJ_FILES} |& tee -a closest_germline.log
    """
}


process mutation_load {
    label 'immcantation'
    //publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{*.log}", overwrite: false
    //publishDir path: "${out_path}", mode: 'copy', pattern: "{*.tsv}", overwrite: false
    cache 'true'

    input:
    path closest_ch // parallelization expected


    output:
    path "*_shm-pass.tsv", emit: mutation_load_ch
    path "*.log", emit: mutation_load_log_ch

    script:
    """
    #!/bin/bash -ue
    FILENAME=\$(basename -- ${closest_ch}) # recover a file name without path
    cp -Lr ${closest_ch} "./TEMPO.tsv" # to have the hard file, not the symlink, because modifications will be performed inside
    chmod 777 TEMPO.tsv
    rm \$FILENAME # remove the initial file to avoid to send it into the channel
    cp -rp TEMPO.tsv "\$FILENAME" # -p for preserve permissions
    rm TEMPO.tsv
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a mutation_load.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a mutation_load.log
    Rscript -e '
        # Clonal assignment and germline sequences reconstruction should have been performed 
        # using the DefineClone.py and CreateGerlines.py in ChangeO
        # A "germline_alignment_d_mask" collumn should be present. 
        # If germline sequences reconstruction has been performed after clonal assignment,
        # a single germline_alignment_d_mask" consensus sequence should be present for each clone.


        args <- commandArgs(trailingOnly = TRUE)  # recover arguments written after the call of the Rscript
        tempo.arg.names <- c("file_name") # objects names exactly in the same order as in the bash code and recovered in args
        if(length(args) != length(tempo.arg.names)){
          tempo.cat <- paste0("======== ERROR IN THE NEXTFLOW EXECUTION OF THE mutation_load PROCESS\\n THE NUMBER OF ELEMENTS IN args (", length(args),") IS DIFFERENT FROM THE NUMBER OF ELEMENTS IN tempo.arg.names (", length(tempo.arg.names),")\nargs:", paste0(args, collapse = ","), "\ntempo.arg.names:", paste0(tempo.arg.names, collapse = ","))
          stop(tempo.cat)
        }
        for(i2 in 1:length(tempo.arg.names)){
          assign(tempo.arg.names[i2], args[i2])
        }

        VDJ_db <- read.table(file_name, sep = "\\t", header = TRUE)

        # Calculate R and S mutation counts
        VDJ_db <- shazam::observedMutations(VDJ_db, sequenceColumn="sequence_alignment",
                                    germlineColumn="germline_alignment_d_mask",
                                    regionDefinition=shazam::IMGT_V,
                                    frequency=FALSE, 
                                    nproc=1)
        # Calculate combined R and S mutation counts
        VDJ_db <- shazam::observedMutations(VDJ_db, sequenceColumn="sequence_alignment",
                                    germlineColumn="germline_alignment_d_mask",
                                    regionDefinition=shazam::IMGT_V,
                                    frequency=FALSE, 
                                    combine=TRUE,
                                    nproc=1)

        # Calculate R and S mutation frequencies
        VDJ_db <- shazam::observedMutations(VDJ_db, sequenceColumn="sequence_alignment",
                                    germlineColumn="germline_alignment_d_mask",
                                    regionDefinition=shazam::IMGT_V,
                                    frequency=TRUE, 
                                    nproc=1)

        # Calculate combined R and S mutation frequencies
        VDJ_db <- shazam::observedMutations(VDJ_db, sequenceColumn="sequence_alignment",
                                    germlineColumn="germline_alignment_d_mask",
                                    regionDefinition=shazam::IMGT_V,
                                    frequency=TRUE, 
                                    combine=TRUE,
                                    nproc=1)

        VDJ_db <- dplyr::arrange(VDJ_db, clone_id)
        write.table(VDJ_db, file = paste0("./", file_name, "_shm-pass.tsv"), row.names = FALSE, sep = "\\t")
    ' "\${FILENAME}" |& tee -a mutation_load.log
    # cp ./tempo_shm-pass.tsv \${FILENAME}_shm-pass.tsv
    # rm tempo_shm-pass.tsv
    """
}

process seq_name_remplacement {
    label 'r_ext'
    cache 'true'

    input:
    path mutation_load_ch // parallelization expected
    path meta_file
    val meta_name_replacement

    output:
    path "*_renamed_seq.tsv", emit: seq_name_remplacement_ch
    path '{metadata2.tsv,NULL}', includeInputs: true, emit: meta_file_ch // either take metadata2.tsv or NULL file. includeInputs: true to authorize the taking of input files, like NULL
    path "seq_name_remplacement.log", emit: seq_name_remplacement_log_ch

    script:
    """
    #!/bin/bash -ue
    FILENAME=\$(basename -- ${mutation_load_ch}) # recover a file name without path
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a seq_name_remplacement.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a seq_name_remplacement.log
    # check first that the data file does not have the second column name starting by "initial_". Otherwise, with create unproper behavior in donut
    if [[ "${meta_file}" == "NULL" ]] ; then
        rm NULL # remove the initial file to avoid to send it into the channel
        echo -n "" > NULL # new hard file that can be sent into the channel
        chmod 777 NULL
        Rscript -e '
            seq <- read.table("./${mutation_load_ch}", sep = "\\t", header = TRUE)
            if(grepl(x = names(seq)[2], pattern = "^initial_")){
                stop(paste0("\\n\\n============\\n\\nERROR IN THE seq_name_remplacement PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS \\"NULL\\", THEN THE SECOND COLUMN OF THE DATA IN THE sample_path PARAMETER CANNOT HAVE THE NAME OF THE SECOND COLUNM STARTING BY \\"initial_\\"\\n\\n============\\n\\n"), call. = FALSE)
            }
        ' |& tee -a seq_name_remplacement.log
    else
        Rscript -e '
            meta <- read.table("./${meta_file}", sep = "\\t", header = TRUE)
            if(names(meta)[1] != "Label"){
                stop(paste0("\\n\\n============\\n\\nERROR IN THE seq_name_remplacement PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS NOT \\"NULL\\", THEN THE FIRST COLUMN OF THE METADATA FILE MUST BE NAMED \\"Label\\" AND BE MADE OF THE NAMES OF THE FASTA SEQUENCES\\n\\n============\\n\\n"), call. = FALSE)
            }
        ' |& tee -a seq_name_remplacement.log
    fi
    if [[ "${meta_file}" != "NULL" && "${meta_name_replacement}" != "NULL" ]] ; then # or [[ "${meta_file}" -ne "NULL" && "${meta_name_replacement}" -ne "NULL" ]], but not !=
        Rscript -e '
            clone_id <- strsplit("${mutation_load_ch}", split = "_")[[1]][1] # because the file name starts by the clone id
            seq <- read.table("./${mutation_load_ch}", sep = "\\t", header = TRUE)
            meta <- read.table("./${meta_file}", sep = "\\t", header = TRUE)
            col_name <- "${meta_name_replacement}"
            if( ! any(names(meta) %in% col_name)){
                stop(paste0("\\n\\n============\\n\\nERROR IN THE seq_name_remplacement PROCESS OF NEXTFLOW\\nIF NOT \\"NULL\\", THE meta_name_replacement PARAMETER MUST BE A COLUMN NAME OF THE meta_path PARAMETER (METADATA FILE): ", col_name, "\\n\\n============\\n\\n"), call. = FALSE)
            }
            seq2 <- data.frame(seq[1], tempo = seq[1], seq[2:length(seq)]) # the second column is created to keep the initial sequence names, before replacement
            for(i2 in 1:nrow(meta)){
                if(sum(seq2[ , 2] %in% meta[i2, 1]) > 1){
                    stop(paste0("\\n\\n============\\n\\nERROR IN THE seq_name_remplacement PROCESS OF NEXTFLOW\\nIN THE METADATA FILE, A SEQUENCE NAME CANNOT BELONG TO SEVERAL VALUES OF THE meta_name_replacement PARAMETER COLUMN NAME OF THE meta_path PARAMETER\\nTHE METAFILE IS: ${meta_file}\\nTHE COLUM NAME IS: ", col_name, "\\nTHE PROBLEMATIC REPLACEMENT NAME IN THE METAFILE IS: ", paste(meta[i2, 1], collapse = " "), "\\n\\n============\\n\\n"), call. = FALSE)
                }else if(any(seq2[ , 2] == meta[i2, 1])){
                    if( ! (meta[i2, col_name] == "" | is.na(meta[i2, col_name]))){
                        seq2[seq2[ , 2] == meta[i2, 1], 1] <- meta[i2, col_name] # remplacement of the name in column 1
                    }
                }
            }
            names(seq2)[2] <- paste0("initial_", names(seq)[1])
            names(seq2)[1] <- names(seq)[1]
            write.table(seq2, file = paste0("./", clone_id, "_renamed_seq.tsv"), row.names = FALSE, col.names = TRUE, sep = "\\t")
            # modification of the metadata file for the correct use of ggtree::"%<+%" in tree_vizu.R that uses the column name "Label" for that 
            meta <- data.frame(meta, initial_label = meta[ , 1])
            meta[ , 1] <- meta[ , col_name]
            write.table(meta, file = "./metadata2.tsv", row.names = FALSE, col.names = TRUE, sep = "\\t")
            # end modification of the metadata file for the correct use of ggtree::"%<+%" in tree_vizu.R that uses the column name "Label" for that 
        ' |& tee -a seq_name_remplacement.log
    else
        IFS='_' read -r -a TEMPO <<< "\${FILENAME}" # string split into array
        cat ${mutation_load_ch} > ./\${TEMPO[0]}_renamed_seq.tsv |& tee -a seq_name_remplacement.log
        if [[ "${meta_file}" != "NULL" && "${meta_name_replacement}" == "NULL" ]] ; then
            cat ${meta_file} > ./metadata2.tsv |& tee -a seq_name_remplacement.log
        fi
    fi
    """
}



// add the distance to the data file to return
process file_assembly {
    label 'immcantation'
    publishDir path: "${out_path}", mode: 'copy', pattern: "{all_passed_seq.tsv}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{file_assembly.log}", overwrite: false
    cache 'true'

    input:
    path seq_name_remplacement_ch2 // no parallelization
    path distToNearest_ch
    path failed_clone_ch

    output:
    path "all_passed_seq.tsv", emit: file_assembly_ch
    path "file_assembly.log"

    script:
    """
    #!/bin/bash -ue
    Rscript -e '
        db <- read.table("${seq_name_remplacement_ch2}", header = TRUE, sep = "\\t")
        dtn <- read.table("${distToNearest_ch}", header = TRUE, sep = "\\t")

        # replace TRUE and FALSE of the read.table() conversion by the initial T and F
        # tempo.log <- sapply(db, FUN = function(x){class(x) == "logical"})
        # db[tempo.log] <- lapply(db[tempo.log], FUN = function(x){substr(as.character(x), 1, 1)})
        # end replace TRUE and FALSE of the read.table() conversion by the initial T and F

        if(file.info("${failed_clone_ch}")\$size > 0L){ # check for empty file
            dtn.fail <- read.table("${failed_clone_ch}", header = TRUE, sep = "\\t")
        }else{
            cat("\\n\\nNOTE: EMPTY ${failed_clone_ch} FILE\\n\\n")
            dtn.fail <- data.frame()
        }

        if(nrow(dtn) != (nrow(db) + nrow(dtn.fail))){
            tempo.cat <- paste0("\\n\\n========\\n\\nINTERNAL CODE ERROR IN THE NEXTFLOW EXECUTION OF THE file_assembly PROCESS\\ndtn SHOULD HAVE THE SAME NUMBER OF ROWS AS db + dtn.fail. HERE:\\ndb: ", nrow(db), "\\ndtn: ", nrow(dtn), "\\ndtn.fail: ", nrow(dtn.fail), "\\n\\n========\\n\\n")
            stop(tempo.cat)
        }
        if( ! any(c("sequence_id", "initial_sequence_id") %in% names(db))){
            tempo.cat <- paste0("\\n\\n========\\n\\nERROR IN THE NEXTFLOW EXECUTION OF THE file_assembly PROCESS\\ndb SHOULD HAVE \\"sequence_id\\" AND ALSO POTENTIALLY \\"initial_sequence_id\\" AS COLUMN NAME. HERE:\\nNAMES: ", paste(names(db), collapse = " "), "\\n\\n========\\n\\n")
            stop(tempo.cat)
        }else{
            tempo.col.name <- ifelse("initial_sequence_id" %in% names(db), "initial_sequence_id", "sequence_id")
        }
        if(all(c("sequence_id", "initial_sequence_id") %in% names(db))){
            if(all(db\$sequence_id == db\$initial_sequence_id)){
                tempo.cat <- paste0("\\n\\n========\\n\\nERROR IN THE file_assembly PROCESS OF NEXTFLOW\\nTHE meta_path AND meta_name_replacement PARAMETERS ARE NOT \\"NULL\\" BUT NO SEQUENCE NAMES HAVE BEEN REPLACED WHEN USING THE meta_name_replacement COLUMN\\n\\n========\\n\\n")
                stop(tempo.cat)
            }
        }
        if(any(is.na(match(db[, tempo.col.name], dtn\$sequence_id)))){
            tempo.cat <- paste0("\\n\\n========\\n\\nINTERNAL CODE ERROR IN THE NEXTFLOW EXECUTION OF THE file_assembly PROCESS\\nNO NA SHOULD APPEAR AT THAT STAGE WITH match()\\n\\n========\\n\\n")
            stop(tempo.cat)
        }
        if(sum(db[, tempo.col.name] %in% dtn\$sequence_id, na.rm = TRUE) != nrow(db)){
            tempo.cat <- paste0("\\n\\n========\\n\\nINTERNAL CODE ERROR IN THE NEXTFLOW EXECUTION OF THE file_assembly PROCESS\\nsum(db[, tempo.col.name] %in% dtn\$sequence_id, na.rm = TRUE) SHOULD BE EQUAL TO nrow(db)\\n\\n========\\n\\n")
            stop(tempo.cat)
        }
        if(nrow(dtn) < nrow(db)){
            tempo.cat <- paste0("\\n\\n========\\n\\nINTERNAL CODE ERROR IN THE NEXTFLOW EXECUTION OF THE file_assembly PROCESS\\ndtn CANNOT HAVE LESS ROWS THAN db. HERE:\\ndb: ", nrow(db), "\\ndtn: ", nrow(dtn), "\\n\\n========\\n\\n")
            stop(tempo.cat)
        }
        # selection of rows and ordering of dtn
        tempo.dtn <- dtn[dtn\$sequence_id %in% db[, tempo.col.name], ] # nb rows of dtn reduced to the one of db
        if(nrow(tempo.dtn) != nrow(db)){
            tempo.cat <- paste0("\\n\\n========\\n\\nINTERNAL CODE ERROR IN THE NEXTFLOW EXECUTION OF THE file_assembly PROCESS\\ntempo.dtn SHOULD HAVE THE SAME NUMBER OF ROWS AS db. HERE:\\ndb: ", nrow(db), "\\ntempo.dtn: ", nrow(tempo.dtn), "\\n\\n========\\n\\n")
            stop(tempo.cat)
        }
        tempo.dtn <- tempo.dtn[match(db[, tempo.col.name], tempo.dtn\$sequence_id), ]
        if( ! all(db[, tempo.col.name] == tempo.dtn\$sequence_id)){
            tempo.cat <- paste0("\\n\\n========\\n\\nERROR IN THE NEXTFLOW EXECUTION OF THE file_assembly PROCESS\\n", tempo.col.name, " COLUMNS OF db AND sequence_id COLUMN OF dtn SHOULD BE IDENTICAL. HERE THEY ARE\\ndb\$", tempo.col.name, ":\\n", paste0(db[, tempo.col.name], collapse = ","), "\\ndtn\$sequence_id:\\n", paste0(dtn\$sequence_id, collapse = ","), "\\n\\n========\\n\\n")
            stop(tempo.cat)
        }
        db3 <- data.frame(db, dist_nearest = tempo.dtn\$dist_nearest)
        write.table(db3, file = paste0("./all_passed_seq.tsv"), row.names = FALSE, col.names = TRUE, sep = "\\t")
    ' |& tee -a file_assembly.log
    """
}

process metadata_check { // cannot be in tree_vizu because I have to use the all_passed_seq.tsv file for the check
    label 'immcantation'
    cache 'true'

    input:
    path file_assembly_ch
    path meta_file
    val meta_name_replacement

    script:
    """
    if [[ "${meta_file}" != "NULL" && "${meta_name_replacement}" != "NULL" ]] ; then
        Rscript -e '
            df <- read.table("${file_assembly_ch}", header = TRUE, sep = "\\t")
            meta <- read.table("${meta_file}", header = TRUE, sep = "\\t")
            if(names(meta)[1] != "Label"){
                stop(paste0("\\n\\n============\\n\\nERROR IN THE metadata_check PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS NOT \\"NULL\\", THEN THE FIRST COLUMN OF THE METADATA FILE MUST BE NAMED \\"Label\\" AND BE MADE OF THE NAMES OF THE FASTA SEQUENCES\\n\\n============\\n\\n"), call. = FALSE)
            }
            if( ! any(meta\$${meta_name_replacement} %in% df[ , 1])){
                stop(paste0("\\n\\n============\\n\\nERROR IN THE metadata_check PROCESS OF NEXTFLOW\\nTHE meta_file AND meta_name_replacement PARAMETERS OF THE nextflow.config FILE ARE NOT NULL BUT NO NAME REPLACEMENT PERFORMED\\nPROBABLY THAT THE \\"Label\\" COLUMN OF THE FILE INDICATED IN THE meta_path PARAMETER IS NOT MADE OF NAMES OF FASTA FILES INDICATED IN THE sample_path PARAMETER\\nFIRST ELEMENTS OF THE LABEL COLUMN OF THE META DATA FILE ARE:\\n", paste(head(meta\$Label, 20), collapse = "\\n"), "\\nFIRST FASTA FILES NAMES ARE:\\n", paste(head(df[ , 1], 20), collapse = "\\n"), "\\n\\n\\n============\\n\\n"), call. = FALSE)
            }
        '
    fi
    """
}


process repertoire {
    label 'r_ext'
    publishDir path: "${out_path}", mode: 'copy', pattern: "{repertoire.pdf}", overwrite: false
    publishDir path: "${out_path}/png", mode: 'copy', pattern: "{*.png}", overwrite: false
    publishDir path: "${out_path}/svg", mode: 'copy', pattern: "{*.svg}", overwrite: false
    publishDir path: "${out_path}/repertoires", mode: 'copy', pattern: "{rep_*.tsv}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{repertoire.log}", overwrite: false
    cache 'true'

    input:
    val igblast_database_path
    path file_assembly_ch
    path igblast_data_check_ch
    path cute_file

    output:
    path "repertoire.pdf"
    path "*.svg"
    path "*.png"
    path "rep_*.tsv"
    path "repertoire.log"

    script:
    """
    #!/bin/bash -ue
    repertoire.R \
"${igblast_database_path}" \
"${file_assembly_ch}" \
"${igblast_data_check_ch}" \
"${cute_file}" \
"repertoire.log"
    """
}



process get_tree {
    label 'immcantation_10cpu'
    publishDir path: "${out_path}/RData", mode: 'copy', pattern: "{*_get_tree_cloneID.RData}", overwrite: false
    cache 'true'

    input:
    path seq_name_remplacement_ch // parallelization expected
    path meta_file // just to determine if metadata have been provided (TRUE means NULL) meta_file_ch not required here
    path cute_file
    val clone_nb_seq
    val tree_duplicate_seq
    val igphylm_exe_path // warning : here val and not path because we do not want the igphyml file to be imported in the work dir

    output:
    path "*_get_tree_cloneID.RData", emit: rdata_tree_ch, optional: true
    path "tree_dismissed_seq.tsv", emit: no_tree_ch
    path "seq_for_tree.tsv", emit: tree_ch
    path "tree_dismissed_clone_id.tsv", emit: no_cloneID_ch
    path "tree_clone_id.tsv", emit: cloneID_ch
    path "get_tree.log", emit: get_tree_log_ch
    //path "HLP10_tree_parameters.tsv"

    script:
    """
    #!/bin/bash -ue
    if [[ ! -s ${seq_name_remplacement_ch} ]]; then
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nEMPTY ${seq_name_remplacement_ch} FILE AS INPUT OF THE mutation_load PROCESS\\nCHECK THE mutation_load.log IN THE report FOLDER\\n\\n========\\n\\n"
        exit 1
    fi
    FILENAME=\$(basename -- ${seq_name_remplacement_ch}) # recover a file name without path
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a get_tree.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a get_tree.log
    get_tree.R \
"${seq_name_remplacement_ch}" \
"${meta_file}" \
"${clone_nb_seq}" \
"${tree_duplicate_seq}" \
"${igphylm_exe_path}" \
"${cute_file}" \
"get_tree.log"
    """
}


process tree_vizu {
    label 'r_ext'
    publishDir path: "${out_path}", mode: 'copy', pattern: "{trees.pdf}", overwrite: false
    publishDir path: "${out_path}/png", mode: 'copy', pattern: "{*.png}", overwrite: false
    publishDir path: "${out_path}/svg", mode: 'copy', pattern: "{*.svg}", overwrite: false
    publishDir path: "${out_path}/RData", mode: 'copy', pattern: "{all_trees.RData}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{tree_vizu.log}", overwrite: false
    cache 'true'

    input:
    path rdata_tree_ch2 // no more parallelization
    val tree_kind
    val clone_nb_seq
    val tree_duplicate_seq
    val tree_leaf_color
    val tree_leaf_shape
    val tree_leaf_size
    val tree_label_size
    val tree_label_hjust
    val tree_label_rigth
    val tree_label_outside
    val tree_right_margin
    val tree_legend
    path meta_file_ch
    val meta_legend
    path cute_file

    output:
    path "*.RData", optional: true
    path "*.pdf"
    path "*.png", emit: tree_vizu_ch // png plot (but sometimes empty) sustematically returned
    path "*.svg"
    path "*tree_seq_not_displayed.tsv", emit: tree_seq_not_displayed_ch
    path "tree_vizu.log"
    //path "HLP10_tree_parameters.tsv"

    script:
    """
    #!/bin/bash -ue
    tree_vizu.R \
"${tree_kind}" \
"${clone_nb_seq}" \
"${tree_duplicate_seq}" \
"${tree_leaf_color}" \
"${tree_leaf_shape}" \
"${tree_leaf_size}" \
"${tree_label_size}" \
"${tree_label_hjust}" \
"${tree_label_rigth}" \
"${tree_label_outside}" \
"${tree_right_margin}" \
"${tree_legend}" \
"${meta_file_ch}" \
"${meta_legend}" \
"${cute_file}" \
"tree_vizu.log"
    """
}



process donut {
    label 'r_ext'
    //publishDir path: "${out_path}", mode: 'copy', pattern: "{*.tsv}", overwrite: false
    //publishDir path: "${out_path}", mode: 'copy', pattern: "{*.pdf}", overwrite: false
    publishDir path: "${out_path}/png", mode: 'copy', pattern: "{*.png}", overwrite: false
    publishDir path: "${out_path}/svg", mode: 'copy', pattern: "{*.svg}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{*.log}", overwrite: false
    cache 'true'

    input:
    tuple val(kind), path(data) // 3 parallelization expected
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

    output:
    path "*.tsv", emit: donut_tsv_ch, optional: true
    path "*.pdf", emit: donut_pdf_ch, optional: true
    path "*.png"
    path "*.svg"
    path "*.log"

    script:
    """
    #!/bin/bash -ue
    FILENAME=\$(basename -- ${data}) # recover a file name without path
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a ${kind}_donut.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a ${kind}_donut.log
    donut.R \
"${data}" \
"${kind}" \
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
"${kind}_donut.log"
    """
}

process donut_assembly {
    label 'r_ext'
    publishDir path: "${out_path}", mode: 'copy', pattern: "{donuts.pdf}", overwrite: false
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
    Rscript -e '
        qpdf::pdf_combine(input = list.files(path = ".", pattern = ".pdf\$"), output = "./donuts.pdf")
    ' |& tee -a donut_assembly.log
    """
}




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
    echo -e "full .nextflow.log is in: ${launchDir}\nThe one in the result folder is not complete (miss the end)" > Log_info.txt
    """
}


//////// End Processes

//////// Workflow



workflow {

    //////// Options of nextflow run

    print("\n\nINITIATION TIME: ${workflow.start}")

    //////// end Options of nextflow run


    //////// Options of nextflow run

    // --modules (it is just for the process workflowParam)
    params.modules = "" // if --module is used, this default value will be overridden
    // end --modules (it is just for the process workflowParam)

    //////// end Options of nextflow run


    //////// Variables

    modules = params.modules // remove the dot -> can be used in bash scripts
    config_file = workflow.configFiles[0] // better to use this than config_file = file("${projectDir}/repertoire_profiler.config") because the latter is not good if -c option of nextflow run is used
    log_file = file("${launchDir}/.nextflow.log")

    //////// end Variables


    //////// Checks

    // tbi = file("${sample_path}.tbi") does not need .tbi 

    if( ! (sample_path in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN repertoire_profiler.config FILE:\n${sample_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (file(sample_path).exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN repertoire_profiler.config FILE (DOES NOT EXIST): ${sample_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    if( ! (igblast_database_path in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_database_path PARAMETER IN repertoire_profiler.config FILE:\n${igblast_database_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (igblast_organism in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_organism PARAMETER IN repertoire_profiler.config FILE:\n${igblast_organism}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (igblast_organism == "mouse" || igblast_organism == "human" || igblast_organism == "rabbit" || igblast_organism == "rat" || igblast_organism == "rhesus_monkey") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_organism PARAMETER IN repertoire_profiler.config FILE:\n${igblast_organism}\nMUST BE EITHER \"mouse\", \"human\", \"rabbit\", \"rat\" OR \"rhesus_monkey\"\n\n========\n\n"
    }
    if( ! (igblast_loci in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_loci PARAMETER IN repertoire_profiler.config FILE:\n${igblast_loci}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (igblast_loci == "ig" || igblast_loci == "tr") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_loci PARAMETER IN repertoire_profiler.config FILE:\n${igblast_loci}\nMUST BE EITHER \"ig\" OR \"tr\"\n\n========\n\n"
    }
    if( ! (igblast_aa in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_aa PARAMETER IN repertoire_profiler.config FILE:\n${igblast_aa}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (igblast_aa == "false" || igblast_aa == "true") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_aa PARAMETER IN repertoire_profiler.config FILE:\n${igblast_aa}\nMUST BE EITHER \"true\" OR \"false\"\n\n========\n\n"
    }
    if( ! (igblast_files in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_files PARAMETER IN repertoire_profiler.config FILE:\n${igblast_files}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (clone_nb_seq in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_nb_seq PARAMETER IN repertoire_profiler.config FILE:\n${clone_nb_seq}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ( ! (clone_nb_seq =~/^\d+$/)) || clone_nb_seq.toInteger() < 3 ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_nb_seq PARAMETER IN repertoire_profiler.config FILE:\n${clone_nb_seq}\nMUST BE A POSITIVE INTEGER VALUE EQUAL OR GREATER TO 3\n\n========\n\n"
    }
    if( ! (clone_model in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_model PARAMETER IN repertoire_profiler.config FILE:\n${clone_model}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (clone_model == "ham" || clone_model == "aa" || clone_model == "hh_s1f" || clone_model == "hh_s5f" || clone_model == "mk_rs1nf" || clone_model == "mk_rs5nf" || clone_model == "m1n_compat" || clone_model == "hs1f_compat") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_model PARAMETER IN repertoire_profiler.config FILE:\n${clone_model}\nMUST BE EITHER \"ham\", \"aa\", \"hh_s1f\", \"hh_s5f\", \"mk_rs1nf\", \"mk_rs5nf\", \"m1n_compat\", \"hs1f_compat\"\n\n========\n\n"
    }
    if( ! (clone_normalize in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_normalize PARAMETER IN repertoire_profiler.config FILE:\n${clone_normalize}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (clone_normalize == "len" || clone_normalize == "none") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_normalize PARAMETER IN repertoire_profiler.config FILE:\n${clone_normalize}\nMUST BE EITHER \"len\" OR \"none\"\n\n========\n\n"
    }
    if( ! (clone_distance in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_distance PARAMETER IN repertoire_profiler.config FILE:\n${clone_distance}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (clone_distance =~ /^((1)|(0)|(0\.[0-9]*))$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_distance PARAMETER IN repertoire_profiler.config FILE:\n${clone_distance}\nMUST BE A POSITIVE PROPORTION VALUE\n\n========\n\n"
    }
    if( ! (meta_path in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID meta_path PARAMETER IN repertoire_profiler.config FILE:\n${meta_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if(meta_path != "NULL"){
        if( ! (file(meta_path).exists()) ){
            error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID meta_path PARAMETER IN repertoire_profiler.config FILE (DOES NOT EXIST): ${meta_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
        }
        print("\n\nWARNING: THE 1st COLUMN OF THE METADATA FILE MUST BE NAMED \"Label\" AND MUST CONTAIN NAMES OF THE FILES IN THE FOLDER OF THE sample_path PARAMETER")
    }
    if( ! (meta_name_replacement in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID meta_name_replacement PARAMETER IN repertoire_profiler.config FILE:\n${meta_name_replacement}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (meta_legend in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID meta_legend PARAMETER IN repertoire_profiler.config FILE:\n${meta_legend}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (tree_kind in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_kind PARAMETER IN repertoire_profiler.config FILE:\n${tree_kind}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (tree_duplicate_seq in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_duplicate_seq PARAMETER IN repertoire_profiler.config FILE:\n${tree_duplicate_seq}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (tree_duplicate_seq == "TRUE" || tree_duplicate_seq == "FALSE") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_duplicate_seq PARAMETER IN repertoire_profiler.config FILE:\n${tree_duplicate_seq}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if( ! (tree_leaf_color in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_leaf_color PARAMETER IN repertoire_profiler.config FILE:\n${tree_leaf_color}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (tree_leaf_shape in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_leaf_shape PARAMETER IN repertoire_profiler.config FILE:\n${tree_leaf_shape}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (tree_leaf_shape =~  /^[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_leaf_shape PARAMETER IN repertoire_profiler.config FILE:\n${tree_leaf_shape}\nMUST BE A POSITIVE INTEGER VALUE\n\n========\n\n"
    }
    if( ! (tree_leaf_size in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_leaf_size PARAMETER IN repertoire_profiler.config FILE:\n${tree_leaf_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (tree_leaf_size =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_leaf_size PARAMETER IN repertoire_profiler.config FILE:\n${tree_leaf_size}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (tree_label_size in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_label_size PARAMETER IN repertoire_profiler.config FILE:\n${tree_label_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (tree_label_size =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_label_size PARAMETER IN repertoire_profiler.config FILE:\n${tree_label_size}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (tree_label_hjust in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_label_hjust PARAMETER IN repertoire_profiler.config FILE:\n${tree_label_hjust}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (tree_label_hjust =~  /^\-{0,1}[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_label_hjust PARAMETER IN repertoire_profiler.config FILE:\n${tree_label_hjust}\nMUST BE A NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (tree_label_rigth in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_label_rigth PARAMETER IN repertoire_profiler.config FILE:\n${tree_label_rigth}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (tree_label_rigth == "TRUE" || tree_label_rigth == "FALSE") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_label_rigth PARAMETER IN repertoire_profiler.config FILE:\n${tree_label_rigth}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if( ! (tree_label_outside in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_label_outside PARAMETER IN repertoire_profiler.config FILE:\n${tree_label_outside}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (tree_label_outside == "TRUE" || tree_label_outside == "FALSE") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_label_outside PARAMETER IN repertoire_profiler.config FILE:\n${tree_label_outside}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if( ! (tree_right_margin in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_right_margin PARAMETER IN repertoire_profiler.config FILE:\n${tree_right_margin}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (tree_right_margin =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_right_margin PARAMETER IN repertoire_profiler.config FILE:\n${tree_right_margin}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (tree_legend in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_legend PARAMETER IN repertoire_profiler.config FILE:\n${tree_legend}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (tree_legend == "TRUE" || tree_legend == "FALSE") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_legend PARAMETER IN repertoire_profiler.config FILE:\n${tree_legend}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }

    if( ! (donut_palette in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_palette PARAMETER IN repertoire_profiler.config FILE:\n${donut_palette}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (donut_hole_size in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_hole_size PARAMETER IN repertoire_profiler.config FILE:\n${donut_hole_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_hole_size =~  /^((1)|(0)|(0\.[0-9]*))$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_hole_size PARAMETER IN repertoire_profiler.config FILE:\n${donut_hole_size}\nMUST BE A POSITIVE PROPORTION VALUE\n\n========\n\n"
    }
    if( ! (donut_hole_text in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_hole_text PARAMETER IN repertoire_profiler.config FILE:\n${tree_legend}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_hole_text == "TRUE" || donut_hole_text == "FALSE") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_hole_text PARAMETER IN repertoire_profiler.config FILE:\n${donut_hole_text}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if( ! (donut_hole_text_size in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_hole_text_size PARAMETER IN repertoire_profiler.config FILE:\n${donut_hole_text_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_hole_text_size =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_hole_text_size PARAMETER IN repertoire_profiler.config FILE:\n${donut_hole_text_size}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (donut_border_color in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_border_color PARAMETER IN repertoire_profiler.config FILE:\n${donut_border_color}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (donut_border_size in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_border_size PARAMETER IN repertoire_profiler.config FILE:\n${donut_border_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (donut_annotation_distance in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_distance PARAMETER IN repertoire_profiler.config FILE:\n${donut_annotation_distance}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_annotation_distance =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_distance PARAMETER IN repertoire_profiler.config FILE:\n${donut_annotation_distance}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (donut_annotation_size in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_size PARAMETER IN repertoire_profiler.config FILE:\n${donut_annotation_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_annotation_size =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_size PARAMETER IN repertoire_profiler.config FILE:\n${donut_annotation_size}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (donut_annotation_force in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_force PARAMETER IN repertoire_profiler.config FILE:\n${donut_annotation_force}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_annotation_force =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_force PARAMETER IN repertoire_profiler.config FILE:\n${donut_annotation_force}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (donut_annotation_force_pull in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_force_pull PARAMETER IN repertoire_profiler.config FILE:\n${donut_annotation_force_pull}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_annotation_force_pull =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_force_pull PARAMETER IN repertoire_profiler.config FILE:\n${donut_annotation_force_pull}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (donut_legend_width in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_width PARAMETER IN repertoire_profiler.config FILE:\n${donut_legend_width}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_legend_width =~  /^((1)|(0)|(0\.[0-9]*))$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_width PARAMETER IN repertoire_profiler.config FILE:\n${donut_legend_width}\nMUST BE A POSITIVE PROPORTION VALUE\n\n========\n\n"
    }
    if( ! (donut_legend_text_size in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_text_size PARAMETER IN repertoire_profiler.config FILE:\n${donut_legend_text_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_legend_text_size =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_text_size PARAMETER IN repertoire_profiler.config FILE:\n${donut_legend_text_size}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (donut_legend_box_size in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_box_size PARAMETER IN repertoire_profiler.config FILE:\n${donut_legend_box_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_legend_box_size =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_box_size PARAMETER IN repertoire_profiler.config FILE:\n${donut_legend_box_size}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (donut_legend_box_space in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_box_space PARAMETER IN repertoire_profiler.config FILE:\n${donut_legend_box_space}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_legend_box_space =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_box_space PARAMETER IN repertoire_profiler.config FILE:\n${donut_legend_box_space}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (donut_legend_limit in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_limit PARAMETER IN repertoire_profiler.config FILE:\n${donut_legend_limit}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (donut_legend_width ==  "NULL") ){
        if( ! (donut_legend_width =~  /^((1)|(0)|(0\.[0-9]*))$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_limit PARAMETER IN repertoire_profiler.config FILE:\n${donut_legend_limit}\nMUST BE A POSITIVE PROPORTION VALUE IF NOT \"NULL\"\n\n========\n\n"
        }
    }
    if( ! (cute_path in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID cute_path PARAMETER IN repertoire_profiler.config FILE:\n${cute_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (file(cute_path).exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID cute_path PARAMETER IN repertoire_profiler.config FILE (DOES NOT EXIST): ${cute_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }

    if( ! (igphylm_exe_path in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igphylm_exe_path PARAMETER IN repertoire_profiler.config FILE:\n${igphylm_exe_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }


    // below : those variable are already used in the config file. Thus, to late to check them. And not possible to check inside the config file
    // out_ini
    print("\n\nRESULT DIRECTORY: ${out_path}")
    //print("\n\nWARNING: PARAMETERS ALREADY INTERPRETED IN THE .config FILE:")
    //print("    system_exec: ${system_exec}")
    //print("    out_path: ${out_path_ini}")
    if("${system_exec}" != "local"){
        print("    queue: ${queue}")
        print("    qos: ${qos}")
        print("    add_options: ${add_options}")
    }
    if(igblast_files =~ /^.*IGKV.*$/){
        print("\n\nWARNING:\nLIGHT CHAIN DETECTED IN THE igblast_files parameter.\nBUT CLONAL GROUPING IS GENERALLY RESTRICTED TO HEAVY CHAIN SEQUENCES, AS THE DIVERSITY OF LIGHT CHAINS IS NOT SUFFICIENT TO DISTINGUISH CLONES WITH REASONABLE CERTAINTY")
    }
    print("\n\nWARNING:\nTO MAKE THE REPERTOIRES, THE SCRIPT CURRENTLY TAKES THE FIRST ANNOTATION OF THE IMGT ANNOTATION IF SEVERAL ARE PRESENTS IN THE v_call OR j_call COLUMN OF THE all_passed_seq.tsv FILE")
    print("\n\n")



    //////// end Checks


    //////// Variable modification


    //////// end Variable modification


    //////// Channels

    fs_ch = Channel.fromPath("${sample_path}/*.*", checkIfExists: false) // in channel because many files

    //////// end Channels


    //////// files import

    meta_file = file(meta_path) // in variable because a single file. If "NULL", will create a empty file, present in work folders, but that cannot be correctly linked. Thus, if the file has to be redirected into a channel inside a process, it will not work. Thus, in the first process using meta_file, I hard copy the NULL file if required (see below)
    cute_file = file(cute_path) // in variable because a single file

    //////// end files import


    //////// Main

    workflowParam(
        modules
    )


    igblast_data_check(
        igblast_database_path, 
        igblast_files
    )

    igblast(
        fs_ch, 
        igblast_database_path, 
        igblast_files, 
        igblast_organism, 
        igblast_loci, 
        igblast_aa
    )
    igblast.out.tsv_ch1.count().subscribe{ n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\n0 ANNOTATION SUCCEEDED BY THE igblast PROCESS\n\nCHECK THAT THE igblast_organism, igblast_loci AND igblast_files ARE CORRECTLY SET IN THE repertoire_profiler.config FILE\n\n========\n\n"}}
    tsv_ch2 = igblast.out.tsv_ch1.collectFile(name: "all_igblast_seq.tsv", skip: 1, keepHeader: true)
    igblast.out.aligned_seq_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\n0 ALIGNED SEQ FILES RETURNED BY THE igblast PROCESS\n\n========\n\n"}}
    igblast.out.aligned_seq_ch.collectFile(name: "igblast_aligned_seq.tsv").subscribe{it -> it.copyTo("${out_path}")}
    igblast.out.unaligned_seq_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\n0 UNALIGNED SEQ FILES RETURNED BY THE igblast PROCESS\n\n========\n\n"}}
    igblast.out.unaligned_seq_ch.collectFile(name: "igblast_unaligned_seq.tsv").subscribe{it -> it.copyTo("${out_path}")}
    nb1 = igblast.out.aligned_seq_ch.collectFile().countLines()
    nb2 = igblast.out.unaligned_seq_ch.collectFile().countLines()
    fs_ch.count().combine(nb1).combine(nb2).subscribe{n,n1,n2 -> if(n != n1 + n2){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nTHE NUMBER OF FILES IN THE igblast_aligned_seq.tsv AND igblast_unaligned_seq.tsv IS NOT EQUAL TO THE NUMBER OF SUBMITTED FASTA FILES\n========\n\n"}}
    igblast.out.log_ch.collectFile(name: "igblast_report.log").subscribe{it -> it.copyTo("${out_path}/reports")}


    parseDb_filtering(
        tsv_ch2,
        igblast_aa
    )
    // parseDb_filtering.out.unproductive_ch.count().subscribe{n -> if ( n == 0 ){print "\n\nWARNING: EMPTY unproductive_seq.tsv FILE RETURNED FOLLOWING THE parseDb_filtering PROCESS\n\n"}else{it -> it.copyTo("${out_path}/unproductive_seq.tsv")}} // see https://www.nextflow.io/docs/latest/script.html?highlight=copyto#copy-files
    parseDb_filtering.out.select_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE parseDb_filtering PROCESS\n\n========\n\n"}}



    translation(
        parseDb_filtering.out.select_ch,
        igblast_aa
    )
    translation.out.translation_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE translation PROCESS\n\n========\n\n"}}



    distToNearest(
        translation.out.translation_ch,
        igblast_aa,
        clone_model,
        clone_normalize
    )
    distToNearest.out.distToNearest_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE distToNearest PROCESS\n\n========\n\n"}}



    distance_hist(
        distToNearest.out.distToNearest_ch,
        cute_file, 
        clone_model,
        clone_normalize,
        clone_distance
    )
    distance_hist.out.distance_hist_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE distance_hist PROCESS\n\n========\n\n"}}



    histogram_assembly(
        distance_hist.out.histogram_pdf_ch.collect()
    )



    clone_assignment(
        translation.out.translation_ch,
        clone_model,
        clone_normalize,
        clone_distance
    )
    // clone_assignment.out.clone_ch.view()



    split_by_clones(
        clone_assignment.out.clone_ch
    )



    closest_germline(
        split_by_clones.out.clone_split_ch.flatten(), // flatten split the list into several objects (required for parallelization)
        igblast_database_path, 
        igblast_files
    )
    closest_germline.out.closest_log_ch.collectFile(name: "closest_germline.log").subscribe{it -> it.copyTo("${out_path}/reports")} // 



    mutation_load(
        closest_germline.out.closest_ch
    )
    mutation_load.out.mutation_load_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE mutation_load PROCESS\n\n========\n\n"}}
    mutation_load.out.mutation_load_log_ch.collectFile(name: "mutation_load.log").subscribe{it -> it.copyTo("${out_path}/reports")} // 



    seq_name_remplacement(
        mutation_load.out.mutation_load_ch,
        meta_file,
        meta_name_replacement
    )
    seq_name_remplacement.out.seq_name_remplacement_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE seq_name_remplacement PROCESS\n\n========\n\n"}}
    seq_name_remplacement_ch2 = seq_name_remplacement.out.seq_name_remplacement_ch.collectFile(name: "replacement.tsv", skip: 1, keepHeader: true)
    //seq_name_remplacement_ch2.subscribe{it -> it.copyTo("${out_path}")}
    // tuple_seq_name_remplacement = new Tuple("all", seq_name_remplacement_ch2) # warning: this is not a channel but a variable now
    seq_name_remplacement.out.seq_name_remplacement_log_ch.collectFile(name: "seq_name_remplacement.log").subscribe{it -> it.copyTo("${out_path}/reports")}

    seq_name_remplacement.out.meta_file_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE seq_name_remplacement PROCESS\n\n========\n\n"}}



    file_assembly(
        seq_name_remplacement_ch2, 
        distToNearest.out.distToNearest_ch, 
        clone_assignment.out.failed_clone_ch
    )
    file_assembly.out.file_assembly_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE file_assembly PROCESS\n\n========\n\n"}}


    metadata_check(
        file_assembly.out.file_assembly_ch,
        meta_file,
        meta_name_replacement
    )


    repertoire(
        igblast_database_path,
        file_assembly.out.file_assembly_ch,
        igblast_data_check.out.igblast_data_check_ch,
        cute_file
    )



    get_tree(
        seq_name_remplacement.out.seq_name_remplacement_ch,
        meta_file, // first() because get_tree process is a parallele one and because meta_file is single
        cute_file, 
        clone_nb_seq,
        tree_duplicate_seq,
        igphylm_exe_path
    )
    get_tree.out.rdata_tree_ch.count().subscribe { n -> if ( n == 0 ){print("\n\nWARNING: EMPTY OUTPUT FOLLOWING THE get_tree PROCESS -> NO TREE RETURNED\n\n")}}
    rdata_tree_ch2 = get_tree.out.rdata_tree_ch.collect()
    //rdata_tree_ch2.view()

    get_tree.out.no_tree_ch.count().subscribe { n -> if ( n == 0 ){print("\n\nWARNING: ALL SEQUENCES IN TREES FOLLOWING THE get_tree PROCESS -> EMPTY tree_dismissed_seq.tsv FILE RETURNED\n\n")}}
    no_tree_ch2 = get_tree.out.no_tree_ch.collectFile(name: "tree_dismissed_seq.tsv", skip: 1, keepHeader: true)
    no_tree_ch2.subscribe{it -> it.copyTo("${out_path}")}

    get_tree.out.tree_ch.count().subscribe { n -> if ( n == 0 ){print("\n\nWARNING: NO SEQUENCES IN TREES FOLLOWING THE get_tree PROCESS -> EMPTY tree_seq.tsv FILE RETURNED\n\n")}}
    tree_ch2 = get_tree.out.tree_ch.collectFile(name: "tree_seq.tsv", skip: 1, keepHeader: true)
    tree_ch2.subscribe{it -> it.copyTo("${out_path}")}

    get_tree.out.no_cloneID_ch.count().subscribe { n -> if ( n == 0 ){print("\n\nWARNING: ALL SEQUENCES IN CLONAL GROUP FOLLOWING THE get_tree PROCESS -> EMPTY tree_dismissed_clone_id.tsv FILE RETURNED\n\n")}}
    no_cloneID_ch2 = get_tree.out.no_cloneID_ch.collectFile(name: "tree_dismissed_clone_id.tsv")
    no_cloneID_ch2.subscribe{it -> it.copyTo("${out_path}")}

    get_tree.out.cloneID_ch.count().subscribe { n -> if ( n == 0 ){print("\n\nWARNING: NO CLONAL GROUP FOLLOWING THE get_tree PROCESS -> EMPTY tree_clone_id.tsv and trees.pdf FILES RETURNED\n\n")}}
    cloneID_ch2 = get_tree.out.cloneID_ch.collectFile(name: "tree_clone_id.tsv")
    cloneID_ch2.subscribe{it -> it.copyTo("${out_path}")}

    get_tree.out.get_tree_log_ch.collectFile(name: "get_tree.log").subscribe{it -> it.copyTo("${out_path}/reports")} // 



    tree_vizu(
        rdata_tree_ch2,
        tree_kind,
        clone_nb_seq,
        tree_duplicate_seq,
        tree_leaf_color,
        tree_leaf_shape,
        tree_leaf_size,
        tree_label_size,
        tree_label_hjust,
        tree_label_rigth,
        tree_label_outside,
        tree_right_margin,
        tree_legend,
        seq_name_remplacement.out.meta_file_ch.first(), // first() because seq_name_remplacement process is a parallele one. Possible that it prevent the cache to work, depending on the order in 
        meta_legend,
        cute_file
    )
   tree_vizu.out.tree_vizu_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE donut_assembly PROCESS\n\n========\n\n"}}
    tree_vizu.out.tree_seq_not_displayed_ch.count().subscribe { n -> if ( n == 0 ){print("\n\nWARNING: -> NO tree_seq_not_displayed.tsv FILE RETURNED\n\n")}}
    tree_seq_not_displayed_ch2 = tree_vizu.out.tree_seq_not_displayed_ch.flatten().collectFile(name: "tree_seq_not_displayed.tsv", skip: 1, keepHeader: true) // flatten split the list into several objects which is required by collectFile()
    tree_seq_not_displayed_ch2.subscribe{it -> it.copyTo("${out_path}")}



    tempo1_ch = Channel.of("all", "annotated", "tree") // 1 channel with 3 values (not list)
    tempo2_ch = seq_name_remplacement_ch2.mix(seq_name_remplacement_ch2.mix(tree_ch2)) // 1 channel with 3 paths (do not use flatten() -> not list)
    tempo3_ch = tempo1_ch.merge(tempo2_ch) // 3 lists


    donut(
        tempo3_ch,
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
        cute_file
    )
    donut.out.donut_pdf_ch.count().subscribe { n -> if ( n == 0 ){print("\n\nWARNING: EMPTY OUTPUT FOLLOWING THE donut PROCESS -> NO DONUT RETURNED\n\n")}}
    donut_pdf_ch2 = donut.out.donut_pdf_ch.collect()

    donut.out.donut_tsv_ch.count().subscribe { n -> if ( n == 0 ){print("\n\nWARNING: -> NO donut_stats.tsv FILE RETURNED\n\n")}}
    donut_tsv_ch2 = donut.out.donut_tsv_ch.collectFile(name: "donut_stats.tsv", skip: 1, keepHeader: true)
    donut_tsv_ch2.subscribe{it -> it.copyTo("${out_path}")}



    donut_assembly(
        donut_pdf_ch2
    )
   donut_assembly.out.donut_assembly_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE donut_assembly PROCESS\n\n========\n\n"}}



    backup(
        config_file, 
        log_file
    )
}

    //////// end Main


//////// end Processes
