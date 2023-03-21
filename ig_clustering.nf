nextflow.enable.dsl=2
/*
#########################################################################
##                                                                     ##
##     ig_clustering.nf                                                ##
##                                                                     ##
##     Gael A. Millot                                                  ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################
*/



//////// Processes


process workflowParam { // create a file with the workflow parameters in out_path
    label 'bash' // see the withLabel: bash in the nextflow config file 
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


process repertoire_names {
    // cannot be repertoire_names outside of process because the files are present in the docker
    publishDir path: "${out_path}", mode: 'copy', pattern: "{*.tsv}", overwrite: false
    label 'immcantation' // see the withLabel: bash in the nextflow config file 
    cache 'true'

    input:
    val igblast_database_path
    val igblast_files

    output:
    path "*.tsv", emit: repertoire_names_ch


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
    label 'immcantation' // see the withLabel: bash in the nextflow config file 
    //publishDir path: "${out_path}", mode: 'copy', overwrite: false
    cache 'true'

    input:
    path fs_ch // parallelization expected (for each fasta file)
    val igblast_database_path
    val igblast_files
    val igblast_organism
    val igblast_loci

    output:
    path "*.tsv", emit: tsv_ch1, optional: true
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
    # end variables

    # checks
    if [[ ! "\${FILE_EXTENSION}" =~ fasta|fa|fna|txt|seq ]] ; then
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nINVALID FILE EXTENSION IN THE sample_path PARAMETER OF THE ig_clustering.config FILE:\\n${fs_ch}\\n\${FILENAME}\\nMUST BE fasta|fs_ch|txt|seq\\n\\n========\\n\\n"
        exit 1
    fi
    sed 's/\\r\$//' ${fs_ch} > tempo_file.fasta #remove carriage returns
    awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' tempo_file.fasta > \${FILE}.fa # remove \\n in the middle of the sequence # \${FILENAME}.fa is a trick to do not use ${fs_ch} and modify the initial file due to the link in the nextflow work folder
    TEMPO=\$(wc -l \${FILE}.fa | cut -f1 -d' ')
    if read -n 1 char <"\${FILE}.fa"; [[ \$char != ">" || \$TEMPO != 2 ]]; then
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nINVALID FASTA FILE IN THE sample_path OF THE ig_clustering.config FILE:\\n${fs_ch}\\nMUST BE A FASTA FILE (\'>\' AS FIRST CHARATER) MADE OF A SINGLE SEQUENCE\\n\\n========\\n\\n"
        exit 1
    fi
    # end checks

    # Alignment <-> annotate sequence using VDJ info
    # See https://changeo.readthedocs.io/en/stable/tools/AssignGenes.html for the details
    AssignGenes.py igblast -s \${FILE}.fa -b /usr/local/share/igblast --organism ${igblast_organism} --loci ${igblast_loci} --format blast |& tee -a igblast_report.log

    # convert to tsv
    # Also convert data from the web interface IMGT/HighV-QUEST
    MakeDb.py igblast -i ./\${FILE}_igblast.fmt7 -s ./\${FILE}.fa -r \${VDJ_FILES} --extended |& tee -a igblast_report.log
    
    # printing if no tsv file made
    if [[ ! -f ./\${FILE}_igblast_db-pass.tsv ]] ; then
        echo -e "MakeDb.py igblast FAIL FOR \${FILENAME}" |& tee -a igblast_report.log
    # else
        # echo -e "" |& tee -a igblast_report.log
    fi
    """
    // write ${} between "" to make a single argument when the variable is made of several values separated by a space. Otherwise, several arguments will be considered
}


process parseDb_filtering {
    label 'immcantation' // see the withLabel: bash in the nextflow config file 
    publishDir path: "${out_path}", mode: 'copy', pattern: "{*_productive-F.tsv}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{ParseDb_filtering.log}", overwrite: false
    cache 'true'

    input:
    path tsv_ch2 // no more parallelization

    output:
    path "*_parse-select.tsv", emit: select_ch
    path "*_productive-F.tsv", emit: unproductive_ch, optional: true
    path "*.log"

    script:
    """
    #!/bin/bash -ue
    if [[ ! -s ${tsv_ch2} ]]; then
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nEMPTY FILE GENERATED BY THE igblast PROCESS\\nCHECK THE ParseDb_filtering.log FILE IN THE report FOLDER\\n\\n========\\n\\n"
        exit 1
    else
        ParseDb.py select -d ${tsv_ch2} -f productive -u T |& tee -a ParseDb_filtering.log
        ParseDb.py split -d ${tsv_ch2} -f productive |& tee -a ParseDb_filtering.log
    fi
    """
}


process clone_assignment {
    label 'immcantation' // see the withLabel: bash in the nextflow config file 
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{*.log}", overwrite: false
    //publishDir path: "${out_path}", mode: 'copy', pattern: "{*_clone-pass.tsv}", overwrite: false
    cache 'true'

    input:
    path select_ch
    val clone_distance

    output:
    path "*_clone-pass.tsv", emit: clone_ch
    path "*.log"

    script:
    """
    #!/bin/bash -ue
    if [[ ! -s ${select_ch} ]]; then
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nEMPTY FILE GENERATED BY THE parseDb_filtering PROCESS\\nCHECK THE clone_assignment.log AND *_productive-F.tsv FILES IN THE report FOLDER\\n\\n========\\n\\n"
        exit 1
    else
        DefineClones.py -d ${select_ch} --act set --model ham --norm len --dist ${clone_distance} |& tee -a clone_assignment.log
    fi
    """
}

process split_by_clones { // split the file into multiple files according to the clone_id column
    label 'immcantation' // see the withLabel: bash in the nextflow config file 
    cache 'true'

    input:
    path clone_ch

    output:
    path "*clone-pass.tsv", emit: clone_split_ch // multiple files -> parall expected

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
        if (NR == 1){header=\$0 ; next}else{print header > \$NF"_"var1 ; print \$0 > \$NF"_"var1}
    }' \$FILENAME
    # \$NF is the last column of the file. The value in the last column is used as name
    #  print \$0 > \$NF"_"var1 print the line into the file named \$NF"_"var1
    # Warning: > append in awk, meaning that if the file already exists, a new line is added into this file
    """
}


process closest_germline {
    label 'immcantation' // see the withLabel: bash in the nextflow config file 
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
    # variables
    REPO_PATH="/usr/local/share/${igblast_database_path}" # path where the imgt_human_IGHV.fasta, imgt_human_IGHD.fasta and imgt_human_IGHJ.fasta files are in the docker container
    VDJ_FILES=\$(awk -v var1="${igblast_files}" -v var2="\${REPO_PATH}" 'BEGIN{ORS=" " ; split(var1, array1, " ") ; for (key in array1) {print var2"/"array1[key]}}')
    # end variables
    CreateGermlines.py -d \$FILENAME -g dmask --cloned -r \${VDJ_FILES} |& tee -a closest_germline.log
    """
}


process mutation_load {
    label 'immcantation' // see the withLabel: bash in the nextflow config file 
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
    Rscript -e '
        # Clonal assignment and germline sequences reconstruction should have been performed 
        # using the DefineClone.py and CreateGerlines.py in ChangeO
        # A "germline_alignment_d_mask" collumn should be present. 
        # If germline sequences reconstruction has been performed after clonal assignment,
        # a single germline_alignment_d_mask" consensus sequence should be present for each clone.


        args <- commandArgs(trailingOnly = TRUE)  # recover arguments written after the call of the Rscript
        tempo.arg.names <- c("file_name") # objects names exactly in the same order as in the bash code and recovered in args
        if(length(args) != length(tempo.arg.names)){
          tempo.cat <- paste0("======== ERROR: THE NUMBER OF ELEMENTS IN args (", length(args),") IS DIFFERENT FROM THE NUMBER OF ELEMENTS IN tempo.arg.names (", length(tempo.arg.names),")\nargs:", paste0(args, collapse = ","), "\ntempo.arg.names:", paste0(tempo.arg.names, collapse = ","))
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

process get_tree {
    label 'immcantation' // see the withLabel: bash in the nextflow config file 
    publishDir path: "${out_path}/RData", mode: 'copy', pattern: "{*.RData}", overwrite: false
    cache 'true'

    input:
    path mutation_load_ch
    val nb_seq_per_clone
    val tree_duplicate_seq

    output:
    path "*.RData", emit: rdata_tree_ch, optional: true
    path "tree_dismissed_seq.tsv", emit: no_tree_ch, optional: true
    path "seq_for_tree.tsv", emit: tree_ch, optional: true
    path "tree_dismissed_clone_id.tsv", emit: no_cloneID_ch, optional: true
    path "tree_clone_id.tsv", emit: cloneID_ch, optional: true
    path "get_tree.log", emit: get_tree_log_ch
    //path "HLP10_tree_parameters.tsv"

    script:
    """
    #!/bin/bash -ue
    if [[ ! -s ${mutation_load_ch} ]]; then
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nEMPTY ${mutation_load_ch} FILE AS INPUT OF THE mutation_load PROCESS\\nCHECK THE mutation_load.log IN THE report FOLDER\\n\\n========\\n\\n"
        exit 1
    fi
    FILENAME=\$(basename -- ${mutation_load_ch}) # recover a file name without path
    echo -e "\${FILENAME}:" |& tee -a get_tree.log
    LINE_NB=\$((\$(cat ${mutation_load_ch} | wc -l) - 1))
    if [[ "\$LINE_NB" -ge "${nb_seq_per_clone}" ]] ; then # the -gt operator can compare strings and means "greater than"
        cat ${mutation_load_ch} > seq_for_tree.tsv
        Rscript -e '
            db <- alakazam::readChangeoDb("${mutation_load_ch}")
            db <- tibble::add_column(db, collapsed = db[[1]]) # add a column to detected collapsed identical sequences names
            clones <- dowser::formatClones(
                data = db, 
                seq = "sequence_alignment", 
                clone = "clone_id", 
                minseq=1, 
                heavy = NULL, 
                text_fields = "collapsed", # column collapsed
                dup_singles = FALSE, # Duplicate sequences in singleton clones to include them as trees? Always use FALSE. Otherwise, it will create an artificial duplicated sequences in thr tree building so that we will have two instead of one sequence in the tree, after identical seq removal. See https://rdrr.io/cran/dowser/src/R/Clones.R
                trait = if(${tree_duplicate_seq}){names(db)[1]}else{NULL} # control the removal of identical sequences but different sequence name
            )
            trees <- dowser::getTrees(clones, build="igphyml", exec="/usr/local/share/igphyml/src/igphyml", nproc = 10)
            # assign(paste0("c", trees\$clone_id, "_trees"), trees)
            save(list = c("trees", "db", "clones"), file = paste0("./", trees\$clone_id, "_trees.RData"))
            write.table(trees\$clone_id, file = paste0("./tree_clone_id.tsv"), row.names = FALSE, col.names = FALSE)
        ' |& tee -a get_tree.log
    else
        echo -e "LESS THAN ${nb_seq_per_clone} SEQUENCES FOR THE CLONAL GROUP: NO TREE COMPUTED" |& tee -a get_tree.log 
        IFS='_' read -r -a TEMPO <<< "\${FILENAME}" # string split into array
        echo \${TEMPO[0]} > tree_dismissed_clone_id.tsv
        cat ${mutation_load_ch} > tree_dismissed_seq.tsv
    fi
    """
}


process tree_vizu {
    label 'r_ext' // see the withLabel: bash in the nextflow config file 
    publishDir path: "${out_path}", mode: 'copy', pattern: "{trees.pdf}", overwrite: false
    publishDir path: "${out_path}/png", mode: 'copy', pattern: "{*.png}", overwrite: false
    publishDir path: "${out_path}/svg", mode: 'copy', pattern: "{*.svg}", overwrite: false
    publishDir path: "${out_path}/RData", mode: 'copy', pattern: "{all_trees.RData}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{tree_vizu.log}", overwrite: false
    cache 'true'

    input:
    path rdata_tree_ch2 // no more parallelization
    val tree_kind
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
    path meta_file
    val tree_meta_path_names
    path cute_file

    output:
    path "*.RData", optional: true
    path "*.pdf", optional: true
    path "*.png", optional: true
    path "*.svg", optional: true
    path "*.tsv", emit: tree_seq_not_displayed_ch, optional: true
    path "tree_vizu.log"
    //path "HLP10_tree_parameters.tsv"

    script:
    """
    #!/bin/bash -ue
    tree_vizu.R \
"${tree_kind}" \
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
"${meta_file}" \
"${tree_meta_path_names}" \
"${cute_file}" \
"tree_vizu.log"
    """
    // Warning: $workflow.projectDir/bin/ is the only way to have the execution rights of a .R file in the bin directory when the gitlab repo is pulled into /pasteur/sonic/homes/gmillot/.nextflow/assets/. See https://github.com/nextflow-io/nextflow/issues/698. Otherwise, the following message can appear: Fatal error: cannot open file '/pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot/14985_loot/bin/plot_fivep_filtering.R': No such file or directory
}



process donut {
    label 'r_ext' // see the withLabel: bash in the nextflow config file 
    //publishDir path: "${out_path}", mode: 'copy', pattern: "{*.tsv}", overwrite: false
    //publishDir path: "${out_path}", mode: 'copy', pattern: "{*.pdf}", overwrite: false
    publishDir path: "${out_path}/png", mode: 'copy', pattern: "{*.png}", overwrite: false
    publishDir path: "${out_path}/svg", mode: 'copy', pattern: "{*.svg}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{*.log}", overwrite: false
    cache 'true'

    input:
    tuple val(kind), path(data) // 2 parallelization expected
    val donut_hole_size
    val donut_colors
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
    donut.R \
"${data}" \
"${kind}" \
"${donut_hole_size}" \
"${donut_colors}" \
"${cute_file}" \
"${kind}_donut.log"
    """
}

process donut_assembly {
    label 'r_ext' // see the withLabel: bash in the nextflow config file 
    publishDir path: "${out_path}", mode: 'copy', pattern: "{donuts.pdf}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{donut_assembly.log}", overwrite: false
    cache 'true'

    input:
    path donut_pdf_ch2

    output:
    path "donuts.pdf"
    path "donut_assembly.log"

    script:
    """
    #!/bin/bash -ue
    Rscript -e '
        qpdf::pdf_combine(input = list.files(path = ".", pattern = ".pdf\$"), output = "./donuts.pdf")
    ' |& tee -a donut_assembly.log
    """
}

process repertoire {
    label 'immcantation' // see the withLabel: bash in the nextflow config file 
    publishDir path: "${out_path}/RData", mode: 'copy', pattern: "{*.RData}", overwrite: false
    cache 'true'

    input:
    path mutation_load_ch2

    output:
    path "tree_dismissed_seq.tsv", emit: no_tree_ch, optional: true

    script:
    """
    #!/bin/bash -ue
    Rscript -e '
        qpdf::pdf_combine(input = list.files(path = ".", pattern = ".pdf\$"), output = "./donuts.pdf")
    ' |& tee -a donut_assembly.log
    """
}



process backup {
    label 'bash' // see the withLabel: bash in the nextflow config file 
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
    config_file = workflow.configFiles[0] // better to use this than config_file = file("${projectDir}/ig_clustering.config") because the latter is not good if -c option of nextflow run is used
    log_file = file("${launchDir}/.nextflow.log")

    //////// end Variables


    //////// Checks

    // tbi = file("${sample_path}.tbi") does not need .tbi 

    if( ! sample_path in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN ig_clustering.config FILE:\n${sample_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! file(sample_path).exists()){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN ig_clustering.config FILE (DOES NOT EXIST): ${sample_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    if( ! igblast_database_path in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_database_path PARAMETER IN ig_clustering.config FILE:\n${igblast_database_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! igblast_organism in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_organism PARAMETER IN ig_clustering.config FILE:\n${igblast_organism}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (igblast_organism == "mouse" || igblast_organism == "human" || igblast_organism == "rabbit" || igblast_organism == "rat" || igblast_organism == "rhesus_monkey")){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_organism PARAMETER IN ig_clustering.config FILE:\n${igblast_organism}\nMUST BE EITHER \"mouse\", \"human\", \"rabbit\", \"rat\" OR \"rhesus_monkey\"\n\n========\n\n"
    }
    if( ! igblast_loci in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_loci PARAMETER IN ig_clustering.config FILE:\n${igblast_loci}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (igblast_loci == "ig" || igblast_loci == "tr")){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_loci PARAMETER IN ig_clustering.config FILE:\n${igblast_loci}\nMUST BE EITHER \"ig\" OR \"tr\"\n\n========\n\n"
    }
    if( ! igblast_files in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_files PARAMETER IN ig_clustering.config FILE:\n${igblast_files}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! nb_seq_per_clone in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID nb_seq_per_clone PARAMETER IN ig_clustering.config FILE:\n${nb_seq_per_clone}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! nb_seq_per_clone =~  /^[0-9]*$/){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID nb_seq_per_clone PARAMETER IN ig_clustering.config FILE:\n${nb_seq_per_clone}\nMUST BE A POSITIVE INTEGER VALUE\n\n========\n\n"
    }
    if( ! clone_distance in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_distance PARAMETER IN ig_clustering.config FILE:\n${clone_distance}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! clone_distance =~  /^[01]{1}\.*[0-9]*$/){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_distance PARAMETER IN ig_clustering.config FILE:\n${clone_distance}\nMUST BE A POSITIVE PROPORTION VALUE\n\n========\n\n"
    }
    if( ! tree_kind in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_kind PARAMETER IN ig_clustering.config FILE:\n${tree_kind}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! tree_leaf_color in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_leaf_color PARAMETER IN ig_clustering.config FILE:\n${tree_leaf_color}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! tree_duplicate_seq in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_duplicate_seq PARAMETER IN ig_clustering.config FILE:\n${tree_duplicate_seq}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (tree_duplicate_seq == "TRUE" || tree_duplicate_seq == "FALSE")){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_duplicate_seq PARAMETER IN ig_clustering.config FILE:\n${tree_duplicate_seq}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if( ! tree_leaf_shape in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_leaf_shape PARAMETER IN ig_clustering.config FILE:\n${tree_leaf_shape}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! tree_leaf_shape =~  /^[0-9]*$/){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_leaf_shape PARAMETER IN ig_clustering.config FILE:\n${tree_leaf_shape}\nMUST BE A POSITIVE INTEGER VALUE\n\n========\n\n"
    }
    if( ! tree_leaf_size in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_leaf_size PARAMETER IN ig_clustering.config FILE:\n${tree_leaf_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! tree_leaf_size =~  /^[0-9]*$/){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_leaf_size PARAMETER IN ig_clustering.config FILE:\n${tree_leaf_size}\nMUST BE A POSITIVE INTEGER VALUE\n\n========\n\n"
    }
    if( ! tree_label_size in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_label_size PARAMETER IN ig_clustering.config FILE:\n${tree_label_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! tree_label_size =~  /^[0-9]*$/){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_label_size PARAMETER IN ig_clustering.config FILE:\n${tree_label_size}\nMUST BE A POSITIVE INTEGER VALUE\n\n========\n\n"
    }
    if( ! tree_label_hjust in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_label_hjust PARAMETER IN ig_clustering.config FILE:\n${tree_label_hjust}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! tree_label_hjust =~  /^\-{0,1}[0-9]+\.*[0-9]*$/){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_label_hjust PARAMETER IN ig_clustering.config FILE:\n${tree_label_hjust}\nMUST BE A NUMERIC VALUE\n\n========\n\n"
    }
    if( ! tree_label_rigth in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_label_rigth PARAMETER IN ig_clustering.config FILE:\n${tree_label_rigth}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (tree_label_rigth == "TRUE" || tree_label_rigth == "FALSE")){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_label_rigth PARAMETER IN ig_clustering.config FILE:\n${tree_label_rigth}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if( ! tree_label_outside in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_label_outside PARAMETER IN ig_clustering.config FILE:\n${tree_label_outside}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (tree_label_outside == "TRUE" || tree_label_outside == "FALSE")){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_label_outside PARAMETER IN ig_clustering.config FILE:\n${tree_label_outside}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if( ! tree_right_margin in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_right_margin PARAMETER IN ig_clustering.config FILE:\n${tree_right_margin}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! tree_right_margin =~  /^[0-9]+\.*[0-9]*$/){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_right_margin PARAMETER IN ig_clustering.config FILE:\n${tree_right_margin}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! tree_legend in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_legend PARAMETER IN ig_clustering.config FILE:\n${tree_legend}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (tree_legend == "TRUE" || tree_legend == "FALSE")){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_legend PARAMETER IN ig_clustering.config FILE:\n${tree_legend}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if( ! tree_meta_path in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_meta_path PARAMETER IN ig_clustering.config FILE:\n${tree_meta_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if(tree_meta_path != "NULL"){
        if( ! file(tree_meta_path).exists()){
            error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_meta_path PARAMETER IN ig_clustering.config FILE. IF DOES NOT EXIST): ${tree_meta_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
        }
    }
    if( ! tree_meta_path_names in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID tree_meta_path_names PARAMETER IN ig_clustering.config FILE:\n${tree_meta_path_names}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! donut_hole_size in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_hole_size PARAMETER IN ig_clustering.config FILE:\n${donut_hole_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! donut_hole_size =~  /^[0-9]+\.*[0-9]*$/){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_hole_size PARAMETER IN ig_clustering.config FILE:\n${donut_hole_size}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! donut_colors in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_colors PARAMETER IN ig_clustering.config FILE:\n${donut_colors}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! cute_path in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID cute_path PARAMETER IN ig_clustering.config FILE:\n${cute_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! file(cute_path).exists()){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID cute_path PARAMETER IN ig_clustering.config FILE (DOES NOT EXIST): ${cute_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }


    // below : those variable are already used in the config file. Thus, to late to check them. And not possible to check inside the config file
    // out_ini
    print("\n\nRESULT DIRECTORY: ${out_path}")
    print("\n\nWARNING: PARAMETERS ALREADY INTERPRETED IN THE .config FILE:")
    print("    system_exec: ${system_exec}")
    print("    out_path: ${out_path_ini}")
    if("${system_exec}" != "local"){
        print("    queue: ${queue}")
        print("    qos: ${qos}")
        print("    add_options: ${add_options}")
    }
    print("\n\n")

    //////// end Checks


    //////// Variable modification


    //////// end Variable modification


    //////// Channels

    fs_ch = Channel.fromPath("${sample_path}/*.*", checkIfExists: false) // in channel because many files

    //////// end Channels


    //////// files import

    meta_file = file(tree_meta_path) // in variable because a single file. If "NULL", will create a file named NULL.
    cute_file=file(cute_path) // in variable because a single file

    //////// end files import


    //////// Main

    workflowParam(
        modules
    )

    repertoire_names(
        igblast_database_path, 
        igblast_files
    )

    igblast(
        fs_ch, 
        igblast_database_path, 
        igblast_files, 
        igblast_organism, 
        igblast_loci
    )

    igblast.out.tsv_ch1.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\n0 ANNOTATION SUCCEEDED BY THE igblast PROCESS\n\nCHECK THAT THE igblast_organism, igblast_loci AND igblast_files ARE CORRECTLY SET IN THE ig_clustering.config FILE\n\n========\n\n"}}
    tsv_ch2 = igblast.out.tsv_ch1.collectFile(name: "all_igblast_seq.tsv", skip: 1, keepHeader: true)
    igblast.out.log_ch.collectFile(name: "igblast_report.log").subscribe{it -> it.copyTo("${out_path}/reports")}


    parseDb_filtering(
        tsv_ch2
    )

    parseDb_filtering.out.unproductive_ch.count().subscribe { n -> if ( n == 0 ){print "\n\nWARNING: EMPTY OUTPUT FOLLOWING THE parseDb_filtering PROCESS -> NO unproductive_seq.tsv FILE RETURNED\n\n"}else{it -> it.copyTo("${out_path}/unproductive_seq.tsv")}} // see https://www.nextflow.io/docs/latest/script.html?highlight=copyto#copy-files


    clone_assignment(
        parseDb_filtering.out.select_ch,
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
    mutation_load_ch2 = mutation_load.out.mutation_load_ch.collectFile(name: "productive_seq.tsv", skip: 1, keepHeader: true)
    mutation_load_ch2.subscribe{it -> it.copyTo("${out_path}")}

    mutation_load.out.mutation_load_log_ch.collectFile(name: "mutation_load.log").subscribe{it -> it.copyTo("${out_path}/reports")} // 
    tuple_mutation_load_ch2 = new Tuple("all", mutation_load_ch2)

    get_tree(
        mutation_load.out.mutation_load_ch,
        nb_seq_per_clone,
        tree_duplicate_seq
    )


    get_tree.out.rdata_tree_ch.count().subscribe { n -> if ( n == 0 ){print("\n\nWARNING: EMPTY OUTPUT FOLLOWING THE get_tree PROCESS -> NO TREE RETURNED\n\n")}}
    rdata_tree_ch2 = get_tree.out.rdata_tree_ch.collect()
    //rdata_tree_ch2.view()

    get_tree.out.no_tree_ch.count().subscribe { n -> if ( n == 0 ){print("\n\nWARNING: ALL SEQUENCES IN TREES FOLLOWING THE get_tree PROCESS -> NO tree_dismissed_seq.tsv FILE RETURNED\n\n")}}
    no_tree_ch2 = get_tree.out.no_tree_ch.collectFile(name: "tree_dismissed_seq.tsv", skip: 1, keepHeader: true)
    no_tree_ch2.subscribe{it -> it.copyTo("${out_path}")}

    get_tree.out.tree_ch.count().subscribe { n -> if ( n == 0 ){print("\n\nWARNING: NO SEQUENCES IN TREES FOLLOWING THE get_tree PROCESS -> NO tree_seq.tsv FILE RETURNED\n\n")}}
    tree_ch2 = get_tree.out.tree_ch.collectFile(name: "tree_seq.tsv", skip: 1, keepHeader: true)
    tree_ch2.subscribe{it -> it.copyTo("${out_path}")}

    get_tree.out.no_cloneID_ch.count().subscribe { n -> if ( n == 0 ){print("\n\nWARNING: ALL SEQUENCES IN CLONAL GROUP FOLLOWING THE get_tree PROCESS -> NO tree_dismissed_clone_id.tsv FILE RETURNED\n\n")}}
    no_cloneID_ch2 = get_tree.out.no_cloneID_ch.collectFile(name: "tree_dismissed_clone_id.tsv")
    no_cloneID_ch2.subscribe{it -> it.copyTo("${out_path}")}

    get_tree.out.cloneID_ch.count().subscribe { n -> if ( n == 0 ){print("\n\nWARNING: NO CLONAL GROUP FOLLOWING THE get_tree PROCESS -> NO tree_clone_id.tsv FILE RETURNED\n\n")}}
    cloneID_ch2 = get_tree.out.cloneID_ch.collectFile(name: "tree_clone_id.tsv")
    cloneID_ch2.subscribe{it -> it.copyTo("${out_path}")}

    get_tree.out.get_tree_log_ch.collectFile(name: "get_tree.log").subscribe{it -> it.copyTo("${out_path}/reports")} // 

    tree_vizu(
        rdata_tree_ch2,
        tree_kind,
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
        meta_file,
        tree_meta_path_names,
        cute_file
    )

    tree_vizu.out.tree_seq_not_displayed_ch.count().subscribe { n -> if ( n == 0 ){print("\n\nWARNING: -> NO tree_seq_not_displayed.tsv FILE RETURNED\n\n")}}
    tree_seq_not_displayed_ch2 = tree_vizu.out.tree_seq_not_displayed_ch.collectFile(name: "tree_seq_not_displayed.tsv", skip: 1, keepHeader: true)
    tree_seq_not_displayed_ch2.subscribe{it -> it.copyTo("${out_path}")}


    tempo1_ch = Channel.of("all", "tree") // 1 channel with 2 values (not list)
    tempo2_ch = mutation_load_ch2.mix(tree_ch2) // 1 channel with 2 paths (flatten() -> not list)
    tempo3_ch = tempo1_ch.merge(tempo2_ch) // 2 lists


    donut(
        tempo3_ch, 
        donut_hole_size, 
        donut_colors,
        cute_file
    )

    donut.out.donut_pdf_ch.count().subscribe { n -> if ( n == 0 ){print("\n\nWARNING: EMPTY OUTPUT FOLLOWING THE donut PROCESS -> NO DONUT RETURNED\n\n")}}
    donut_pdf_ch2 = donut.out.donut_pdf_ch.collect()

    donut.out.donut_tsv_ch.count().subscribe { n -> if ( n == 0 ){print("\n\nWARNING: -> NO tree_seq_not_displayed.tsv FILE RETURNED\n\n")}}
    donut_tsv_ch2 = donut.out.donut_tsv_ch.collectFile(name: "donut_stats.tsv", skip: 1, keepHeader: true)
    donut_tsv_ch2.subscribe{it -> it.copyTo("${out_path}")}

    donut_assembly(
        donut_pdf_ch2
    )

    backup(
        config_file, 
        log_file
    )
}

    //////// end Main


//////// end Processes
