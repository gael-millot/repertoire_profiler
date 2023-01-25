nextflow.enable.dsl=2
/*
#########################################################################
##                                                                     ##
##     ig_clustering.nf                                                ##
##                                                                     ##
##     Gael A. Millot                                                  ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Computational Biology Department                                ##
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
        echo "loaded modules (according to specification by the user thanks to the --modules argument of main.nf)": ${modules} >> Run_info.txt
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


process tests {
    // cannot be test outside of process because the files are present in the docker
    label 'immcantation' // see the withLabel: bash in the nextflow config file 
    cache 'true'

    input:
    val igblast_database_path
    val igblast_files

    script:
    """
    #!/bin/bash -ue
    REPO_PATH="/usr/local/share/${igblast_database_path}" # path where the imgt_human_IGHV.fasta, imgt_human_IGHD.fasta and imgt_human_IGHJ.fasta files are in the docker container
    VDJ_FILES=\$(awk -v var1="${igblast_files}" -v var2="\${REPO_PATH}" 'BEGIN{ORS=" " ; split(var1, array1, " ") ; for (key in array1) {print var2"/"array1[key]}}')
    for i1 in \$VDJ_FILES ; do
        if [[ ! -e \${i1} ]] ; then
            echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nFILE DOES NOT EXISTS:\\n\${i1}\\n\\nINDICATED PATH:\\n\${REPO_PATH}\\n\\nCONTAINS:\\n"
            ls -la -w 1 \${REPO_PATH}
            echo -e "\\n\\n========\\n\\n"
            exit 1
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
    val igblast_aa

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
    echo 
    if [[ ${igblast_aa} == "false" ]] ; then
        AssignGenes.py igblast -s \${FILE}.fa -b /usr/local/share/igblast --organism ${igblast_organism} --loci ${igblast_loci} --format blast
    else
        # WARNING: does not work if the fasta file contains \\r (CR carriage return, instead or in addition of \\n, LF line feed) but ok here because removed above
        awk -v var1=\${FILENAME} '{lineKind=(NR-1)%2;}lineKind==0{record=\$0 ; next}lineKind==1{if(\$0~/^[NATGC]*\$/){print "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nFASTA FILE\\n"var1"\\nMUST BE AN AMINO ACID SEQUENCE IF THE igblast_aa PARAMETER IS SET TO true\\nHERE IT SEEMS ONLY MADE OF NUCLEOTIDES:\\n"\$0"\\n\\n========\\n\\n" ; exit 1}}' \${FILE}.fa
        AssignGenes.py igblast-aa -s \${FILE}.fa -b /usr/local/share/igblast --organism ${igblast_organism} --loci ${igblast_loci}
    fi

    # convert to tsv
    # Also convert data from the web interface IMGT/HighV-QUEST
    if [[ ${igblast_aa} == "false" ]] ; then
        MakeDb.py igblast -i ./\${FILE}_igblast.fmt7 -s ./\${FILE}.fa -r \${VDJ_FILES} --extended
    else
        MakeDb.py igblast-aa -i ./\${FILE}_igblast.fmt7 -s ./\${FILE}.fa -r \${VDJ_FILES} --extended
    fi

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
    path tsv_ch2

    output:
    path "*_parse-select.tsv", emit: select_ch
    path "*_productive-F.tsv", optional: true
    path "*.log"

    script:
    """
    #!/bin/bash -ue
    if [[ ! -s ${tsv_ch2} ]]; then
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nEMPTY FILE GENERATED BY THE igblast PROCESS\\nCHECK THE igblast_report.log FILE IN THE report FOLDER\\n\\n========\\n\\n"
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
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nEMPTY FILE GENERATED BY THE parseDb_filtering PROCESS\\nCHECK THE ParseDb_filtering.log AND *_productive-F.tsv FILES IN THE report FOLDER\\n\\n========\\n\\n"
        exit 1
    else
        DefineClones.py -d ${select_ch} --act set --model ham --norm len --dist ${clone_distance} |& tee -a clone_assignment.log
    fi
    """
}

process closest_germline {
    label 'immcantation' // see the withLabel: bash in the nextflow config file 
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{*.log}", overwrite: false
    cache 'true'

    input:
    path clone_ch
    val igblast_database_path
    val igblast_files

    output:
    path "*.tsv", emit: closest_ch
    path "*.log"

    script:
    """
    #!/bin/bash -ue
    # variables
    REPO_PATH="/usr/local/share/${igblast_database_path}" # path where the imgt_human_IGHV.fasta, imgt_human_IGHD.fasta and imgt_human_IGHJ.fasta files are in the docker container
    VDJ_FILES=\$(awk -v var1="${igblast_files}" -v var2="\${REPO_PATH}" 'BEGIN{ORS=" " ; split(var1, array1, " ") ; for (key in array1) {print var2"/"array1[key]}}')
    # end variables
    CreateGermlines.py -d ${clone_ch} -g dmask --cloned -r \${VDJ_FILES} |& tee -a closest_germline.log
    """
}

process igphyml {
    label 'immcantation' // see the withLabel: bash in the nextflow config file 
    publishDir path: "${out_path}", mode: 'copy', pattern: "{*.tab}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{*.log}", overwrite: false
    cache 'true'

    input:
    path closest_ch

    output:
    path "phy_igphyml-pass.tab", emit: tree_ch
    path "igphyml.log"

    script:
    """
    #!/bin/bash -ue
    BuildTrees.py -d ${closest_ch} --outname phy --log igphyml.log --collapse --sample -1 --igphyml --clean all --nproc 10
    """
}

process tree_vizu {
    label 'r_ext' // see the withLabel: bash in the nextflow config file 
    publishDir path: "${out_path}", mode: 'copy', overwrite: false
    cache 'true'

    input:
    path tree_ch

    output:
    path "*_tree.pdf"
    path "HLP10_tree_parameters.tsv"

    script:
    """
    #!/usr/bin/env Rscript
    suppressMessages(library(alakazam))
    suppressMessages(library(igraph))
    suppressMessages(library(ape))

    # Plot basic ugly lineage tree
    db <- alakazam::readIgphyml("./${tree_ch}")
    pdf(file = "./ugly_tree.pdf")
    igraph::plot.igraph(db\$trees[[1]], layout = layout_as_tree, edge.arrow.mode = 0, vertex.frame.color = "black", vertex.size = 4,  edge.label = NA, vertex.label.dist = 1, vertex.label.cex = 0.5)
    # see library/igraph/html/plot.common.html
    # add vertex. or edge. before the described arguments
    # see library/igraph/html/plot.igraph.html

    pdf(file = "./nice_tree.pdf")
    db2 <- alakazam::readIgphyml("./${tree_ch}", format = "phylo")
    ape::plot.phylo(ape::ladderize(db2\$trees[[1]]), cex = 0.7, no.margin = TRUE)

    graphics.off()

    # Show HLP10 parameters
    write.table(x = t(db\$param[1,]), file = "./HLP10_tree_parameters.tsv")
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
    echo -e "full .nextflow.log is in: ${launchDir}\nThe one in the result folder is not complete (miss the end)" > Log_info.txt
    """
}








workflow {

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


    //////// Channels


    //////// end Channels


    //////// Checks

    // tbi = file("${sample_path}.tbi") does not need .tbi 

    def file_exists1 = file(sample_path).exists()
    if( ! file_exists1){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN ig_clustering.config FILE (DOES NOT EXIST): ${sample_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }else{
        fs_ch = Channel.fromPath("${sample_path}/*.*", checkIfExists: false)
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
    if( ! igblast_aa in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_aa PARAMETER IN ig_clustering.config FILE:\n${igblast_aa}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (igblast_aa == "true" || igblast_aa == "false")){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_aa PARAMETER IN ig_clustering.config FILE:\n${igblast_aa}\nMUST BE EITHER \"true\" OR \"false\"\n\n========\n\n"
    }
    if( ! clone_distance in String ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_distance PARAMETER IN ig_clustering.config FILE:\n${clone_distance}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! clone_distance =~  /^[01]{1}\.*[0-9]*$/){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_distance PARAMETER IN ig_clustering.config FILE:\n${clone_distance}\nMUST BE A POSITIVE PROPORTION VALUE\n\n========\n\n"
    }

    // below : those variable are already used in the config file. Thus, to late to check them. And not possible to check inside the config file
    // system_exec
    // out_ini
    print("\n\nRESULT DIRECTORY: ${out_path}")
    if(igblast_aa == "true"){
        print("\n\nWARNING: ANNOTATIONS OF SEQUENCES AT THE AA LEVEL, ACCORDING TO THE igblast_aa PARAMETER")
    }
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


    //////// Main

    workflowParam(
        modules
    )

    tests(
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

    tsv_ch2 = igblast.out.tsv_ch1.collectFile(name: "merge.tsv", skip: 1, keepHeader: true) // concatenate all the cov_report.txt files in channel cov_report_ch into a single file published into ${out_path}/reports
    // or igblast.out[0].collectFile() if output of igblast process is not named using emit.
    //the following test does not work because tested at the beginning of the nexflow run, not after the igblast process, as shown with the print("coucou")
    //if(tsv_ch2.toList().size() == 0){
        //print("coucou_1")
        //error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION: MakeDb.py igblast FAIL FOR ALL THE FILES DURING THE igblast PROCESS\nEMPTY FILE GENERATED\n\n\n\n========\n\n"
    //}else{
        //print("coucou_2")
        //
    //}
    tsv_ch2.subscribe{it -> it.copyTo("${out_path}")}
    igblast.out.log_ch.collectFile(name: "igblast_report.log").subscribe{it -> it.copyTo("${out_path}/reports")} // concatenate all the cov_report.txt files in channel cov_report_ch into a single file published into ${out_path}/reports


    parseDb_filtering(
        tsv_ch2
    )

    clone_assignment(
        parseDb_filtering.out.select_ch,
        clone_distance
    )

    closest_germline(
        clone_assignment.out.clone_ch, 
        igblast_database_path, 
        igblast_files
    )

    igphyml(
        closest_germline.out.closest_ch
    )

    tree_vizu(
        igphyml.out.tree_ch
    )

    backup(
        config_file, 
        log_file
    )
}

    //////// end Main

//////// end Processes
