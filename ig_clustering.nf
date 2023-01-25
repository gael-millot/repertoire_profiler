nextflow.enable.dsl=1
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


//////// Arguments of nextflow run

params.modules = ""

//////// end Arguments of nextflow run


//////// Variables

config_file = file("${projectDir}/ig_clustering.config")
log_file = file("${launchDir}/.nextflow.log")
modules = params.modules // remove the dot -> can be used in bash scripts

fs = file(sample_path)

//////// end Variables


//////// Channels


//////// end Channels


//////// Checks

// tbi = file("${sample_path}.tbi") does not need .tbi 

def file_exists1 = fs.exists()
if( ! file_exists1){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN ig_clustering.config FILE (DOES NOT EXIST): ${sample_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
}else{
    fs_ch = Channel.fromPath("${sample_path}/*.*", checkIfExists: false)
}
if( ! igblast_database_path in String ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_database_path PARAMETER IN ig_clustering.config FILE:\n${igblast_database_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! clone_distance in String ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_distance PARAMETER IN ig_clustering.config FILE:\n${clone_distance}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}

// below : those variable are already used in the config file. Thus, to late to check them. And not possible to check inside the config file
// system_exec
// out_ini
print("\n\nRESULT DIRECTORY: ${out_path}")
print("\n\nWARNING: PARAMETERS ALREADY INTERPRETED IN THE .config FILE:")
print("    system_exec: ${system_exec}")
print("    out_path: ${out_path_ini}")
print("    queue: ${queue}")
print("    qos: ${qos}")
print("    add_options: ${add_options}")
print("\n\n")


//////// end Checks


//////// Variable modification


//////// end Variable modification


//////// Processes


process WorkflowVersion { // create a file with the workflow version in out_path
    label 'bash' // see the withLabel: bash in the nextflow config file 
    publishDir "${out_path}/reports", mode: 'copy', overwrite: false
    cache 'false'

    output:
    file "Run_info.txt"

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
    file fs from fs_ch // parallelization expected (for each fasta file)
    val igblast_database_path
    val igblast_files
    val igblast_organism
    val igblast_loci

    output:
    file "*.tsv" optional true into tsv_ch1
    file "*.log" optional true into log_ch

    script:
    """
    #!/bin/bash -ue
    # variables
    REPO_PATH="/usr/local/share/${igblast_database_path}" # path where the imgt_human_IGHV.fasta, imgt_human_IGHD.fasta and imgt_human_IGHJ.fasta files are in the docker container
    VDJ_FILES=\$(awk -v var1="${igblast_files}" -v var2="\${REPO_PATH}" 'BEGIN{ORS=" " ; split(var1, array1, " ") ; for (key in array1) {print var2"/"array1[key]}}')
    FILENAME=\$(basename -- "${fs}") # recover a file name without path
    FILE=\${FILENAME%.*} # file name without extension
    FILE_EXTENSION="\${FILENAME##*.}" #  ## means "delete the longest regex starting at the beginning of the tested string". If nothing, delete nothing. Thus ##*. means delete the longest string finishing by a dot. Use # instead of ## for "delete the shortest regex starting at the beginning of the tested string"
    # end variables
    # checks
    if [[ ! "\${FILE_EXTENSION}" =~ fasta|fa|fna|txt|seq ]] ; then
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nINVALID FILE EXTENSION IN THE sample_path PARAMETER OF THE ig_clustering.config FILE:\\n${fs}\\n\${FILENAME}\\nMUST BE fasta|fs|txt|seq\\n\\n========\\n\\n"
        exit 1
    fi
    sed 's/\\r\$//' ${fs} > tempo_file.fasta #remove carriage returns
    awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' tempo_file.fasta > \${FILE}.fa # remove \\n in the middle of the sequence # \${FILENAME}.fa is a trick to do not use ${fs} and modify the initial file due to the link in the nextflow work folder
    TEMPO=\$(wc -l \${FILE}.fa | cut -f1 -d' ')
    if read -n 1 char <"\${FILE}.fa"; [[ \$char != ">" || \$TEMPO != 2 ]]; then
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nINVALID FASTA FILE IN THE sample_path OF THE ig_clustering.config FILE:\\n${fs}\\nMUST BE A FASTA FILE (\'>\' AS FIRST CHARATER) MADE OF A SINGLE SEQUENCE\\n\\n========\\n\\n"
        exit 1
    fi
    # Alignment <-> annotate sequence using VDJ info
    AssignGenes.py igblast -s \${FILE}.fa -b /usr/local/share/igblast --organism ${igblast_organism} --loci ${igblast_loci} --format blast
    # convert to tsv
    MakeDb.py igblast -i ./\${FILE}_igblast.fmt7 -s ./\${FILE}.fa -r \${VDJ_FILES} --extended
    # printing if no tsv file made
    if [[ ! -f ./\${FILE}_igblast_db-pass.tsv ]] ; then
        echo -e "MakeDb.py igblast FAIL FOR \${FILENAME}" |& tee -a igblast_report.log
    fi
    """
    // write ${} between "" to make a single argument when the variable is made of several values separated by a space. Otherwise, several arguments will be considered
}

tsv_ch1.collectFile(name: "merge.tsv", skip: 1, keepHeader: true).into{tsv_ch2 ; tsv_ch3 ; tsv_ch4 ; tsv_ch5} // concatenate all the cov_report.txt files in channel cov_report_ch into a single file published into ${out_path}/reports
//tsv_ch5.toList().size().view()
if(tsv_ch2.toList().size() == 0){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION: MakeDb.py igblast FAIL FOR ALL THE FILES DURING THE igblast PROCESS\nEMPTY FILE GENERATED\n\n\n\n========\n\n"
}else{
    //print("COUCOU_1")
    tsv_ch3.subscribe{it -> it.copyTo("${out_path}")}
}
//print("COUCOU_3")
log_ch.collectFile(name: "igblast_report.log").subscribe{it -> it.copyTo("${out_path}/reports")} // concatenate all the cov_report.txt files in channel cov_report_ch into a single file published into ${out_path}/reports



process ParseDb_filtering {
    label 'immcantation' // see the withLabel: bash in the nextflow config file 
    publishDir path: "${out_path}", mode: 'copy', pattern: "{*_productive-F.tsv}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{ParseDb_filtering.log}", overwrite: false
    cache 'true'

    input:
    file tsv from tsv_ch4

    output:
    file "*_parse-select.tsv" into select_ch
    file "*_productive-F.tsv" optional true
    file "*.log"

    script:
    """
    #!/bin/bash -ue
    ParseDb.py select -d ${tsv} -f productive -u T |& tee -a ParseDb_filtering.log
    ParseDb.py split -d ${tsv} -f productive |& tee -a ParseDb_filtering.log
    """
}


process clone_assignment {
    label 'immcantation' // see the withLabel: bash in the nextflow config file 
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{*.log}", overwrite: false
    cache 'true'

    input:
    file tsv from select_ch
    val clone_distance

    output:
    file "*_clone-pass.tsv" into clone_ch
    file "*.log"

    script:
    """
    #!/bin/bash -ue
    DefineClones.py -d ${tsv} --act set --model ham --norm len --dist ${clone_distance} |& tee -a clone_assignment.log
    """
}

process closest_germline {
    label 'immcantation' // see the withLabel: bash in the nextflow config file 
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{*.log}", overwrite: false
    cache 'true'

    input:
    file tsv from clone_ch
    val igblast_database_path
    val igblast_files

    output:
    file "*.tsv" into closest_ch
    file "*.log"

    script:
    """
    #!/bin/bash -ue
    # variables
    REPO_PATH="/usr/local/share/${igblast_database_path}" # path where the imgt_human_IGHV.fasta, imgt_human_IGHD.fasta and imgt_human_IGHJ.fasta files are in the docker container
    VDJ_FILES=\$(awk -v var1="${igblast_files}" -v var2="\${REPO_PATH}" 'BEGIN{ORS=" " ; split(var1, array1, " ") ; for (key in array1) {print var2"/"array1[key]}}')
    # end variables
    CreateGermlines.py -d ${tsv} -g dmask --cloned -r \${VDJ_FILES} |& tee -a closest_germline.log
    """
}

process igphyml {
    label 'immcantation' // see the withLabel: bash in the nextflow config file 
    publishDir path: "${out_path}", mode: 'copy', pattern: "{*.tab}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{*.log}", overwrite: false
    cache 'true'

    input:
    file tsv from closest_ch

    output:
    file "phy_igphyml-pass.tab" into tree_ch
    file "igphyml.log"

    script:
    """
    #!/bin/bash -ue
    BuildTrees.py -d ${tsv} --outname phy --log igphyml.log --collapse --sample -1 --igphyml --clean all --nproc 10
    """
}

process tree_vizu {
    label 'r_ext' // see the withLabel: bash in the nextflow config file 
    publishDir path: "${out_path}", mode: 'copy', overwrite: false
    cache 'true'

    input:
    file tab from tree_ch

    output:
    file "tree.pdf"
    file "HLP10_tree_parameters.tsv"

    script:
    """
    #!/usr/bin/env Rscript
    library(alakazam)
    library(igraph)
    db <- alakazam::readIgphyml("./${tab}")
    # Plot largest lineage tree
    pdf(file = "./tree.pdf")
    plot(db\$trees[[1]], layout=layout_as_tree)
    graphics.off()
    # Show HLP10 parameters
    write.table(x = t(db\$param[1,]), file = "./HLP10_tree_parameters.tsv")
    """
}




process Backup {
    label 'bash' // see the withLabel: bash in the nextflow config file 
    publishDir "${out_path}/reports", mode: 'copy', overwrite: false // since I am in mode copy, all the output files will be copied into the publishDir. See \\wsl$\Ubuntu-20.04\home\gael\work\aa\a0e9a739acae026fb205bc3fc21f9b
    cache 'false'

    input:
    file config_file
    file log_file

    output:
    file "${config_file}" // warning message if we use file config_file
    file "${log_file}" // warning message if we use file log_file
    file "Log_info.txt"

    script:
    """
    echo -e "full .nextflow.log is in: ${launchDir}\nThe one in the result folder is not complete (miss the end)" > Log_info.txt
    """
}


//////// end Processes
