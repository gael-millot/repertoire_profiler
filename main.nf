nextflow.enable.dsl=2
/*
#########################################################################
##                                                                     ##
##     main.nf of repertoire profiler                                  ##
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


// Compares input fasta files with reference files to determine which cassettes are present
// The input fasta files can be nucleotidic or amino acid sequences (igblast_aa in nextflow.config)
// Output :
// - tsv_ch1 : tsv files with the sequences aligned with the reference files (annotated by igblast)
//             file used in the next process
// - igblast_aligned_seq.tsv: names of the sequences aligned with the reference files
// - igblast_unaligned_seq.tsv: names of the sequences not aligned with the reference files
process igblast {
    label 'immcantation'
    //publishDir path: "${out_path}", mode: 'copy', overwrite: false
    cache 'true'

    input:
    path fs_ch // parallelization expected (for each fasta file)
    val igblast_variable_ref_files
    val igblast_organism
    val igblast_loci
    val igblast_aa

    output:
    path "*_igblast_db-pass.tsv", emit: tsv_ch1, optional: true
    path "igblast_unaligned_seq_name.tsv", emit: unaligned_seq_ch
    path "igblast_aligned_seq_name.tsv", emit: aligned_seq_ch
    path "*.log", emit: log_ch, optional: true

    script:
    """
    #!/bin/bash -ue

    # variables
    REPO_PATH="/usr/local/share/germlines/imgt/${igblast_organism}/vdj" # path where the imgt_human_IGHV.fasta, imgt_human_IGHD.fasta and imgt_human_IGHJ.fasta files are in the docker container
    VDJ_FILES=\$(awk -v var1="${igblast_variable_ref_files}" -v var2="\${REPO_PATH}" 'BEGIN{ORS=" " ; split(var1, array1, " ") ; for (key in array1) {print var2"/"array1[key]}}')
    FILENAME=\$(basename -- "${fs_ch}") # recover a file name without path
    FILE=\${FILENAME%.*} # file name without extension
    MODIF_FILE="\${FILE// /_}" # spaces in the name file replaced by _
    FILE="\${MODIF_FILE}"
    FILE_EXTENSION="\${FILENAME##*.}" #  ## means "delete the longest regex starting at the beginning of the tested string". If nothing, delete nothing. Thus ##*. means delete the longest string finishing by a dot. Use # instead of ## for "delete the shortest regex starting at the beginning of the tested string"
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a igblast_report.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a igblast_report.log
    # end variables

    # checks
    if [[ ! "\${FILE_EXTENSION}" =~ fasta|fa|fas|fna|txt|seq ]] ; then
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nINVALID FILE EXTENSION IN THE sample_path PARAMETER OF THE repertoire_profiler.config FILE:\\n${fs_ch}\\n\${FILENAME}\\nMUST BE fasta|fa|fas|fna|txt|seq|faa\\n\\n========\\n\\n"
        exit 1
    fi
    sed 's/\\r\$//' ${fs_ch} > tempo_file.fasta # remove carriage returns
    awk 'BEGIN{ORS=""}{if(\$0~/^>.*/){gsub(/[ \t.-]+/, "_", \$0) ; if(NR>1){print "\\n"} ; print \$0"\\n"} else {print \$0 ; next}}END{print "\\n"}' tempo_file.fasta > \${FILE}.fa # remove \\n in the middle of the sequence # gsub(/[ \t]+/, "_", \$0) replace spaces and tabs in the first line by a single underscore# \${FILENAME}.fa is a trick to do not use ${fs_ch} and modify the initial file due to the link in the nextflow work folder
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
        echo "\${FILENAME}" | cat > igblast_unaligned_seq_name.tsv
        echo -n "" | cat > igblast_aligned_seq_name.tsv
    else
        echo "\${FILENAME}" | cat > igblast_aligned_seq_name.tsv
        echo -n "" | cat > igblast_unaligned_seq_name.tsv
    fi
    """
    // write ${} between "" to make a single argument when the variable is made of several values separated by a space. Otherwise, several arguments will be considered
}

// Parse the tsv files generated by the igblast process to select the productive and unproductive sequences
// Output :
// - productive_seq_init.tsv : tsv file with the productive sequences
// - unproductive_seq.tsv : tsv file with the unproductive sequences
process parseDb_filtering {
    label 'immcantation'
    cache 'true'

    input:
    path tsv_ch1 // parallelization expected
    val igblast_aa

    output:
    path "productive_seq_init.tsv", emit: select_ch, optional: true
    path "unproductive_seq.tsv", emit: unselect_ch
    path "ParseDb_filtering.log", emit: parseDb_filtering_log_ch

    script:
    """
    #!/bin/bash -ue
    FILENAME=\$(basename -- "${tsv_ch1}") # recover a file name without path
    FILE=\${FILENAME%.*} # file name without extension
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a ParseDb_filtering.log
    if [[ "${igblast_aa}" == "true" ]]; then # if igblast_aa is true, then the productive column is empty because aa sequences means productive ones
        ParseDb.py select -d ${tsv_ch1} -f v_call j_call -u ".+" --regex --logic any |& tee -a ParseDb_filtering.log #means look inside the -f v_call j_call fields of the input and return any lines that are non empty for at least one field(--logic any) # should be identical to cp ${tsv_ch1} "\${FILE}_parse-select.tsv" |& tee -a ParseDb_filtering.log
        echo -n "" > unproductive_seq.tsv |& tee -a ParseDb_filtering.log
    elif [[ -s ${tsv_ch1} ]]; then # -s means "exists and non empty". Thus, return FALSE is the file does not exists or is empty
        ParseDb.py select -d ${tsv_ch1} -f productive -u T |& tee -a ParseDb_filtering.log
        ParseDb.py split -d ${tsv_ch1} -f productive |& tee -a ParseDb_filtering.log
        if [ -f *_parse-select.tsv ]; then
            cp *_parse-select.tsv productive_seq_init.tsv |& tee -a ParseDb_filtering.log # can be empty file (only header)
        else
            echo -e "\\n\\n========\\n\\nINTERNAL ERROR IN NEXTFLOW EXECUTION\\n\\nNO *_parse-select.tsv FILE GENERATED BY THE igblast PROCESS\\nCHECK THE ParseDb_filtering.log FILE IN THE report FOLDER INSIDE THE OUTPUT FOLDER\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n"
            exit 1
        fi
        if [ -s *_productive-F.tsv ]; then # see above for -s
            cp *_productive-F.tsv unproductive_seq_init.tsv |& tee -a ParseDb_filtering.log
        else
            echo -e "\n\nWARNING: EMPTY unproductive_seq.tsv FILE RETURNED FOLLOWING THE parseDb_filtering PROCESS\n\n" |& tee -a ParseDb_filtering.log
            # echo -n "" | cat > unproductive_seq.tsv
            head -1 productive_seq_init.tsv | cat > unproductive_seq.tsv # only header in the file
        fi
    else
        echo -e "\\n\\n========\\n\\nINTERNAL ERROR IN NEXTFLOW EXECUTION\\n\\nEMPTY FILE GENERATED BY THE igblast PROCESS\\nCHECK THE ParseDb_filtering.log FILE IN THE report FOLDER INSIDE THE OUTPUT FOLDER\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n" # I said empty because the existence of the file has been checked after the igblast process
        exit 1
    fi
    """
}


// Translate the nucleotidic sequences into amino acid sequences
process translation {
    label 'seqkit'
    publishDir path: "${out_path}", mode: 'copy', pattern: "aa/*.fasta", overwrite: false
    publishDir path: "${out_path}", mode: 'copy', pattern: "aligned_seq/*.fasta", overwrite: false
    cache 'true'

    input:
    path select_ch // parallelization expected
    val igblast_aa

    output:
    path "translation.tsv", optional: true, emit: translation_ch // productive file with column sequence_alignment_aa added
    path "aligned_seq/*.*", optional: true
    path "aa.tsv", optional: true, emit: aa_tsv_ch
    path "aa/*.*", optional: true
    path "translation.log", optional: true, emit: translation_log_ch

    script:
    """
    #!/bin/bash -ue
    translation.sh ${select_ch} ${igblast_aa} # |& tee -a translation.log not used because in translation.sh 
    """
}




process distToNearest {
    label 'immcantation'
    //publishDir path: "${out_path}", mode: 'copy', pattern: "{nearest_distance.tsv}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{distToNearest.log}", overwrite: false
    cache 'true'

    input:
    path translation_ch2 // no parallelization
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
         # WEIRD stuf: if db alone is returned, and if distToNearest_ch is used for the clone_assignment process and followings, everything is fine. But if db3 is returned with db3 <- data.frame(db, dist_nearest = db2\$dist_nearest) or db3 <- data.frame(db, caca = db2\$dist_nearest) or data.frame(db, caca = db\$sequence_id) or db3 <- data.frame(db, caca = as.numeric(db2\$dist_nearest)) or db3 <- data.frame(db[1:3], caca = db\$sequence_id, db[4:length(db)]), the get_germ_tree process cannot make trees, while the productive.tsv seem identical at the end, between the use of db or db3, except that the clone_id order is not the same
        db <- read.table("${translation_ch2}", header = TRUE, sep = "\\t")
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
    publishDir path: "${out_path}/files", mode: 'copy', pattern: "{seq_distance.pdf}", overwrite: false
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
    publishDir path: "${out_path}/files", mode: 'copy', pattern: "non_clone_assigned_sequence.tsv", overwrite: false
    cache 'true'

    input:
    path productive_ch // no parallelization
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
    if [[ -s ${productive_ch} ]]; then # see above for -s
        DefineClones.py -d ${productive_ch} --act set --model ${clone_model} --norm ${clone_normalize} --dist ${clone_distance} --fail |& tee -a clone_assignment.log
        if [ -s *_clone-fail.tsv ]; then # see above for -s
            cp *_clone-fail.tsv non_clone_assigned_sequence.tsv |& tee -a clone_assignment.log
        else
            echo -e "\n\nNOTE: EMPTY non_clone_assigned_sequence.tsv FILE RETURNED FOLLOWING THE clone_assignment PROCESS\n\n" |& tee -a clone_assignment.log
            echo -n "" | cat > non_clone_assigned_sequence.tsv
        fi
    else
        echo -e "\\n\\n========\\n\\nINTERNAL ERROR IN NEXTFLOW EXECUTION\\n\\nEMPTY FILE GENERATED BY THE translation PROCESS\\nCHECK THE clone_assignment.log AND *_productive-F.tsv FILES IN THE report FOLDER INSIDE THE OUTPUT FOLDER\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n"
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


process closest_germline {
    label 'immcantation'
    //publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{*.log}", overwrite: false
    //publishDir path: "${out_path}", mode: 'copy', pattern: "{*.tsv}", overwrite: false
    cache 'true'

    input:
    path clone_split_ch // parallelization expected
    val igblast_organism
    val igblast_variable_ref_files

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

    REPO_PATH="/usr/local/share/germlines/imgt/${igblast_organism}/vdj" # path where the imgt_human_IGHV.fasta, imgt_human_IGHD.fasta and imgt_human_IGHJ.fasta files are in the docker container
    VDJ_FILES=\$(awk -v var1="${igblast_variable_ref_files}" -v var2="\${REPO_PATH}" 'BEGIN{ORS=" " ; split(var1, array1, " ") ; for (key in array1) {print var2"/"array1[key]}}')
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
          tempo.cat <- paste0("\\n\\n========\\n\\nINTERNAL ERROR IN THE NEXTFLOW EXECUTION OF THE mutation_load PROCESS\\n THE NUMBER OF ELEMENTS IN args (", length(args),") IS DIFFERENT FROM THE NUMBER OF ELEMENTS IN tempo.arg.names (", length(tempo.arg.names),")\\nargs:", paste0(args, collapse = ","), "\\ntempo.arg.names:", paste0(tempo.arg.names, collapse = ","), "\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n")
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
        # write.table(VDJ_db, file = paste0("./", file_name, "_shm-pass.tsv"), row.names = FALSE, sep = "\\t")

        # Reposition the column starting with "initial_" to the second position
        initial_col <- grep("^initial_", names(VDJ_db), value = TRUE)
        if (length(initial_col) == 0) {
          stop("No column starting with 'initial_' found")
        }
        tempo <- VDJ_db[[initial_col]]
        VDJ_db <- VDJ_db[ , !(names(VDJ_db) %in% initial_col)]
        VDJ_db <- data.frame(VDJ_db[1], tempo, VDJ_db[2:length(VDJ_db)])
        names(VDJ_db)[2] <- initial_col

        # Check if the column "kd" exists and rename it to "KD"
        if ("kd" %in% names(VDJ_db)) {
          names(VDJ_db)[names(VDJ_db) == "kd"] <- "KD"
        }

        write.table(VDJ_db, file = paste0("./", file_name, "_shm-pass.tsv"), row.names = FALSE, col.names = TRUE, sep = "\\t")

    ' "\${FILENAME}" |& tee -a mutation_load.log
    # cp ./tempo_shm-pass.tsv \${FILENAME}_shm-pass.tsv
    # rm tempo_shm-pass.tsv
    """
}

// Rename the sequence names in the data file
// The sequence names are replaced by the values of the column meta_name_replacement (specified in .config) of the metadata file
process seq_name_replacement {
    label 'r_ext'
    cache 'true'

    input:
    path translation_ch // parallelization expected
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
    FILENAME=\$(basename -- ${translation_ch}) # recover a file name without path
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a seq_name_replacement.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a seq_name_replacement.log
    # check first that the data file does not have the second column name starting by "initial_". Otherwise, with create unproper behavior in donut
    if [[ "${meta_file}" == "NULL" ]] ; then
        rm NULL # remove the initial file to avoid to send it into the channel
        echo -n "" > NULL # new hard file that can be sent into the channel
        chmod 777 NULL
        Rscript -e '
            seq <- read.table("./${translation_ch}", sep = "\\t", header = TRUE)
            if(grepl(x = names(seq)[2], pattern = "^initial_")){
                stop(paste0("\\n\\n============\\n\\nERROR IN THE seq_name_replacement PROCESS OF NEXTFLOW\\nIF THE meta_path PARAMETER IS \\"NULL\\", THEN THE SECOND COLUMN OF THE DATA IN THE sample_path PARAMETER CANNOT HAVE THE NAME OF THE SECOND COLUNM STARTING BY \\"initial_\\"\\n\\n============\\n\\n"), call. = FALSE)
            }
        ' |& tee -a seq_name_replacement.log
        IFS='_' read -r -a TEMPO <<< "\${FILENAME}" # string split into array
        cat ${translation_ch} > ./\${TEMPO[0]}_renamed_seq.tsv |& tee -a seq_name_replacement.log
    else
        # if [[ "${meta_file}" != "NULL" && "${meta_name_replacement}" != "NULL" ]] ; then # or [[ "${meta_file}" -ne "NULL" && "${meta_name_replacement}" -ne "NULL" ]], but not !=
        Rscript -e '
            meta <- read.table("./${meta_file}", sep = "\\t", header = TRUE)
            seq <- read.table("./${translation_ch}", sep = "\\t", header = TRUE)
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
            write.table(seq, file = paste0("./", id, "_renamed_seq.tsv"), row.names = FALSE, col.names = TRUE, sep = "\\t")
            # modification of the metadata file for the correct use of ggtree::"%<+%" in germ_tree_vizu.R that uses the column name meta_seq_names for that 
            # meta <- data.frame(meta, initial_label = meta[ , "${meta_seq_names}"])
            # meta[ , "${meta_seq_names}"] <- meta[ , "${meta_name_replacement}"]
            # write.table(meta, file = "./metadata2.tsv", row.names = FALSE, col.names = TRUE, sep = "\\t")
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
    publishDir path: "${out_path}/files", mode: 'copy', pattern: "{productive_seq.tsv}", overwrite: false
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
        # remove allele info
        tempo_v <- strsplit(db3\$v_call, ",")
        sub_v <- sapply(X = tempo_v, FUN = function(x){y <- sub(pattern = "\\\\*.*", replacement = "", x = x) ; paste0(unique(y), collapse = ",")})
        tempo_j <- strsplit(db3\$j_call, ",")
        sub_j <- sapply(X = tempo_j, FUN = function(x){y <- sub(pattern = "\\\\*.*", replacement = "", x = x) ; paste0(unique(y), collapse = ",")})
        tempo_subclass <- strsplit(db3\$c_call, ",")
        subclass <- sapply(X = tempo_subclass, FUN = function(x){y <- sub(pattern = "\\\\*.*", replacement = "", x = x) ; paste0(unique(y), collapse = ",")})
        class <- sapply(X = tempo_subclass, FUN = function(x){y <- sub(pattern = "\\\\*.*", replacement = "", x = x) ; y <- substr(y, 1, 4) ; paste0(unique(y), collapse = ",")})
        db4 <- data.frame(db3, v_gene = sub_v, j_gene = sub_j, isotype_class = class, c_gene = subclass)
        # end remove allele info
        write.table(db4, file = paste0("./productive_seq.tsv"), row.names = FALSE, col.names = TRUE, sep = "\\t")
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
    label 'r_ext'
    publishDir path: "${out_path}/files", mode: 'copy', pattern: "{*_repertoire.pdf}", overwrite: false
    publishDir path: "${out_path}/png", mode: 'copy', pattern: "{*.png}", overwrite: false
    publishDir path: "${out_path}/svg", mode: 'copy', pattern: "{*.svg}", overwrite: false
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
    path "*.png", emit: repertoire_png
    path "rep_*.tsv"
    path "repertoire.log"

    script:
    """
    #!/bin/bash -ue
    repertoire.R \
"${seq_name_replacement_ch2}" \
"${igblast_data_check_ch}" \
"${cute_file}" \
"repertoire.log"
    """
}



process get_germ_tree {
    label 'immcantation_10cpu'
    publishDir path: "${out_path}/RData", mode: 'copy', pattern: "{*_get_germ_tree_cloneID.RData}", overwrite: false
    cache 'true'

    input:
    path mutation_load_ch // parallelization expected
    path meta_file // just to determine if metadata have been provided (TRUE means NULL) meta_file_ch not required here
    path cute_file
    val clone_nb_seq
    val germ_tree_duplicate_seq
    val igphylm_exe_path // warning : here val and not path because we do not want the igphyml file to be imported in the work dir

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
"${clone_nb_seq}" \
"${germ_tree_duplicate_seq}" \
"${igphylm_exe_path}" \
"${cute_file}" \
"get_germ_tree.log"
    """
}


process germ_tree_vizu {
    label 'r_ext'
    publishDir path: "${out_path}/files", mode: 'copy', pattern: "{germ_tree.pdf}", overwrite: false
    publishDir path: "${out_path}/files", mode: 'copy', pattern: "{germ_no_tree.pdf}", overwrite: false
    publishDir path: "${out_path}/png", mode: 'copy', pattern: "{*.png}", overwrite: false
    publishDir path: "${out_path}/svg", mode: 'copy', pattern: "{*.svg}", overwrite: false
    publishDir path: "${out_path}/RData", mode: 'copy', pattern: "{all_trees.RData}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{germ_tree_vizu.log}", overwrite: false
    cache 'true'

    input:
    path rdata_germ_tree_ch2 // no more parallelization
    val germ_tree_kind
    val clone_nb_seq
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
    germ_tree_vizu.R \
"${germ_tree_kind}" \
"${clone_nb_seq}" \
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



process donut {
    label 'r_ext'
    //publishDir path: "${out_path}", mode: 'copy', pattern: "{*.tsv}", overwrite: false
    //publishDir path: "${out_path}", mode: 'copy', pattern: "{*.pdf}", overwrite: false
    publishDir path: "${out_path}/png", mode: 'copy', pattern: "{*.png}", overwrite: false
    publishDir path: "${out_path}/svg", mode: 'copy', pattern: "{*.svg}", overwrite: false
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

    output:
    path "*.tsv", emit: donut_tsv_ch, optional: true
    path "*.pdf", emit: donut_pdf_ch, optional: true
    path "*.png", emit: donuts_png
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
"${kind}_donut.log"
    """
}

process donut_assembly {
    label 'r_ext'
    publishDir path: "${out_path}/files", mode: 'copy', pattern: "{donuts.pdf}", overwrite: false
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
    echo -e "full .nextflow.log is in: ${launchDir}\nThe one in the result folder is not complete (miss the end)" > Log_info.txt
    """
}



// Converts a tsv file containing several amino acid sequences into fasta file
// Input : tsv file (columns sequence_id (CL4184329_VH) and sequence_alignment (QVQLQQSGAXELARPGASVK...))
// Output : same data converted into a fasta file
process Reformat{
    
    input:
    path aatsv
    
    output:
    path "*.fasta", emit : aa_fasta_ch
    
    script:
    """
    cat $aatsv | tail -n+2 | awk '{print ">"\$1; print \$2}' > ${aatsv.baseName}.fasta
    """
    
}

process Align {
    label 'abalign'
    
    input:
    path fasta
    val phylo_tree_heavy
    
    output:
    path "*_align.fasta", emit : aligned_groups_ch
    
    script:
    parms="-al"
    if(phylo_tree_heavy){parms="-ah"}
    """
    /bin/Abalign_V2_Linux_Term/Abalign -i $fasta ${parms} ${fasta.baseName}_align.fasta -sp MU || true
    """
}

// Groups input sequences (amino-acid ones) into same vj combination
// Input : single fasta file containing several sequences & a tsv file containing V|J info to corresponding sequence id
// Output : one fasta file PER group, containing said group sequences
process DefineGroups {
    publishDir path: "${out_path}/phylo", mode: 'copy'
    
    label 'goalign'
    
    input:
    path fasta
    path productivetsv
    
    output:
    path "sequences_*.fasta", emit : groups_ch
    
    script:
    """
    # On prend toutes les combinaisons existantes de V|J
    defineGroups.pl $productivetsv > groups.txt
    # for each defined group
    for vj in \$(cut -f 2 groups.txt | sort -u)
    do
	# we take from the fasta the sequences corresponding to that group
        goalign subset --unaligned -i aa.fasta -f <(grep -e "\$vj" groups.txt | cut -f 1) > sequences_\${vj}.fasta
    done
    """
}

// Extract the number of sequences from the input file
process NbSequences {
    label 'goalign'

    input:
    path fasta

    output:
    tuple stdout, path(fasta), emit : nb_out

    script:
    """
    printf \$(goalign stats nseq --unaligned -i $fasta)
    """
}

process Tree {
    publishDir path: "${out_path}/files", mode: 'copy'

    label 'iqtree'

    input:
    path msa
    path phylo_tree_model_file

    output:
    path "*.treefile", emit : tree_file
    path "*.log"

    script:
    """
    iqtree -nt ${task.cpus} -s ${msa} -m ${phylo_tree_model_file}+I+R6 --seed 123456789
    """
}

process ProcessMeta {
    label 'tabletoitol'

    input:
    path meta

    output:
    path "iTOL*", emit : itol_out

    script:
    """
    table2itol.R -i ${meta_seq_names} ${meta}
    """
}


process ITOL{
    publishDir path: "${out_path}/files", mode: 'copy'

    label 'gotree'

    input:
    path tree
    path meta
    val phylo_tree_itolkey

    output:
    path "*itol_url.txt"

    script:
    """
    gotree upload itol --project gotree_uploads --user-id $phylo_tree_itolkey -i $tree $meta > ${tree.baseName}_itol_url.txt 2>&1
    """



}

// The render function only creates the html file with the rmardown
// Inputs :
//      template_rmd : rmardown template file used to create the html (file path is defined in config.nextflow)
//      nb_input : number of fasta sequences in initial input
//      nb_seq_aligned : number of sequences that igblast could analyse
//      nb_seq_not_aligned : number of sequences that igblast could not alayse
//      donuts_png, repertoire_png : to make sure the report process is not called before processes that produces images it needs
// Outputs :
//      "report.html" : finalized html report for a specific run
//      "print_report.log" : will contain any error or warning messages produced by rmardown::render
process print_report{
    label 'r_ext'

    publishDir path: "${out_path}", mode: 'copy', pattern: "{*.html}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{print_report.log}", overwrite: false
    cache 'false'

    input:
    file template_rmd
    val nb_input
    val nb_seq_aligned
    val nb_seq_not_aligned
    val nb_productive
    val nb_unproductive
    val nb_clone_assigned
    val nb_failed_clone
    path donuts_png
    path repertoire_png

    output:
    file "report.html"
    file "print_report.log"

    script:
    """
    #!/bin/bash -ue
    cp ${template_rmd} report_file.rmd
    cp -r "${out_path}/png" .
    Rscript -e '
        rmarkdown::render(
        input = "report_file.rmd",
        output_file = "report.html",
        # list of the variables waiting to be replaced in the rmd file :
        params =    list(nb_input = ${nb_input},
                        nb_seq_aligned = ${nb_seq_aligned}, 
                        nb_seq_not_aligned = ${nb_seq_not_aligned},
                        nb_productive = ${nb_productive},
                        nb_unproductive = ${nb_unproductive},
                        nb_clone_assigned = ${nb_clone_assigned},
                        nb_failed_clone = ${nb_failed_clone}),
        # output_dir = ".",
        # intermediates_dir = "./",
        # knit_root_dir = "./",
        run_pandoc = TRUE,
        quiet = TRUE,
        clean = TRUE
        )
    ' |& tee -a print_report.log
    """
}

//////// End Processes

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
    config_file = workflow.configFiles[0] // better to use this than config_file = file("${projectDir}/repertoire_profiler.config") because the latter is not good if -c option of nextflow run is used
    log_file = file("${launchDir}/.nextflow.log")

    //////// end Variables


    //////// Checks

    //// check of the bin folder
    if( ! (file("${projectDir}/bin/donut.R").exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nTHE donut.R FILE MUST BE PRESENT IN THE ./bin FOLDER, WHERE THE main.nf file IS PRESENT\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    // AB_model not trested because in parameters
    if( ! (file("${projectDir}/bin/germ_tree_vizu.R").exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nTHE germ_tree_vizu.R FILE MUST BE PRESENT IN THE ./bin FOLDER, WHERE THE main.nf file IS PRESENT\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    if( ! (file("${projectDir}/bin/get_germ_tree.R").exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nTHE get_germ_tree.R FILE MUST BE PRESENT IN THE ./bin FOLDER, WHERE THE main.nf file IS PRESENT\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    if( ! (file("${projectDir}/bin/get_germ_tree.R").exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nTHE get_germ_tree.R FILE MUST BE PRESENT IN THE ./bin FOLDER, WHERE THE main.nf file IS PRESENT\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    if( ! (file("${projectDir}/bin/repertoire.R").exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nTHE repertoire.R FILE MUST BE PRESENT IN THE ./bin FOLDER, WHERE THE main.nf file IS PRESENT\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    if( ! (file("${projectDir}/bin/translation.sh").exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nTHE translation.sh FILE MUST BE PRESENT IN THE ./bin FOLDER, WHERE THE main.nf file IS PRESENT\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    if( ! (file("${projectDir}/bin/defineGroups.pl").exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nTHE defineGroups.pl FILE MUST BE PRESENT IN THE ./bin FOLDER, WHERE THE main.nf file IS PRESENT\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    //// end check of the bin folder
    if( ! (sample_path in String || sample_path in GString) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN repertoire_profiler.config FILE:\n${sample_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (file(sample_path).exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN repertoire_profiler.config FILE (DOES NOT EXIST): ${sample_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    if( ! (igblast_organism in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_organism PARAMETER IN repertoire_profiler.config FILE:\n${igblast_organism}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (igblast_organism =~ /^(mouse|human|rabbit|rat|rhesus_monkey)$/) ){
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
    if( ! (igblast_variable_ref_files in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_variable_ref_files PARAMETER IN repertoire_profiler.config FILE:\n${igblast_variable_ref_files}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    } // path not checked here but below because inside a container
    if( ! (igblast_constant_ref_files in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_constant_ref_files PARAMETER IN repertoire_profiler.config FILE:\n${igblast_constant_ref_files}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    } // path not checked here but below because inside a container
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
    if( ! (meta_path in String || meta_path in GString) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID meta_path PARAMETER IN repertoire_profiler.config FILE:\n${meta_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if(meta_path != "NULL"){
        if( ! (file(meta_path).exists()) ){
            error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID meta_path PARAMETER IN repertoire_profiler.config FILE (DOES NOT EXIST): ${meta_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
        }
    }
    if( ! (meta_seq_names in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID meta_seq_names PARAMETER IN repertoire_profiler.config FILE:\n${meta_seq_names}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (meta_name_replacement in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID meta_name_replacement PARAMETER IN repertoire_profiler.config FILE:\n${meta_name_replacement}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (meta_legend in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID meta_legend PARAMETER IN repertoire_profiler.config FILE:\n${meta_legend}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (germ_tree_kind in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_kind PARAMETER IN repertoire_profiler.config FILE:\n${germ_tree_kind}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (germ_tree_duplicate_seq in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_duplicate_seq PARAMETER IN repertoire_profiler.config FILE:\n${germ_tree_duplicate_seq}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (germ_tree_duplicate_seq == "TRUE" || germ_tree_duplicate_seq == "FALSE") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_duplicate_seq PARAMETER IN repertoire_profiler.config FILE:\n${germ_tree_duplicate_seq}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if( ! (germ_tree_leaf_color in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_leaf_color PARAMETER IN repertoire_profiler.config FILE:\n${germ_tree_leaf_color}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (germ_tree_leaf_shape in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_leaf_shape PARAMETER IN repertoire_profiler.config FILE:\n${germ_tree_leaf_shape}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (germ_tree_leaf_shape =~  /^[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_leaf_shape PARAMETER IN repertoire_profiler.config FILE:\n${germ_tree_leaf_shape}\nMUST BE A POSITIVE INTEGER VALUE\n\n========\n\n"
    }
    if( ! (germ_tree_leaf_size in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_leaf_size PARAMETER IN repertoire_profiler.config FILE:\n${germ_tree_leaf_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (germ_tree_leaf_size =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_leaf_size PARAMETER IN repertoire_profiler.config FILE:\n${germ_tree_leaf_size}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (germ_tree_label_size in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_size PARAMETER IN repertoire_profiler.config FILE:\n${germ_tree_label_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (germ_tree_label_size =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_size PARAMETER IN repertoire_profiler.config FILE:\n${germ_tree_label_size}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (germ_tree_label_hjust in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_hjust PARAMETER IN repertoire_profiler.config FILE:\n${germ_tree_label_hjust}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (germ_tree_label_hjust =~  /^\-{0,1}[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_hjust PARAMETER IN repertoire_profiler.config FILE:\n${germ_tree_label_hjust}\nMUST BE A NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (germ_tree_label_rigth in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_rigth PARAMETER IN repertoire_profiler.config FILE:\n${germ_tree_label_rigth}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (germ_tree_label_rigth == "TRUE" || germ_tree_label_rigth == "FALSE") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_rigth PARAMETER IN repertoire_profiler.config FILE:\n${germ_tree_label_rigth}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if( ! (germ_tree_label_outside in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_outside PARAMETER IN repertoire_profiler.config FILE:\n${germ_tree_label_outside}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (germ_tree_label_outside == "TRUE" || germ_tree_label_outside == "FALSE") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_outside PARAMETER IN repertoire_profiler.config FILE:\n${germ_tree_label_outside}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if( ! (germ_tree_right_margin in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_right_margin PARAMETER IN repertoire_profiler.config FILE:\n${germ_tree_right_margin}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (germ_tree_right_margin =~  /^[0-9]+\.*[0-9]*$/) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_right_margin PARAMETER IN repertoire_profiler.config FILE:\n${germ_tree_right_margin}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
    }
    if( ! (germ_tree_legend in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_legend PARAMETER IN repertoire_profiler.config FILE:\n${germ_tree_legend}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (germ_tree_legend == "TRUE" || germ_tree_legend == "FALSE") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_legend PARAMETER IN repertoire_profiler.config FILE:\n${germ_tree_legend}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
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
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_hole_text PARAMETER IN repertoire_profiler.config FILE:\n${germ_tree_legend}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
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
    if( ! (phylo_tree_heavy in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID phylo_tree_heavy PARAMETER IN repertoire_profiler.config FILE:\n${phylo_tree_heavy}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (phylo_tree_heavy == "true" || phylo_tree_heavy == "false") ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID phylo_tree_heavy PARAMETER IN repertoire_profiler.config FILE:\n${phylo_tree_heavy}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    if( ! (phylo_tree_model_path in String || phylo_tree_model_path in GString) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID phylo_tree_model_path PARAMETER IN repertoire_profiler.config FILE:\n${phylo_tree_model_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }else if( ! (file(phylo_tree_model_path).exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID phylo_tree_model_path PARAMETER IN repertoire_profiler.config FILE (DOES NOT EXIST): ${phylo_tree_model_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
    if( ! (phylo_tree_itolkey in String) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID phylo_tree_itolkey PARAMETER IN repertoire_profiler.config FILE:\n${phylo_tree_itolkey}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
    }
    if( ! (cute_path in String || cute_path in GString) ){
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
    if(igblast_variable_ref_files =~ /^.*IG(K|L)V.*$/){
        print("\n\nWARNING:\nLIGHT CHAIN DETECTED IN THE igblast_variable_ref_files parameter.\nBUT CLONAL GROUPING IS GENERALLY RESTRICTED TO HEAVY CHAIN SEQUENCES, AS THE DIVERSITY OF LIGHT CHAINS IS NOT SUFFICIENT TO DISTINGUISH CLONES WITH REASONABLE CERTAINTY")
    }
    print("\n\nWARNING:\nTO MAKE THE REPERTOIRES AND DONUTS, THE SCRIPT CURRENTLY TAKES THE FIRST ANNOTATION OF THE IMGT ANNOTATION IF SEVERAL ARE PRESENTS IN THE v_call, j_call OR c_call COLUMN OF THE productive_seq.tsv FILE")
    print("\n\n")



    //////// end Checks


    //////// Variable modification


    //////// end Variable modification


    //////// Channels

    // fs_ch define below because can be a .zip file

    //////// end Channels


    //////// files import

    meta_file = file(meta_path) // in variable because a single file. If "NULL", will create a empty file, present in work folders, but that cannot be correctly linked. Thus, if the file has to be redirected into a channel inside a process, it will not work. Thus, in the first process using meta_file, I hard copy the NULL file if required (see below)
    cute_file = file(cute_path) // in variable because a single file
    phylo_tree_model_file  = file(phylo_tree_model_path)
    template_rmd = file(template_rmd_path)

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

    file("${out_path}/files").mkdirs()

    igblast_data_check(
        igblast_organism, 
        igblast_variable_ref_files,
        igblast_constant_ref_files
    )

    igblast(
        fs_ch, 
        igblast_variable_ref_files, 
        igblast_organism, 
        igblast_loci, 
        igblast_aa
    )
    igblast.out.tsv_ch1.count().subscribe{ n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\n0 ANNOTATION SUCCEEDED BY THE igblast PROCESS\n\nCHECK THAT THE igblast_organism, igblast_loci AND igblast_variable_ref_files ARE CORRECTLY SET IN THE repertoire_profiler.config FILE\n\n========\n\n"}}
    tsv_ch2 = igblast.out.tsv_ch1.collectFile(name: "all_igblast_seq.tsv", skip: 1, keepHeader: true)
    igblast.out.aligned_seq_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\n0 ALIGNED SEQ FILES RETURNED BY THE igblast PROCESS\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n"}}
    aligned_seq_ch2 = igblast.out.aligned_seq_ch.collectFile(name: "igblast_aligned_seq_name.tsv")
    aligned_seq_ch2.subscribe{it -> it.copyTo("${out_path}/files")}
    igblast.out.unaligned_seq_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\n0 UNALIGNED SEQ FILES RETURNED BY THE igblast PROCESS\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n"}}
    unaligned_seq_ch2 = igblast.out.unaligned_seq_ch.collectFile(name: "igblast_unaligned_seq_name.tsv") // because an empty file must be present
    unaligned_seq_ch2.subscribe{it -> it.copyTo("${out_path}/files")}
    nb1 = aligned_seq_ch2.countLines() 
    nb2 =  unaligned_seq_ch2.countLines()
    // nb1.view()
    //nb1.view()
    //fs_ch.count().view()
    fs_ch.count().combine(nb1).combine(nb2).subscribe{n,n1,n2 -> if(n != n1 + n2){error "\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nTHE NUMBER OF FILES IN THE igblast_aligned_seq.tsv (${n1}) AND igblast_unaligned_seq_name.tsv (${n2}) IS NOT EQUAL TO THE NUMBER OF SUBMITTED FASTA FILES (${n})\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n"}}
    igblast.out.log_ch.collectFile(name: "igblast_report.log").subscribe{it -> it.copyTo("${out_path}/reports")}


    parseDb_filtering(
        igblast.out.tsv_ch1.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE igblast PROCESS\n\n========\n\n"},
        igblast_aa
    )
    // parseDb_filtering.out.unproductive_ch.count().subscribe{n -> if ( n == 0 ){print "\n\nWARNING: EMPTY unproductive_seq.tsv FILE RETURNED FOLLOWING THE parseDb_filtering PROCESS\n\n"}else{it -> it.copyTo("${out_path}/unproductive_seq.tsv")}} // see https://www.nextflow.io/docs/latest/script.html?highlight=copyto#copy-files
    parseDb_filtering.out.select_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nO PRODUCTIVE SEQUENCE FILES FOLLOWING THE parseDb_filtering PROCESS\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n"}}
    select_ch2 = parseDb_filtering.out.select_ch.collectFile(name: "productive_seq_init.tsv", skip: 1, keepHeader: true)
    parseDb_filtering.out.unselect_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nO UNPRODUCTIVE SEQUENCE FILES FOLLOWING THE parseDb_filtering PROCESS\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n"}} // because an empty file must be present
    unselect_ch2 = parseDb_filtering.out.unselect_ch.collectFile(name: "unproductive_seq.tsv", skip: 1, keepHeader: true)
    unselect_ch2.subscribe{it -> it.copyTo("${out_path}/files")}
    nb1_b = select_ch2.countLines()
    nb2_b = unselect_ch2.countLines()
    nb1_b.subscribe { n -> if ( n == 1 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nO PRODUCTIVE SEQUENCE FOLLOWING THE parseDb_filtering PROCESS\nSEE THE EMPTY productive_seq.tsv FILE AND THE unproductive_seq.tsv FILE IN THE OUTPUT FOLDER\n\n========\n\n"}}
    //nb1_b.view()
    //nb2_b.view()
    //tsv_ch2.countLines().view()
    tsv_ch2.countLines().combine(nb1_b).combine(nb2_b).subscribe{n,n1,n2 -> if(n != n1 + n2 - 1){error "\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nTHE NUMBER OF FILES IN THE productive_seq.tsv (${n1}) AND unproductive_seq.tsv (${n2} - 1) IS NOT EQUAL TO THE NUMBER OF FILES IN all_igblast_seq.tsv (${n})\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n"}} // n1 + n2 - 1 because header counted in both n1 and n2 while only one header in n
    parseDb_filtering.out.parseDb_filtering_log_ch.collectFile(name: "ParseDb_filtering.log").subscribe{it -> it.copyTo("${out_path}/reports")}


    translation(
        parseDb_filtering.out.select_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE parseDb_filtering PROCESS\n\n========\n\n"},
        igblast_aa
    )
    translation.out.translation_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE translation PROCESS\n\n========\n\n"}
    //translation.out.translation_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE translation PROCESS\n\n========\n\n"}}
    translation_ch2 = translation.out.translation_ch.collectFile(name: "translation.tsv", skip: 1, keepHeader: true) // productive file with column sequence_alignment_aa added
    aa_tsv_ch2 = translation.out.aa_tsv_ch.collectFile(name: "aa.tsv", skip: 1, keepHeader: true)
    aa_tsv_ch2.subscribe{it -> it.copyTo("${out_path}/files")}
    translation.out.translation_log_ch.collectFile(name: "translation.log").subscribe{it -> it.copyTo("${out_path}/reports")}



    distToNearest(
        translation_ch2.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE translation PROCESS\n\n========\n\n"},
        igblast_aa,
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
        translation.out.translation_ch,
        meta_file,
        meta_seq_names, 
        meta_name_replacement,
        meta_legend
    )
    seq_name_replacement.out.seq_name_replacement_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE seq_name_replacement PROCESS\n\n========\n\n"}}
    seq_name_replacement_ch2 = seq_name_replacement.out.seq_name_replacement_ch.collectFile(name: "replacement.tsv", skip: 1, keepHeader: true)
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
        clone_distance
    )
    // clone_assignment.out.clone_ch.view()

    nb_failed_clone = clone_assignment.out.failed_clone_ch.countLines() // Assuming there are only the sequence names and no column names on the first line

    split_by_clones(
        clone_assignment.out.clone_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE clone_assignment PROCESS\n\n========\n\n"}
    )



    closest_germline(
        split_by_clones.out.clone_split_ch.flatten(), // flatten split the list into several objects (required for parallelization)
        igblast_organism, 
        igblast_variable_ref_files
    )
    closest_germline.out.closest_log_ch.collectFile(name: "closest_germline.log").subscribe{it -> it.copyTo("${out_path}/reports")} // 



    mutation_load(
        closest_germline.out.closest_ch
    )
    //mutation_load.out.mutation_load_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE mutation_load PROCESS\n\n========\n\n"}}
    mutation_load.out.mutation_load_log_ch.collectFile(name: "mutation_load.log").subscribe{it -> it.copyTo("${out_path}/reports")} // 

    clone_assigned_seq = mutation_load.out.mutation_load_ch.collectFile(name: "clone_assigned_seq.tsv", skip: 1, keepHeader: true)
    nb_clone_assigned = clone_assigned_seq.countLines() - 1 // Minus 1 because 1st line = column names
    clone_assigned_seq.subscribe{it -> it.copyTo("${out_path}/files")}



    get_germ_tree(
        mutation_load.out.mutation_load_ch,
        meta_file, // first() because get_germ_tree process is a parallele one and because meta_file is single
        cute_file, 
        clone_nb_seq,
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
    no_germ_tree_ch2.subscribe{it -> it.copyTo("${out_path}/files")}

    get_germ_tree.out.germ_tree_ch.count().subscribe { n -> if ( n == 0 ){
        print("\n\nWARNING: NO SEQUENCES IN TREES FOLLOWING THE get_germ_tree PROCESS -> EMPTY germ_tree_seq.tsv FILE RETURNED\n\n")
    }}
    germ_tree_ch2 = get_germ_tree.out.germ_tree_ch.collectFile(name: "germ_tree_seq.tsv", skip: 1, keepHeader: true)
    germ_tree_ch2.subscribe{it -> it.copyTo("${out_path}/files")}

    get_germ_tree.out.no_cloneID_ch.count().subscribe { n -> if ( n == 0 ){
        print("\n\nWARNING: ALL SEQUENCES IN CLONAL GROUP FOLLOWING THE get_germ_tree PROCESS -> EMPTY germ_tree_dismissed_clone_id.tsv FILE RETURNED\n\n")
    }}
    no_cloneID_ch2 = get_germ_tree.out.no_cloneID_ch.collectFile(name: "germ_tree_dismissed_clone_id.tsv")
    no_cloneID_ch2.subscribe{it -> it.copyTo("${out_path}/files")}

    get_germ_tree.out.cloneID_ch.count().subscribe { n -> if ( n == 0 ){
        print("\n\nWARNING: NO CLONAL GROUP FOLLOWING THE get_germ_tree PROCESS -> EMPTY germ_tree_clone_id.tsv and germ_tree.pdf FILES RETURNED\n\n")
    }}
    cloneID_ch2 = get_germ_tree.out.cloneID_ch.collectFile(name: "germ_tree_clone_id.tsv")
    cloneID_ch2.subscribe{it -> it.copyTo("${out_path}/files")}

    get_germ_tree.out.get_germ_tree_log_ch.collectFile(name: "get_germ_tree.log").subscribe{it -> it.copyTo("${out_path}/reports")} // 



    germ_tree_vizu(
        rdata_germ_tree_ch2,
        germ_tree_kind,
        clone_nb_seq,
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
        clone_assigned_seq, // may be add .first()
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
    germ_tree_dup_seq_not_displayed_ch2.subscribe{it -> it.copyTo("${out_path}/files")}

    tempo1_ch = Channel.of("all", "annotated", "tree") // 1 channel with 3 values (not list)
    tempo2_ch = data_assembly.out.productive_ch.mix(data_assembly.out.productive_ch.mix(germ_tree_ch2)) // 1 channel with 3 paths (do not use flatten() -> not list)
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
        cute_file
    )
    donut.out.donut_pdf_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE donut PROCESS\n\n========\n\n"}
    donut.out.donut_pdf_ch.count().subscribe { n -> if ( n == 0 ){
        print("\n\nWARNING: EMPTY OUTPUT FOLLOWING THE donut PROCESS -> NO DONUT RETURNED\n\n")}}
    donut_pdf_ch2 = donut.out.donut_pdf_ch.collect()

    donut.out.donut_tsv_ch.count().subscribe { n -> if ( n == 0 ){
        print("\n\nWARNING: -> NO donut_stats.tsv FILE RETURNED\n\n")
    }}
    donut_tsv_ch2 = donut.out.donut_tsv_ch.collectFile(name: "donut_stats.tsv", skip: 1, keepHeader: true)
    donut_tsv_ch2.subscribe{it -> it.copyTo("${out_path}/files")}



    donut_assembly(
        donut_pdf_ch2
    )
    //donut_assembly.out.donut_assembly_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE donut_assembly PROCESS\n\n========\n\n"}
    donut_assembly.out.donut_assembly_ch.count().subscribe { n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE donut_assembly PROCESS\n\n========\n\n"}}


    Reformat(
        aa_tsv_ch2
    )
    fasta = Reformat.out.aa_fasta_ch


    DefineGroups(
        fasta,
        select_ch2
    )
    fastagroups = DefineGroups.out.groups_ch.flatten()


    Align(
        fastagroups,
        phylo_tree_heavy
    )
    align = Align.out.aligned_groups_ch


    NbSequences(
        align
    )
    filtered = NbSequences.out.nb_out.filter{it[0].toInteger()>=3}.map{it->it[1]}

    Tree(
        filtered,
        phylo_tree_model_file
    )
    tree = Tree.out.tree_file


    ProcessMeta(
        meta_file
    )
    itolmeta = ProcessMeta.out.itol_out


    ITOL(
        tree,
        itolmeta,
        phylo_tree_itolkey
    )


    print_report(
        template_rmd,
        nb_input,
        nb1,
        nb2,
        nb1_b - 1, // Minus 1 because the 1st line = the column names
        nb2_b - 1, // Minus 1 because the 1st line = the column names
        nb_clone_assigned,
        nb_failed_clone,
        donut.out.donuts_png.collect(),
        repertoire.out.repertoire_png.collect()
    )

    



    backup(
        config_file, 
        log_file
    )

}

    //////// end Main


//////// end Processes
