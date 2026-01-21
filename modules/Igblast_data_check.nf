

// next process is for repertoire
process Igblast_data_check { // cannot be Igblast_data_check outside of process because the files are present in the docker
    label 'immcantation'
    //publishDir path: "${out_path}", mode: 'copy', pattern: "{*.tsv}", overwrite: false

    cache 'true'

    input:
    val igblast_organism
    val igblast_v_ref_files
    val igblast_d_ref_files
    val igblast_j_ref_files
    val igblast_constant_ref_files

    output:
    path "*.tsv", emit: allele_names_tsv_all_ch // all tsv files with allele names

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    REPO_PATH_VAR="/usr/local/share/germlines/imgt/${igblast_organism}/vdj" # path where the imgt .fasta reference seq files are in the docker container
    REPO_PATH_CONST="/usr/local/share/germlines/imgt/${igblast_organism}/constant" # path where the imgt .fasta reference seq files are in the docker container
    V_FILES=\$(awk -v var1="${igblast_v_ref_files}" -v var2="\${REPO_PATH_VAR}" 'BEGIN{ORS=" " ; split(var1, array1, " ") ; for (key in array1) {print var2"/"array1[key]}}') # assemble files with their path
    if [[ "${igblast_d_ref_files}" != "NULL" ]] ; then
        D_FILES=\$(awk -v var1="${igblast_d_ref_files}" -v var2="\${REPO_PATH_VAR}" 'BEGIN{ORS=" " ; split(var1, array1, " ") ; for (key in array1) {print var2"/"array1[key]}}') # assemble files with their path
    fi
    J_FILES=\$(awk -v var1="${igblast_j_ref_files}" -v var2="\${REPO_PATH_VAR}" 'BEGIN{ORS=" " ; split(var1, array1, " ") ; for (key in array1) {print var2"/"array1[key]}}') # assemble files with their path
    CONST_FILES=\$(awk -v var1="${igblast_constant_ref_files}" -v var2="\${REPO_PATH_CONST}" 'BEGIN{ORS=" " ; split(var1, array1, " ") ; for (key in array1) {print var2"/"array1[key]}}') # assemble files with their path
    for i1 in \$V_FILES ; do
        if [[ ! -e \${i1} ]] ; then
            echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nFILE DOES NOT EXISTS:\\n\${i1}\\n\\nINDICATED PATH:\\n\${REPO_PATH_VAR}\\n\\nCONTAINS:\\n"
            ls -la -w 1 \${REPO_PATH_VAR}
            echo -e "\\n\\n========\\n\\n"
            exit 1
        else
            FILENAME=\$(basename -- "\${i1}") # recover a file name without path
            FILENAME="\${FILENAME%.*}" # remove extension
            grep -E '^>.*\$' \${i1} | cut -f2 -d'|' > \${FILENAME}.tsv.v # detect line starting by > and extract the 2nd field after cutting by |
            cp \${FILENAME}.tsv.v \${FILENAME}.tsv
        fi
    done
    if [[ "${igblast_d_ref_files}" != "NULL" ]] ; then
        for i1 in \$D_FILES ; do
            if [[ ! -e \${i1} ]] ; then
                echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nFILE DOES NOT EXISTS:\\n\${i1}\\n\\nINDICATED PATH:\\n\${REPO_PATH_VAR}\\n\\nCONTAINS:\\n"
                ls -la -w 1 \${REPO_PATH_VAR}
                echo -e "\\n\\n========\\n\\n"
                exit 1
            else
                FILENAME=\$(basename -- "\${i1}") # recover a file name without path
                FILENAME="\${FILENAME%.*}" # remove extension
                grep -E '^>.*\$' \${i1} | cut -f2 -d'|' > \${FILENAME}.tsv.d # detect line starting by > and extract the 2nd field after cutting by |
                cp \${FILENAME}.tsv.d \${FILENAME}.tsv
            fi
        done
    fi
    for i1 in \$J_FILES ; do
        if [[ ! -e \${i1} ]] ; then
            echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nFILE DOES NOT EXISTS:\\n\${i1}\\n\\nINDICATED PATH:\\n\${REPO_PATH_VAR}\\n\\nCONTAINS:\\n"
            ls -la -w 1 \${REPO_PATH_VAR}
            echo -e "\\n\\n========\\n\\n"
            exit 1
        else
            FILENAME=\$(basename -- "\${i1}") # recover a file name without path
            FILENAME="\${FILENAME%.*}" # remove extension
            grep -E '^>.*\$' \${i1} | cut -f2 -d'|' > \${FILENAME}.tsv.j # detect line starting by > and extract the 2nd field after cutting by |
            cp \${FILENAME}.tsv.j \${FILENAME}.tsv
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
            grep -E '^>.*\$' \${i1} | cut -f2 -d'|' > \${FILENAME}.tsv.c # detect line starting by > and extract the 2nd field after cutting by |
            cp \${FILENAME}.tsv.c \${FILENAME}.tsv
        fi
    done
    """
    // write ${} between "" to make a single argument when the variable is made of several values separated by a space. Otherwise, several arguments will be considered
}

