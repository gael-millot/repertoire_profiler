

// This process adds the germline sequences for each v,d,j gene specified in the germline_(v|d|j)_call columns present in the tsv file
// Inputs:
//      - closest_ch: tsv file containing at least germline_(v|d|j)_call column (otherwise error), sequences already regrouped in clonal groups
//      - igblast_organism: value specified in nextflow.config
//  contains germline ref sequences of IMGT database
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
    VDJ_FILES_INI="${igblast_variable_ref_files}"
    VDJ_FILES_INI="\${VDJ_FILES_INI// NULL / }" # remove the string NULL if exists
    VDJ_FILES=\$(awk -v var1="\${VDJ_FILES_INI}" -v var2="\${REPO_PATH}" 'BEGIN{ORS=" " ; split(var1, array1, " ") ; for (key in array1) {print var2"/"array1[key]}}')
    GermlineSequences.py -i \$FILENAME -r \${VDJ_FILES} |& tee -a GermlineSequences.log
    """
}


