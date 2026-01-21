

process Clone_assignment {
    label 'immcantation'
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{*.log}", overwrite: false
    publishDir path: "${out_path}/tsv", mode: 'copy', pattern: "failed_clone_assigned_seq.tsv", overwrite: false
    cache 'true'

    input:
    path wanted_ch // no parallelization
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
    # if [[ -s ${wanted_ch} ]]; then # see above for -s # no need anymore
    FILENAME=\$(basename -- ${wanted_ch}) # recover a file name without path
    FILE=\${FILENAME%.*} # file name without extension
    DefineClones.py -d ${wanted_ch} --act ${clone_strategy} --model ${clone_model} --norm ${clone_normalize} --dist ${clone_distance} --fail |& tee -a clone_assignment.log
    shopt -s nullglob
    # files=(\${FILE}_clone-{pass,fail}.tsv) # expend the two possibilities but assign only the existing ones
    # cp -Lr "\${files[0]}" "tempo.tsv" # copy whatever failed or succeeded
    cat > run.R <<'RSCRIPT'
        options(warning.length = 8170)
        args = commandArgs(trailingOnly=TRUE)
        db <- read.table(args[1], sep = "\\t", header = TRUE)
        wanted <- read.table(args[2], sep = "\\t", header = TRUE)
        output_name <- args[3]
        # put back meta_legend column name as in wanted_seq.tsv, beacause in lowercase
        if("${meta_file}" != "NULL" & "${meta_legend}" != "NULL" ){
            tempo_log <- names(db) == tolower("${meta_legend}")
            if(any(tempo_log, na.rm = TRUE)){
                names(db)[tempo_log] <- "${meta_legend}"
            }
        }
        # end put back meta_legend column names as in wanted_seq.tsv, beacause in lowercase
        # reorder as in wanted_seq.tsv
        if( ! all(names(wanted) %in% names(db))){
            stop(paste0("\\n\\n================\\n\\nERROR IN Clone_assignment PROCESS.\\nNAMES OF wanted_seq.tsv SHOULD ALL BE IN THE OUTPUT .tsv FILE OF DefineClones.py.\\nNAMES OF wanted_seq.tsv:\\n", paste0(sort(names(wanted)), collapse = " "), "\\nNAMES OF THE OUTPUT:\\n", paste0(sort(names(db)), collapse = " "), "\\n\\n================\\n\\n"), call. = FALSE)
        }else{
            tempo_log <- names(db) %in% names(wanted)
            tempo_db1 <- db[tempo_log]
            tempo_db2 <- db[ ! tempo_log]
            tempo_db1 <- tempo_db1[ , match(names(wanted), names(tempo_db1))] # reorder as in wanted_seq.tsv
            db2 <- data.frame(tempo_db1, tempo_db2)
            write.table(db2, file = output_name, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\\t")
        }
        # end reorder as in wanted_seq.tsv
RSCRIPT
# must be at the start of the line (no spaces/tabs) and with no trailing spaces.
    if [ -s \${FILE}_clone-fail.tsv ]; then # see above for -s
        Rscript run.R \${FILE}_clone-fail.tsv ${wanted_ch} "failed_clone_assigned_seq.tsv" |& tee -a clone_assignment.log
    else
        echo -e "\\n\\nNOTE: EMPTY failed_clone_assigned_seq.tsv FILE RETURNED FOLLOWING THE Clone_assignment PROCESS\\n\\n" |& tee -a clone_assignment.log
        head -1 \${FILE}_clone-pass.tsv | cat > failed_clone_assigned_seq.tsv # header kept
    fi
    if [ -s \${FILE}_clone-pass.tsv ] ; then # see above for -s
        Rscript run.R \${FILE}_clone-pass.tsv ${wanted_ch} "\${FILE}_clone-pass.tsv" |& tee -a clone_assignment.log
    else
        echo -e "\\n\\nNOTE: EMPTY *_clone-pass.tsv FILE RETURNED FOLLOWING THE Clone_assignment PROCESS\\n\\n" |& tee -a clone_assignment.log
        head -1 failed_clone_assigned_seq | cat > \${FILE}_clone-pass.tsv # header kept
        # set -o pipefail
        # echo -e "\\n\\n========\\n\\nINTERNAL ERROR IN NEXTFLOW EXECUTION\\n\\nOUTPUT FILE OF DefineClones.py IN THE Clone_assignment PROCESS SHOULD RETURN *_clone-pass.tsv.\\nCHECK THE clone_assignment.log IN THE report FOLDER INSIDE THE OUTPUT FOLDER\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n"
        # exit 1
    fi
    """
}