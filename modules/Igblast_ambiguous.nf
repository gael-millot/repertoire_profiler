// Select the sequences of the igblast_seq.tsv file with more than one value (comma separated, e.g., IGKV1-39*01,IGKV1D-39*01) in the v_call, d_call, j_call or c_call column
// Inputs:
//      igblast_seq_ch: igblast_seq.tsv file (all the sequences annotated by igblast, in a single file)
// Outputs:
//      "ambiguous_igblast_seq.tsv": igblast_seq.tsv_line_nb (number of the line in the igblast_seq.tsv file, header excluded), sequence_id, v_call, d_call, j_call and c_call columns of these sequences. Header only if no such sequence
process Igblast_ambiguous {
    label 'bash'
    publishDir path: "${out_path}/tsv", mode: 'copy', pattern: "ambiguous_igblast_seq.tsv", overwrite: false
    cache 'true'

    input:
    path igblast_seq_ch // no parallelization

    output:
    path "ambiguous_igblast_seq.tsv", emit: ambiguous_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    if [[ ! -s ${igblast_seq_ch} ]]; then # -s means "exists and non empty". Thus, return FALSE is the file does not exists or is empty
        echo -e "\\n\\n========\\n\\nINTERNAL ERROR IN NEXTFLOW EXECUTION\\n\\nEMPTY ${igblast_seq_ch} FILE (NO HEADER) IN THE Igblast_ambiguous PROCESS.\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n" >&2
        exit 1
    fi
    # columns recovered by name in the header (first line). Then, lines with a comma in the v_call, d_call, j_call or c_call column, with their line number in the file, header excluded (NR - 1)
    awk -F '\\t' -v OFS='\\t' '
        NR == 1{
            nb_col = split("sequence_id v_call d_call j_call c_call", col_names, " ")
            for(i1 = 1; i1 <= nb_col; i1++){
                for(i2 = 1; i2 <= NF; i2++){
                    if(\$i2 == col_names[i1]){
                        col_pos[i1] = i2
                    }
                }
                if( ! (i1 in col_pos)){
                    print "\\n\\n========\\n\\nINTERNAL ERROR IN NEXTFLOW EXECUTION\\n\\nNO " col_names[i1] " COLUMN IN THE ${igblast_seq_ch} FILE OF THE Igblast_ambiguous PROCESS.\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n" > "/dev/stderr"
                    exit 1
                }
            }
            print "igblast_seq.tsv_line_nb", "sequence_id", "v_call", "d_call", "j_call", "c_call"
            next
        }
        \$(col_pos[2]) ~ /,/ || \$(col_pos[3]) ~ /,/ || \$(col_pos[4]) ~ /,/ || \$(col_pos[5]) ~ /,/{
            print NR - 1, \$(col_pos[1]), \$(col_pos[2]), \$(col_pos[3]), \$(col_pos[4]), \$(col_pos[5])
        }
    ' ${igblast_seq_ch} > ambiguous_igblast_seq.tsv
    """
}
