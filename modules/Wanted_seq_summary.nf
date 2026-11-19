// Select the main columns of the wanted_seq.tsv file (same rows, same values): the columns used by the other outputs (repertoires, donuts, distance histograms, alignments, trees) and those describing each sequence (name, chain, V(D)J and C annotation, CDRs, sequences)
// Inputs:
//      wanted_ch: wanted_seq.tsv file (all the wanted sequences, in a single file)
//      meta_legend: meta_legend parameter of the nextflow.config file, i.e., name of the metadata column present in wanted_seq.tsv if the meta_path and meta_legend parameters are not "NULL"
// Outputs:
//      "wanted_seq_summary.tsv": sequence_id, initial_sequence_id (optional), meta_legend column (optional), locus, v_call, d_call, j_call, c_call, v_gene, j_gene, c_gene, isotype_class, complete_vdj, v_identity, cdr1_aa, cdr2_aa, cdr3_aa, junction_aa, junction_aa_length, dist_nearest, sequence_alignment_aa and sequence_alignment columns of wanted_seq.tsv, in this order. Optional columns are skipped if absent from wanted_seq.tsv
process Wanted_seq_summary {
    label 'bash'
    publishDir path: "${out_path}/tsv", mode: 'copy', pattern: "wanted_seq_summary.tsv", overwrite: false
    cache 'true'

    input:
    path wanted_ch // no parallelization
    val meta_legend

    output:
    path "wanted_seq_summary.tsv", emit: wanted_summary_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    if [[ ! -s ${wanted_ch} ]]; then # -s means "exists and non empty". Thus, return FALSE is the file does not exists or is empty
        echo -e "\\n\\n========\\n\\nINTERNAL ERROR IN NEXTFLOW EXECUTION\\n\\nEMPTY ${wanted_ch} FILE (NO HEADER) IN THE Wanted_seq_summary PROCESS.\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n" >&2
        exit 1
    fi
    # columns recovered by name in the header (first line), in the order below. initial_sequence_id and the meta_legend column are optional: only present in wanted_seq.tsv if the meta_path and meta_name_replacement (respectively meta_legend) parameters are not "NULL"
    awk -F '\\t' -v OFS='\\t' -v meta_legend="${meta_legend}" '
        NR == 1{
            nb_col = split("sequence_id initial_sequence_id meta_legend locus v_call d_call j_call c_call v_gene j_gene c_gene isotype_class complete_vdj v_identity cdr1_aa cdr2_aa cdr3_aa junction_aa junction_aa_length dist_nearest sequence_alignment_aa sequence_alignment", col_names, " ")
            nb_kept = 0
            for(i1 = 1; i1 <= nb_col; i1++){
                optional = (col_names[i1] == "initial_sequence_id" || col_names[i1] == "meta_legend")
                if(col_names[i1] == "meta_legend"){
                    col_names[i1] = meta_legend # name of the metadata column
                }
                pos = 0
                for(i2 = 1; i2 <= NF; i2++){
                    if(\$i2 == col_names[i1]){
                        pos = i2
                        break
                    }
                }
                if(pos == 0 && ! optional){
                    print "\\n\\n========\\n\\nINTERNAL ERROR IN NEXTFLOW EXECUTION\\n\\nNO " col_names[i1] " COLUMN IN THE ${wanted_ch} FILE OF THE Wanted_seq_summary PROCESS.\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n" > "/dev/stderr"
                    exit 1
                }
                if(pos > 0 && ! (pos in kept)){ # a column is not selected twice (if meta_legend is the name of another selected column)
                    nb_kept++
                    kept[pos] = 1
                    col_pos[nb_kept] = pos
                }
            }
        }
        {
            line = \$(col_pos[1])
            for(i1 = 2; i1 <= nb_kept; i1++){
                line = line OFS \$(col_pos[i1])
            }
            print line
        }
    ' ${wanted_ch} > wanted_seq_summary.tsv
    """
}
