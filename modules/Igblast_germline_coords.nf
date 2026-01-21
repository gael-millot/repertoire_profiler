// recover the coordinates in the germline sequence
process Igblast_germline_coords {
    label 'immcantation'
    //publishDir path: "${out_path}", mode: 'copy', overwrite: false
    cache 'true'

    input:
    tuple path(tsv), path(fasta) // parallelization expected (for each fasta file)
    val igblast_organism
    val igblast_loci

    output:
    path "with_germ_coords.tsv", emit: germline_coords_ch
    path "*.log", emit: log_ch

    script:
    """
    #!/bin/bash -ue
    # See https://changeo.readthedocs.io/en/stable/tools/AssignGenes.html for the details
    AssignGenes.py igblast -s ${fasta} -b /usr/local/share/igblast --organism ${igblast_organism} --loci ${igblast_loci} --format airr |& tee -a igblast_germline_report.log # output is *_igblast.tsv
    # test that the AssignGenes.py igblast --format airr is ok
        expected="sequence_id	sequence	sequence_aa	locus	stop_codon	vj_in_frame	v_frameshift	productive	rev_comp	complete_vdj	d_frame	v_call	d_call	j_call	c_call	sequence_alignment	germline_alignment	sequence_alignment_aa	germline_alignment_aa	v_alignment_start	v_alignment_end	d_alignment_start	d_alignment_end	j_alignment_start	j_alignment_end	c_alignment_start	c_alignment_end	v_sequence_alignment	v_sequence_alignment_aa	v_germline_alignment	v_germline_alignment_aa	d_sequence_alignment	d_sequence_alignment_aa	d_germline_alignment	d_germline_alignment_aa	j_sequence_alignment	j_sequence_alignment_aa	j_germline_alignment	j_germline_alignment_aa	c_sequence_alignment	c_sequence_alignment_aa	c_germline_alignment	c_germline_alignment_aa	fwr1	fwr1_aa	cdr1	cdr1_aa	fwr2	fwr2_aa	cdr2	cdr2_aa	fwr3	fwr3_aa	fwr4	fwr4_aa	cdr3	cdr3_aa	junction	junction_length	junction_aa	junction_aa_length	v_score	d_score	j_score	c_score	v_cigar	d_cigar	j_cigar	c_cigar	v_support	d_support	j_support	c_support	v_identity	d_identity	j_identity	c_identity	v_sequence_start	v_sequence_end	v_germline_start	v_germline_end	d_sequence_start	d_sequence_end	d_germline_start	d_germline_end	j_sequence_start	j_sequence_end	j_germline_start	j_germline_end	c_sequence_start	c_sequence_end	c_germline_start	c_germline_end	fwr1_start	fwr1_end	cdr1_start	cdr1_end	fwr2_start	fwr2_end	cdr2_start	cdr2_end	fwr3_start	fwr3_end	fwr4_start	fwr4_end	cdr3_start	cdr3_end	np1	np1_length	np2	np2_length"
    if [[ -s *_igblast.tsv ]]; then # -s means --format airr has worked
        read -r firstline < *_igblast.tsv
        if [[ "\$firstline" != "\$expected" ]]; then
            echo -e "\\n\\n========\\n\\nINTERNAL ERROR IN NEXTFLOW EXECUTION\\n\\nTHE COMANND AssignGenes.py igblast --format airr DOES NOT RETURN THE EXPECTED COLUMN NAMES\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n"
            echo -e "COLUMN NAMES:\\n \${firstline} \\n\\n"
            echo -e "EXPECTED COLUMN NAMES:\\n \${expected}"
            echo -e "\\n\\n========\\n\\n"
            exit 1
        fi
        # end test that the AssignGenes.py igblast --format airr is ok
    else
        echo -e "\\n\\n========\\n\\nINTERNAL ERROR IN NEXTFLOW EXECUTION\\n\\nTHE COMANND AssignGenes.py igblast --format airr DOES NOT RETURN THE EXPECTED tsv FILE\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n"
        # exit 1
    fi
    Rscript -e '
        args <- commandArgs(trailingOnly = TRUE)  # recover arguments written after the call of the Rscript
        coords <- read.table(args[1], sep = "\\t", header = TRUE)
        tsv <- read.table(args[2], sep = "\\t", header = TRUE)
        col_names <- c("sequence_alignment", "sequence_alignment_aa", "v_sequence_start", "v_sequence_end", "d_sequence_start", "d_sequence_end", "j_sequence_start", "j_sequence_end", "c_sequence_start", "c_sequence_end", "fwr1_start", "fwr1_end", "cdr1_start", "cdr1_end", "fwr2_start", "fwr2_end", "cdr2_start", "cdr2_end", "fwr3_start", "fwr3_end", "cdr3_start", "cdr3_end", "fwr4_start", "fwr4_end")
        db <- coords[, match(col_names, names(coords))]
        names(db) <- sub(x = names(db), pattern = "sequence_", "")
        names(db) <- paste0("clonal_germline_", names(db))
        names(db)[names(db) == "clonal_germline_alignment"] <- "clonal_germline_alignment_igblast_airr"
        names(db)[names(db) == "clonal_germline_alignment_aa"] <- "clonal_germline_alignment_aa_igblast_airr"
        tsv <- data.frame(tsv, db)
        write.table(tsv, file = "with_germ_coords.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\\t")
    ' *_igblast.tsv *_germ-pass.tsv |& tee -a igblast_germline_report.log
    """
    // write ${} between "" to make a single argument when the variable is made of several values separated by a space. Otherwise, several arguments will be considered
}