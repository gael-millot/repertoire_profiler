workflow CheckVariables {

//////// Checks
//// check of the bin folder
tested_files_bin = [
    "AB_model", 
    "abalign_CLI_help.txt", 
    "alignments_viz.rmd", 
    "circos.R", 
    "circos_data_prep.R", 
    "cute_little_R_functions_v12.8.R", 
    "defineGroups.pl", 
    "donut.R", 
    "fields_not_kept.txt", 
    "germ_tree_vizu.R", 
    "GermlineSequences.py", 
    "get_germ_tree.R", 
    "Gff.R", 
    "Gff_imgt.R", 
    "histogram.R", 
    "parse_coordinates.R", 
    "repertoire.R", 
    "repertoire_profiler_template.rmd", 
    "trimtranslate.sh", 
    "Tsv2fasta.R"
    ]
for(i1 in tested_files_bin){
    if( ! (file("${projectDir}/bin/${i1}").exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nTHE ${i1} FILE MUST BE PRESENT IN THE ./bin FOLDER, WHERE THE main.nf file IS PRESENT\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
}
//// end check of the bin folder

//// check of the modules folder
tested_files_modules = [
    "add_dotted_coord.nf", 
    "closest_germline.nf", 
    "data_assembly.nf", 
    "gff.nf", 
    "gff_imgt.nf", 
    "igblast_germline.nf", 
    "igblast_query.nf", 
    "mutation_load_germ_genes.nf", 
    "parseDb_filtering.nf", 
    "print_alignment.nf", 
    "subscribe_helpers.nf", 
    "trim_translate.nf"
    ]
for(i1 in tested_files_modules){
    if( ! (file("${projectDir}/modules/${i1}").exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nTHE ${i1} FILE MUST BE PRESENT IN THE ./modules FOLDER, WHERE THE main.nf file IS PRESENT\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
}
//// end check of the modules folder

//// check of config file parameters
// Data
if( ! (sample_path in String || sample_path in GString) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN nextflow.config FILE:\n${sample_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (file(sample_path).exists()) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID sample_path PARAMETER IN nextflow.config FILE (DOES NOT EXIST): ${sample_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
}
if( ! (meta_path in String || meta_path in GString) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID meta_path PARAMETER IN nextflow.config FILE:\n${meta_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if(meta_path != "NULL"){
    if( ! (file(meta_path).exists()) ){
        error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID meta_path PARAMETER IN nextflow.config FILE (DOES NOT EXIST): ${meta_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
    }
}
if( ! (meta_seq_names in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID meta_seq_names PARAMETER IN nextflow.config FILE:\n${meta_seq_names}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! (meta_name_replacement in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID meta_name_replacement PARAMETER IN nextflow.config FILE:\n${meta_name_replacement}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}

// Ig annotation
if( ! (igblast_organism in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_organism PARAMETER IN nextflow.config FILE:\n${igblast_organism}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (igblast_organism =~ /^(mouse|human|rabbit|rat|rhesus_monkey)$/) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_organism PARAMETER IN nextflow.config FILE:\n${igblast_organism}\nMUST BE EITHER \"mouse\", \"human\", \"rabbit\", \"rat\" OR \"rhesus_monkey\"\n\n========\n\n"
}else if( ! (igblast_organism =~ /^(mouse|human)$/)){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_organism PARAMETER IN nextflow.config FILE: ${igblast_organism}\n\nTHE repertoire PROCESS CURRENTLY ONLY SUPPORTS mouse AND human SPECIES\nTHEREFORE igblast_organism MUST BE EITHER \"mouse\" OR \"human\"\n\n========\n\n"
}
if( ! (igblast_loci in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_loci PARAMETER IN nextflow.config FILE:\n${igblast_loci}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (igblast_loci == "ig" || igblast_loci == "tr") ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_loci PARAMETER IN nextflow.config FILE:\n${igblast_loci}\nMUST BE EITHER \"ig\" OR \"tr\"\n\n========\n\n"
}
if( ! (igblast_heavy_chain in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_heavy_chain PARAMETER IN nextflow.config FILE:\n${igblast_heavy_chain}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (igblast_heavy_chain == "TRUE" || igblast_heavy_chain == "FALSE") ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_heavy_chain PARAMETER IN nextflow.config FILE:\n${igblast_heavy_chain}\nMUST BE EITHER \"TRUE\" OR \"FALSE\".\n\n========\n\n"
}
if( ! (igblast_lambda_chain in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_lambda_chain PARAMETER IN nextflow.config FILE:\n${igblast_lambda_chain}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (igblast_lambda_chain == "TRUE" || igblast_lambda_chain == "FALSE") ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_lambda_chain PARAMETER IN nextflow.config FILE:\n${igblast_lambda_chain}\nMUST BE EITHER \"TRUE\" OR \"FALSE\".\n\n========\n\n"
}
if( ! (igblast_kappa_chain in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_kappa_chain PARAMETER IN nextflow.config FILE:\n${igblast_kappa_chain}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (igblast_kappa_chain == "TRUE" || igblast_kappa_chain == "FALSE") ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igblast_kappa_chain PARAMETER IN nextflow.config FILE:\n${igblast_kappa_chain}\nMUST BE EITHER \"TRUE\" OR \"FALSE\".\n\n========\n\n"
}
// Checking of studied chain coherence
if (igblast_heavy_chain == "TRUE"  && (igblast_lambda_chain == "TRUE" || igblast_kappa_chain == "TRUE" )) {
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID PARAMETERS IN nextflow.config FILE:\nFOR THE MOMENT, HEAVY (igblast_heavy_chain) AND LIGHT (igblast_lambda_chain, igblast_kappa_chain) CHAIN LOCI CANNOT BOTH BE TRUE.\nHERE ARE THEIR CURRENT VALUES: \nigblast_heavy_chain: ${igblast_heavy_chain}\nigblast_lambda_chain: ${igblast_lambda_chain}\nigblast_kappa_chain: ${igblast_kappa_chain}\n\n========\n\n"
}
if (!igblast_heavy_chain && !igblast_lambda_chain && !igblast_kappa_chain) {
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID PARAMETERS IN nextflow.config FILE:\nAT LEAST ONE OF igblast_heavy_chain, igblast_lambda_chain OR igblast_kappa_chain MUST BE TRUE\nHERE ARE THEIR CURRENT VALUES: \nigblast_heavy_chain: ${igblast_heavy_chain}\nigblast_lambda_chain: ${igblast_lambda_chain}\nigblast_kappa_chain: ${igblast_kappa_chain}\n\n========\n\n"
}
//if ([heavy, lambda, kappa].count { it } > 2) {
    //error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID PARAMETERS IN nextflow.config FILE:\nONLY ONE OF igblast_heavy_chain AND LIGHT CHAINS (igblast_lambda_chain, igblast_kappa_chain) CAN BE TRUE AT A TIME\nHERE ARE THEIR CURRENT VALUES: \nigblast_heavy_chain: ${igblast_heavy_chain}\nigblast_lambda_chain: ${igblast_lambda_chain}\nigblast_kappa_chain: ${igblast_kappa_chain}\n\n========\n\n"
//}

// Clonal groups (clustering) and mutation load
if( ! (clone_strategy in String) ){
            error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_strategy PARAMETER IN nextflow.config FILE:\n${clone_strategy}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (clone_strategy == "first" || clone_strategy == "set") ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_strategy PARAMETER IN nextflow.config FILE:\n${clone_strategy}\nMUST BE EITHER \"first\" OR \"set\"\n\n========\n\n"
}
if( ! (clone_model in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_model PARAMETER IN nextflow.config FILE:\n${clone_model}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (clone_model == "ham" || clone_model == "aa" || clone_model == "hh_s1f" || clone_model == "hh_s5f" || clone_model == "mk_rs1nf" || clone_model == "mk_rs5nf" || clone_model == "m1n_compat" || clone_model == "hs1f_compat") ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_model PARAMETER IN nextflow.config FILE:\n${clone_model}\nMUST BE EITHER \"ham\", \"aa\", \"hh_s1f\", \"hh_s5f\", \"mk_rs1nf\", \"mk_rs5nf\", \"m1n_compat\", \"hs1f_compat\"\n\n========\n\n"
}
if( ! (clone_normalize in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_normalize PARAMETER IN nextflow.config FILE:\n${clone_normalize}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (clone_normalize == "len" || clone_normalize == "none") ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_normalize PARAMETER IN nextflow.config FILE:\n${clone_normalize}\nMUST BE EITHER \"len\" OR \"none\"\n\n========\n\n"
}
if( ! (clone_distance in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_distance PARAMETER IN nextflow.config FILE:\n${clone_distance}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (clone_distance =~ /^((1)|(0)|(0\.[0-9]*))$/) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_distance PARAMETER IN nextflow.config FILE:\n${clone_distance}\nMUST BE A POSITIVE PROPORTION VALUE\n\n========\n\n"
}
if( ! (clone_germline_kind in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_germline_kind PARAMETER IN nextflow.config FILE:\n${clone_germline_kind}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (clone_germline_kind =~ /^(dmask|full|vonly)$/) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_germline_kind PARAMETER IN nextflow.config FILE:\n${clone_germline_kind}\nMUST BE EITHER \"dmask\", \"full\", \"vonly\".\n\n========\n\n"
}

if( ! (clone_mut_obs_seq in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_mut_obs_seq PARAMETER IN nextflow.config FILE:\n${clone_mut_obs_seq}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! (clone_mut_germ_seq in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_mut_germ_seq PARAMETER IN nextflow.config FILE:\n${clone_mut_germ_seq}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! (clone_mut_regionDefinition in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_mut_regionDefinition PARAMETER IN nextflow.config FILE:\n${clone_mut_regionDefinition}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (clone_mut_regionDefinition =~ /^(NULL|IMGT_V|IMGT_V_BY_CODONS|IMGT_V_BY_REGIONS|IMGT_V_BY_SEGMENTS|IMGT_VDJ|IMGT_VDJ_BY_REGIONS)$/) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID clone_mut_regionDefinition PARAMETER IN nextflow.config FILE:\n${clone_mut_regionDefinition}\nMUST BE EITHER \"NULL\", \"IMGT_V\", \"IMGT_V_BY_CODONS\", \"IMGT_V_BY_REGIONS\", \"IMGT_V_BY_SEGMENTS\", \"IMGT_VDJ\", \"IMGT_VDJ_BY_REGIONS\".\n\n========\n\n"
}

// Aligments
if( ! (align_clone_nb in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID align_clone_nb PARAMETER IN nextflow.config FILE:\n${align_clone_nb}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ( ! (align_clone_nb =~/^\d+$/)) || align_clone_nb.toInteger() < 2 ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID align_clone_nb PARAMETER IN nextflow.config FILE:\n${align_clone_nb}\nMUST BE A POSITIVE INTEGER VALUE EQUAL OR GREATER TO 2\n\n========\n\n"
}
if( ! (align_soft in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID align_soft PARAMETER IN nextflow.config FILE:\n${align_soft}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (align_soft =~ /^(abalign|mafft)$/) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID align_soft PARAMETER IN nextflow.config FILE:\n${align_soft}\nMUST BE EITHER \"abalign\", OR \"mafft\".\n\n========\n\n"
}
if( ! (align_seq in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID align_seq PARAMETER IN nextflow.config FILE:\n${align_seq}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (align_seq =~ /^(query|igblast_full|trimmed|fwr1|fwr2|fwr3|fwr4|cdr1|cdr2|cdr3|junction|sequence_alignment|v_sequence_alignment|d_sequence_alignment|j_sequence_alignment|c_sequence_alignment|germline_alignment|v_germline_alignment|d_germline_alignment|j_germline_alignment|c_germline_alignment)$/) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID align_seq PARAMETER IN nextflow.config FILE:\n${align_seq}\nMUST BE EITHER \"query\", \"igblast_full\", \"trimmed\", \"fwr1\", \"fwr2\", \"fwr3\", \"fwr4\", \"cdr1\", \"cdr2\", \"cdr3\", \"junction\", \"sequence_alignment\", \"v_sequence_alignment\", \"d_sequence_alignment\", \"j_sequence_alignment\", \"c_sequence_alignment\", \"germline_alignment\", \"v_germline_alignment\", \"d_germline_alignment\", \"j_germline_alignment\", \"c_germline_alignment\".\n\n========\n\n"
}
if( ! (align_abalign_options in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID align_abalign_options PARAMETER IN nextflow.config FILE:\n${align_abalign_options}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! (align_mafft_all_options in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID align_mafft_all_options PARAMETER IN nextflow.config FILE:\n${align_mafft_all_options}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! (align_mafft_clonal_options in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID align_mafft_clonal_options PARAMETER IN nextflow.config FILE:\n${align_mafft_clonal_options}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! (meta_legend in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID meta_legend PARAMETER IN nextflow.config FILE:\n${meta_legend}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! (germ_tree_kind in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_kind PARAMETER IN nextflow.config FILE:\n${germ_tree_kind}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! (germ_tree_duplicate_seq in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_duplicate_seq PARAMETER IN nextflow.config FILE:\n${germ_tree_duplicate_seq}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (germ_tree_duplicate_seq == "TRUE" || germ_tree_duplicate_seq == "FALSE") ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_duplicate_seq PARAMETER IN nextflow.config FILE:\n${germ_tree_duplicate_seq}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
}
if( ! (germ_tree_leaf_color in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_leaf_color PARAMETER IN nextflow.config FILE:\n${germ_tree_leaf_color}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! (germ_tree_leaf_shape in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_leaf_shape PARAMETER IN nextflow.config FILE:\n${germ_tree_leaf_shape}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (germ_tree_leaf_shape =~  /^[0-9]*$/) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_leaf_shape PARAMETER IN nextflow.config FILE:\n${germ_tree_leaf_shape}\nMUST BE A POSITIVE INTEGER VALUE\n\n========\n\n"
}
if( ! (germ_tree_leaf_size in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_leaf_size PARAMETER IN nextflow.config FILE:\n${germ_tree_leaf_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (germ_tree_leaf_size =~  /^[0-9]+\.*[0-9]*$/) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_leaf_size PARAMETER IN nextflow.config FILE:\n${germ_tree_leaf_size}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
}
if( ! (germ_tree_label_size in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_size PARAMETER IN nextflow.config FILE:\n${germ_tree_label_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (germ_tree_label_size =~  /^[0-9]+\.*[0-9]*$/) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_size PARAMETER IN nextflow.config FILE:\n${germ_tree_label_size}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
}
if( ! (germ_tree_label_hjust in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_hjust PARAMETER IN nextflow.config FILE:\n${germ_tree_label_hjust}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (germ_tree_label_hjust =~  /^\-{0,1}[0-9]+\.*[0-9]*$/) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_hjust PARAMETER IN nextflow.config FILE:\n${germ_tree_label_hjust}\nMUST BE A NUMERIC VALUE\n\n========\n\n"
}
if( ! (germ_tree_label_rigth in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_rigth PARAMETER IN nextflow.config FILE:\n${germ_tree_label_rigth}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (germ_tree_label_rigth == "TRUE" || germ_tree_label_rigth == "FALSE") ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_rigth PARAMETER IN nextflow.config FILE:\n${germ_tree_label_rigth}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
}
if( ! (germ_tree_label_outside in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_outside PARAMETER IN nextflow.config FILE:\n${germ_tree_label_outside}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (germ_tree_label_outside == "TRUE" || germ_tree_label_outside == "FALSE") ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_label_outside PARAMETER IN nextflow.config FILE:\n${germ_tree_label_outside}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
}
if( ! (germ_tree_right_margin in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_right_margin PARAMETER IN nextflow.config FILE:\n${germ_tree_right_margin}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (germ_tree_right_margin =~  /^[0-9]+\.*[0-9]*$/) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_right_margin PARAMETER IN nextflow.config FILE:\n${germ_tree_right_margin}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
}
if( ! (germ_tree_legend in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_legend PARAMETER IN nextflow.config FILE:\n${germ_tree_legend}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (germ_tree_legend == "TRUE" || germ_tree_legend == "FALSE") ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID germ_tree_legend PARAMETER IN nextflow.config FILE:\n${germ_tree_legend}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
}



if( ! (donut_palette in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_palette PARAMETER IN nextflow.config FILE:\n${donut_palette}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! (donut_hole_size in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_hole_size PARAMETER IN nextflow.config FILE:\n${donut_hole_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (donut_hole_size =~  /^((1)|(0)|(0\.[0-9]*))$/) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_hole_size PARAMETER IN nextflow.config FILE:\n${donut_hole_size}\nMUST BE A POSITIVE PROPORTION VALUE\n\n========\n\n"
}
if( ! (donut_hole_text in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_hole_text PARAMETER IN nextflow.config FILE:\n${germ_tree_legend}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (donut_hole_text == "TRUE" || donut_hole_text == "FALSE") ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_hole_text PARAMETER IN nextflow.config FILE:\n${donut_hole_text}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
}
if( ! (donut_hole_text_size in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_hole_text_size PARAMETER IN nextflow.config FILE:\n${donut_hole_text_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (donut_hole_text_size =~  /^[0-9]+\.*[0-9]*$/) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_hole_text_size PARAMETER IN nextflow.config FILE:\n${donut_hole_text_size}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
}
if( ! (donut_border_color in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_border_color PARAMETER IN nextflow.config FILE:\n${donut_border_color}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! (donut_border_size in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_border_size PARAMETER IN nextflow.config FILE:\n${donut_border_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! (donut_annotation_distance in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_distance PARAMETER IN nextflow.config FILE:\n${donut_annotation_distance}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (donut_annotation_distance =~  /^[0-9]+\.*[0-9]*$/) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_distance PARAMETER IN nextflow.config FILE:\n${donut_annotation_distance}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
}
if( ! (donut_annotation_size in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_size PARAMETER IN nextflow.config FILE:\n${donut_annotation_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (donut_annotation_size =~  /^[0-9]+\.*[0-9]*$/) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_size PARAMETER IN nextflow.config FILE:\n${donut_annotation_size}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
}
if( ! (donut_annotation_force in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_force PARAMETER IN nextflow.config FILE:\n${donut_annotation_force}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (donut_annotation_force =~  /^[0-9]+\.*[0-9]*$/) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_force PARAMETER IN nextflow.config FILE:\n${donut_annotation_force}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
}
if( ! (donut_annotation_force_pull in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_force_pull PARAMETER IN nextflow.config FILE:\n${donut_annotation_force_pull}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (donut_annotation_force_pull =~  /^[0-9]+\.*[0-9]*$/) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_annotation_force_pull PARAMETER IN nextflow.config FILE:\n${donut_annotation_force_pull}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
}
if( ! (donut_legend_width in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_width PARAMETER IN nextflow.config FILE:\n${donut_legend_width}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (donut_legend_width =~  /^((1)|(0)|(0\.[0-9]*))$/) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_width PARAMETER IN nextflow.config FILE:\n${donut_legend_width}\nMUST BE A POSITIVE PROPORTION VALUE\n\n========\n\n"
}
if( ! (donut_legend_text_size in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_text_size PARAMETER IN nextflow.config FILE:\n${donut_legend_text_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (donut_legend_text_size =~  /^[0-9]+\.*[0-9]*$/) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_text_size PARAMETER IN nextflow.config FILE:\n${donut_legend_text_size}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
}
if( ! (donut_legend_box_size in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_box_size PARAMETER IN nextflow.config FILE:\n${donut_legend_box_size}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (donut_legend_box_size =~  /^[0-9]+\.*[0-9]*$/) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_box_size PARAMETER IN nextflow.config FILE:\n${donut_legend_box_size}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
}
if( ! (donut_legend_box_space in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_box_space PARAMETER IN nextflow.config FILE:\n${donut_legend_box_space}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (donut_legend_box_space =~  /^[0-9]+\.*[0-9]*$/) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_box_space PARAMETER IN nextflow.config FILE:\n${donut_legend_box_space}\nMUST BE A POSITIVE NUMERIC VALUE\n\n========\n\n"
}
if( ! (donut_legend_limit in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_limit PARAMETER IN nextflow.config FILE:\n${donut_legend_limit}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (donut_legend_width ==  "NULL") ){
    if( ! (donut_legend_width =~  /^((1)|(0)|(0\.[0-9]*))$/) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID donut_legend_limit PARAMETER IN nextflow.config FILE:\n${donut_legend_limit}\nMUST BE A POSITIVE PROPORTION VALUE IF NOT \"NULL\"\n\n========\n\n"
    }
}
if( ! (phylo_tree_model_path in String || phylo_tree_model_path in GString) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID phylo_tree_model_path PARAMETER IN nextflow.config FILE:\n${phylo_tree_model_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (file(phylo_tree_model_path).exists()) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID phylo_tree_model_path PARAMETER IN nextflow.config FILE (DOES NOT EXIST): ${phylo_tree_model_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
}
if( ! (phylo_tree_itolkey in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID phylo_tree_itolkey PARAMETER IN nextflow.config FILE:\n${phylo_tree_itolkey}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
if( ! (phylo_tree_itol_subscription in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID phylo_tree_itol_subscription PARAMETER IN nextflow.config FILE:\n${phylo_tree_itol_subscription}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (phylo_tree_itol_subscription == "TRUE" || phylo_tree_itol_subscription == "FALSE") ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID phylo_tree_itol_subscription PARAMETER IN nextflow.config FILE:\n${phylo_tree_itol_subscription}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
}
if( ! (cute_path in String || cute_path in GString) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID cute_path PARAMETER IN nextflow.config FILE:\n${cute_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}else if( ! (file(cute_path).exists()) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID cute_path PARAMETER IN nextflow.config FILE (DOES NOT EXIST): ${cute_path}\nIF POINTING TO A DISTANT SERVER, CHECK THAT IT IS MOUNTED\n\n========\n\n"
}
if( ! (igphylm_exe_path in String) ){
    error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nINVALID igphylm_exe_path PARAMETER IN nextflow.config FILE:\n${igphylm_exe_path}\nMUST BE A SINGLE CHARACTER STRING\n\n========\n\n"
}
//// end check of config file parameters
}
