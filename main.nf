nextflow.enable.dsl=2
/*
#########################################################################
##                                                                     ##
##     main.nf of repertoire profiler                                  ##
##                                                                     ##
##     Gael A. Millot                                                  ##
##     Chloe Taurel                                                    ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Computational Biology Department                                ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################
*/

//////// Processes


// see ./modules folder





    /*
// alakazam::readChangeoDb() ; dowser::formatClones() ; dowser::getTrees()
process get_germ_tree {
    label 'immcantation_10cpu'
    publishDir path: "${out_path}/RData", mode: 'copy', pattern: "{*_get_germ_tree_cloneID.RData}", overwrite: false
    cache 'true'

    input:
    path mutation_load_ch // parallelization expected
    path meta_file // just to determine if metadata have been provided (TRUE means NULL) meta_file_ch not required here
    path cute_file
    val align_clone_nb
    val germ_tree_duplicate_seq
    val igphylm_exe_path // warning: here val and not path because we do not want the igphyml file to be imported in the work dir

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
    set -o pipefail
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
"${align_clone_nb}" \
"${germ_tree_duplicate_seq}" \
"${igphylm_exe_path}" \
"${cute_file}" \
"get_germ_tree.log"
    """
}

process germ_tree_vizu {
    label 'r_ig_clustering'
    publishDir path: "${out_path}/pdf", mode: 'copy', pattern: "{germ_tree.pdf}", overwrite: false
    publishDir path: "${out_path}/pdf", mode: 'copy', pattern: "{germ_no_tree.pdf}", overwrite: false
    publishDir path: "${out_path}/figures/png", mode: 'copy', pattern: "{*.png}", overwrite: false
    publishDir path: "${out_path}/figures/svg", mode: 'copy', pattern: "{*.svg}", overwrite: false
    publishDir path: "${out_path}/RData", mode: 'copy', pattern: "{all_trees.RData}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{germ_tree_vizu.log}", overwrite: false
    cache 'true'

    input:
    path rdata_germ_tree_ch2 // no more parallelization
    val germ_tree_kind
    val align_clone_nb
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
    set -o pipefail
    germ_tree_vizu.R \
"${germ_tree_kind}" \
"${align_clone_nb}" \
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
    */


//////// End Processes


//////// Checking file


//////// end Checking file


//////// Modules

include {CheckVariables} from './conf/CheckVariables.nf'
include {reportEmptyProcess; copyLogFile} from './modules/Functions.nf'
include {Unzip} from './modules/Unzip.nf'
include {Split_fasta} from './modules/Split_fasta.nf'
include {WorkflowParam} from './modules/WorkflowParam.nf'
include {Igblast_data_check} from './modules/Igblast_data_check.nf'
include {Igblast_query} from './modules/Igblast_query.nf'
include {ParseDb_filtering} from './modules/ParseDb_filtering.nf'
include {Igblast_chain_check} from './modules/Igblast_chain_check.nf'
include {TrimTranslate} from './modules/TrimTranslate.nf'
include {DistToNearest} from './modules/DistToNearest.nf'
include {Distance_hist} from './modules/Distance_hist.nf'
include {Histogram_assembly} from './modules/Histogram_assembly.nf'
include {Add_dotted_coord} from './modules/Add_dotted_coord.nf'
include {Translate_with_IMGT_gaps} from './modules/Translate_with_imgt_gaps.nf'
include {Seq_name_replacement} from './modules/Seq_name_replacement.nf'
include {Data_assembly} from './modules/Data_assembly.nf'
include {Metadata_check} from './modules/Metadata_check.nf'
include {Repertoire} from './modules/Repertoire.nf'
include {Clone_assignment} from './modules/Clone_assignment.nf'
include {Split_by_clones} from './modules/Split_by_clones.nf'
include {Closest_germline} from './modules/Closest_germline.nf'
include {Igblast_germline_coords} from './modules/Igblast_germline_coords.nf'
include {AddGermlineSequences} from './modules/AddGermlineSequences.nf'
include {Mutation_load_germ_genes} from './modules/Mutation_load_germ_genes.nf'
include {TranslateGermline} from './modules/TranslateGermline.nf'
include {Clone_id_count} from './modules/Clone_id_count.nf'
include {Donut} from './modules/Donut.nf'
include {Donut_assembly} from './modules/Donut_assembly.nf'
include {Tsv2fasta} from './modules/Tsv2fasta.nf'
include {Gff_imgt} from './modules/Gff_imgt.nf'
include {Mafft_align} from './modules/Mafft_align.nf'
include {Abalign_align_aa} from './modules/Abalign_align_aa.nf'
include {Abalign_rename} from './modules/Abalign_rename.nf'
include {Abalign_align_nuc} from './modules/Abalign_align_nuc.nf'
include {Tree_nuc} from './modules/Tree_nuc.nf'
include {Tree_aa} from './modules/Tree_aa.nf'
include {Meta2Itol} from './modules/Meta2Itol.nf'
include {Itol} from './modules/Itol.nf'
include {Print_warnings} from './modules/Print_warnings.nf'
include {Print_report} from './modules/Print_report.nf'
include {Backup} from './modules/Backup.nf'



//include {CopyLogFile as CopyLogFile_Closest_germline} from './modules/CopyLogFile.nf'
include {Gff as GffNuc} from './modules/Gff.nf'
include {Gff as GffAa} from './modules/Gff.nf'
include {PrintAlignment as PrintAlignmentIMGTnuc} from './modules/Print_alignment.nf'
include {PrintAlignment as PrintAlignmentIMGTaa} from './modules/Print_alignment.nf'
include {PrintAlignment as PrintAlignmentNuc} from './modules/Print_alignment.nf'
include {PrintAlignment as PrintAlignmentAa} from './modules/Print_alignment.nf'

//////// end Modules


//////// Workflow


workflow {

    print("\n\nINITIATION TIME: ${workflow.start}")

    //////// Options of nextflow run

    // --modules (it is just for the process WorkflowParam)
    params.modules = "" // if --module is used, this default value will be overridden
    // end --modules (it is just for the process WorkflowParam)

    //////// end Options of nextflow run


    //////// Variables

    modules = params.modules // remove the dot -> can be used in bash scripts
    config_file = workflow.configFiles[0] // better to use this than config_file = file("${projectDir}/nextflow.config") because the latter is not good if -c option of nextflow run is used
    log_file = file("${launchDir}/.nextflow.log")

    //////// end Variables


    //////// inititation

    print("\n\nRESULT DIRECTORY: ${out_path}")
    if("${system_exec}" == "slurm"){
        print("    queue: ${slurm_queue}")
        print("    qos: ${slurm_qos}")
    }
    if("${system_exec}" != "local"){
        print("    add_options: ${add_options}")
    }
    warning_ch = Channel.empty() // to collect all the warnings
    warn = "\n\nWARNING:\nCURRENTLY VALIDATED FOR 'mouse' OR 'human' SPECIES. POTENTIAL ERRORS WITH OTHER SPECIES IF VDJ STRUCTURATION OF CHAINS IS DIFFERENT."
    print(warn)
    warning_ch = warning_ch.mix(Channel.value(warn))
    warn = "\n\nWARNING:\nTO MAKE THE REPERTOIRES AND DONUTS, THE SCRIPT CURRENTLY TAKES THE FIRST ANNOTATION OF THE IMGT ANNOTATION IF SEVERAL ARE PRESENTS IN THE v_call, j_call OR c_call COLUMN OF THE wanted_seq.tsv FILE."
    print(warn)
    warning_ch = warning_ch.mix(Channel.value(warn))
    print("\n\n")

    //////// end inititation


    //////// Variable modification

    // CONSTRUCTION OF THE igblast REFERENCE FILES PATHS

    igblast_v_ref_files = ""
    igblast_d_ref_files = ""
    igblast_j_ref_files = ""
    igblast_constant_ref_files = ""

    def var_files = []
    def v_files = []
    def d_files = []
    def j_files = []
    def const_files = []

    if (igblast_loci == "ig" && igblast_B_heavy_chain == "TRUE") {
        v_files += [
            "imgt_${igblast_organism}_IGHV.fasta"
        ]
        d_files += [
            "imgt_${igblast_organism}_IGHD.fasta"
        ]
        j_files += [
            "imgt_${igblast_organism}_IGHJ.fasta"
        ]
        const_files += "imgt_${igblast_organism}_IGHC.fasta"
    }
    if (igblast_loci == "ig" && igblast_B_lambda_chain == "TRUE") {
        v_files += [
            "imgt_${igblast_organism}_IGLV.fasta"
        ]
        j_files += [
            "imgt_${igblast_organism}_IGLJ.fasta"
        ]
        const_files += "imgt_${igblast_organism}_IGLC.fasta"
    }
    if (igblast_loci == "ig" && igblast_B_kappa_chain == "TRUE") {
        v_files += [
            "imgt_${igblast_organism}_IGKV.fasta"
        ]
        j_files += [
            "imgt_${igblast_organism}_IGKJ.fasta"
        ]
        const_files += "imgt_${igblast_organism}_IGKC.fasta"
    }
    if (igblast_loci == "tr" && igblast_T_alpha_chain == "TRUE") {
        v_files += [
            "imgt_${igblast_organism}_TRAV.fasta"
        ]
        j_files += [
            "imgt_${igblast_organism}_TRAJ.fasta"
        ]
        const_files += "imgt_${igblast_organism}_TRAC.fasta"
    }
    if (igblast_loci == "tr" && igblast_T_beta_chain == "TRUE") {
        v_files += [
            "imgt_${igblast_organism}_TRBV.fasta"
        ]
        d_files += [
            "imgt_${igblast_organism}_TRBD.fasta"
        ]
        j_files += [
            "imgt_${igblast_organism}_TRBJ.fasta"
        ]
        const_files += "imgt_${igblast_organism}_TRBC.fasta"
    }
    if (igblast_loci == "tr" && igblast_T_gamma_chain == "TRUE") {
        v_files += [
            "imgt_${igblast_organism}_TRGV.fasta"
        ]
        j_files += [
            "imgt_${igblast_organism}_TRGJ.fasta"
        ]
        const_files += "imgt_${igblast_organism}_TRGC.fasta"
    }
    if (igblast_loci == "tr" && igblast_T_delta_chain == "TRUE") {
        v_files += [
            "imgt_${igblast_organism}_TRDV.fasta"
        ]
        d_files += [
            "imgt_${igblast_organism}_TRDD.fasta"
        ]
        j_files += [
            "imgt_${igblast_organism}_TRDJ.fasta"
        ]
        const_files += "imgt_${igblast_organism}_TRDC.fasta"
    }

    // names of the ref files without path
    igblast_v_ref_files = v_files.join(' ').toString()
    igblast_d_ref_files = d_files.join(' ').toString()
    if(igblast_d_ref_files == ""){igblast_d_ref_files = "NULL"} // NULL file created if no D (K and L for instance)
    igblast_j_ref_files = j_files.join(' ').toString()
    igblast_constant_ref_files = const_files.join(' ').toString()
    igblast_variable_ref_files = "${igblast_v_ref_files} ${igblast_d_ref_files} ${igblast_j_ref_files}"
    // end names of the ref files without path

    //////// end Variable modification


    //////// Channels

    // fs_ch define below because can be a .zip file
    // warning_ch = Channel.empty() // already set above
    // for print_report
    nb_productive = Channel.empty()
    nb_unproductive = Channel.empty()
    nb_wanted = Channel.empty()
    nb_unwanted = Channel.empty()
    nb_dist_ignored = Channel.empty()
    nb_clone_assignment = Channel.empty()
    nb_clone_unassignment = Channel.empty()
    nb_clone_germline = Channel.empty()
    nb_clone_ungermline = Channel.empty()
    nb_clone_tot = Channel.empty()
    nb_unclone_tot = Channel.empty()
    distance_hist_ch  = Channel.empty() // Distance_hist.out.distance_hist_ch, 
    donuts_png_ch = Channel.empty() // Donut.out.donuts_png.collect(), 
    repertoire_png_ch = Channel.empty() // Repertoire.out.repertoire_png_ch.collect(), 
    repertoire_constant_ch = Channel.empty()
    repertoire_vj_ch = Channel.empty()
    empty_distance_hist = file('empty_distance_hist') // to deal with empty path channel
    empty_donuts_png = file('empty_donuts_png') // to deal with empty path channel
    empty_repertoire_png = file('empty_repertoire_png') // to deal with empty path channel
    // end for print_report


    //////// end Channels


    //////// files import

    meta_file = file(meta_path) // in variable because a single file. If "NULL", will create a empty file, present in work folders, but that cannot be correctly linked. Thus, if the file has to be redirected into a channel inside a process, it will not work. Thus, in the first process using meta_file, I hard copy the NULL file if required (see below)
    cute_file = file(cute_path) // in variable because a single file
    phylo_tree_model_file  = file(phylo_tree_model_path)
    template_rmd = file(template_rmd_path)
    alignments_viz_rmd = file(alignments_viz_path)
    alignments_viz_html = file(alignments_viz_html_path)

    //////// end files import


    //////// Main

    CheckVariables()

    if(sample_path =~ /.*\.zip$/){
        Unzip( // warning: it is a process defined above
            Channel.fromPath(sample_path),
            sample_path
        ) 
        dir_ch = Unzip.out.unzip_ch.flatten()
    }else{
        dir_ch = Channel.fromPath("${sample_path}", checkIfExists: false) // in channel because many files 
    }


    // is the path a dir or a single file ?
    dir_ch.branch {
            dir: it.isDirectory()
            file: true
        }.set { branched }
    // Handle directories: list contents
    fs_ch_from_dir = branched.dir.flatMap { it.listFiles() } // is it is a dir, then recover all the files
    // Handle files: pass through
    fs_ch_from_file = branched.file
    // Merge back
    fs_ch = fs_ch_from_dir.mix(fs_ch_from_file)
    // end is the path a dir or a single file ?

    fs_ch.toList().branch {
            single: it.size() == 1
                return it[0]
            multiple: true
                return it
        }.set { branched }
    Split_fasta(branched.single)
    fs_ch2 = Split_fasta.out.split_fasta_ch.mix(branched.multiple.flatten()).flatten()

    nb_input = fs_ch2.count()

    WorkflowParam(
        modules
    )

    file("${out_path}/tsv").mkdirs()
    file("${out_path}/pdf").mkdirs()

    Igblast_data_check(
        igblast_organism, 
        igblast_v_ref_files, 
        igblast_d_ref_files, 
        igblast_j_ref_files, 
        igblast_constant_ref_files
    )

    Igblast_query( // module igblast_query.nf
        fs_ch2, 
        igblast_variable_ref_files, 
        igblast_organism, 
        igblast_loci
    )
    igblast_ch1 = Igblast_query.out.db_pass_ch.collectFile(name: "igblast_seq.tsv", skip: 1, keepHeader: true) // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files
    unigblast_ch1 = Igblast_query.out.db_unpass_ch.collectFile(name: "failed_igblast_seq.tsv", skip: 1, keepHeader: true) // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files
    igblast_ch1.subscribe{it -> it.copyTo("${out_path}/tsv")} // nothing if no file
    unigblast_ch1.subscribe{it -> it.copyTo("${out_path}/tsv")}
    nb_igblast = igblast_ch1.countLines() - 1 // -1 for the header
    nb_unigblast = unigblast_ch1.countLines() - 1 // -1 for the header
    check_igblast = igblast_ch1.ifEmpty {'NO_FILE'}   // marker \ue202turn0search2
    check_unigblast = unigblast_ch1.ifEmpty {'NO_FILE'}   // marker \ue202turn0search2
    check_igblast.combine(check_unigblast).subscribe {x,unx -> 
        if(x != 'NO_FILE'){
            n = x.countLines() - 1   // here `x` is a Path, so it exists
            if(n == -1){
                throw new IllegalStateException("\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nIGBLAST FILE IS EMPTY, WHILE IT SHOULD HAVE AT LEAST A HEADER, FOLLOWING THE Igblast_query PROCESS.\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n")
            }
            if(n == 0){
                throw new IllegalStateException("\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\n0 ANNOTATION SUCCEEDED BY THE Igblast_query PROCESS.\n\nCHECK THAT THE\nigblast_organism\nigblast_loci\nigblast_B_heavy_chain\nigblast_B_lambda_chain\nigblast_B_kappa_chain\nigblast_T_alpha_chain\nigblast_T_beta_chain\nigblast_T_gamma_chain\nigblast_T_delta_chain\nARE CORRECTLY SET IN THE nextflow.config FILE\nAND CHECK THE SUBMITTED SEQUENCES (VERY BAD QUALITY).\n\n========\n\n")
            }
        }
        if(unx != 'NO_FILE'){
            n = unx.countLines() - 1   // here `x` is a Path, so it exists
            if(n == -1){
                throw new IllegalStateException("\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nFAILED IGBLAST FILE IS EMPTY, WHILE IT SHOULD HAVE AT LEAST A HEADER, FOLLOWING THE Igblast_query PROCESS.\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n")
            }
        }
    }
    fs_ch2.count().combine(nb_igblast).combine(nb_unigblast).subscribe{n,n1,n2 -> 
        if(n != -1 && n1 != -1 && n2 != -1){
            if(n != n1 + n2){
                throw new IllegalStateException("\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nTHE NUMBER OF LINES IN THE igblast_aligned_seq.tsv (${n1}) AND igblast_unaligned_seq.tsv (${n2}) IS NOT EQUAL TO THE NUMBER OF SUBMITTED FASTA FILES (${n})\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n")
            }
        }
    }
    copyLogFile('igblast_query.log', Igblast_query.out.log_ch, out_path)


    // when: nb_igblast > 0


        ParseDb_filtering(
            Igblast_query.out.db_pass_ch.map{file -> nlines = file.text.readLines().size() ; tuple(file, nlines)}
        )
        productive_ch1 = ParseDb_filtering.out.productive_ch.collectFile(name: "productive_seq.tsv", skip: 1, keepHeader: true) // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files
        unproductive_ch1 = ParseDb_filtering.out.unproductive_ch.collectFile(name: "unproductive_seq.tsv", skip: 1, keepHeader: true) // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files
        productive_ch1.subscribe{it -> it.copyTo("${out_path}/tsv")} // nothing if no file
        unproductive_ch1.subscribe{it -> it.copyTo("${out_path}/tsv")}
        nb_productive = nb_productive.mix(productive_ch1.countLines() - 1) // -1 for the header
        nb_unproductive = nb_unproductive.mix(unproductive_ch1.countLines() - 1) // -1 for the header
        check_productive = productive_ch1.ifEmpty {'NO_FILE'}   // marker \ue202turn0search2
        check_unproductive = unproductive_ch1.ifEmpty {'NO_FILE'}   // marker \ue202turn0search2
        warning_ch = warning_ch.mix(check_productive.combine(check_unproductive).map{x,unx -> 
            if(x != 'NO_FILE'){
                n = x.countLines() - 1   // here `x` is a Path, so it exists
                if(n == -1){
                    throw new IllegalStateException("\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nPRODUCTIVE FILE IS EMPTY, WHILE IT SHOULD HAVE AT LEAST A HEADER, FOLLOWING THE ParseDb_filtering PROCESS.\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n")
                }
                if(n == 0){
                    warn = "\n\nWARNING:\n0 PRODUCTIVE SEQUENCE FOLLOWING THE ParseDb_filtering PROCESS\nWORFLOW ENDED.\nWAIT FOR THE END OF THE Print_report PROCESS AND SEE THE PARTIAL RESULTS IN:\n${out_path}.\n\n"
                    print(warn)
                    return warn // accumulate
                }
            }else{
                return null
            }
            if(unx != 'NO_FILE'){
                n = unx.countLines() - 1   // here `x` is a Path, so it exists
                if(n == -1){
                    throw new IllegalStateException("\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nUNPRODUCTIVE FILE IS EMPTY, WHILE IT SHOULD HAVE AT LEAST A HEADER, FOLLOWING THE ParseDb_filtering PROCESS.\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n")
                }else{
                    return null
                }
            }else{
                return null
            }
        })
        nb_igblast.combine(nb_productive).combine(nb_unproductive).subscribe{n,n1,n2 -> 
            if(n != -1 && n1 != -1 && n2 != -1){
                if(n != n1 + n2){
                    throw new IllegalStateException("\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nTHE NUMBER OF LINES IN THE productive_seq.tsv (${n1}) AND unproductive_seq.tsv (${n2}) IS NOT EQUAL TO THE NUMBER OF LINES IN igblast_seq.tsv (${n})\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n")
                }
            }
        }
        copyLogFile('parseDb_filtering.log', ParseDb_filtering.out.parseDb_filtering_log_ch, out_path)

        // when: nb_productive > 0

            Igblast_chain_check( // module Igblast_chain_check.nf
                ParseDb_filtering.out.productive_ch.map{file -> nlines = file.text.readLines().size() ; tuple(file, nlines)}, 
                Igblast_data_check.out.allele_names_tsv_all_ch, 
                igblast_v_ref_files, 
                igblast_d_ref_files, 
                igblast_j_ref_files, 
                igblast_constant_ref_files, 
                cute_file
            )
            check_ch1 = Igblast_chain_check.out.checked_tsv_ch.collectFile(name: "wanted_seq.tsv", skip: 1, keepHeader: true) // not printed here. Later // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files
            uncheck_ch1 = Igblast_chain_check.out.not_checked_tsv_ch.collectFile(name: "unwanted_seq.tsv", skip: 1, keepHeader: true) // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files
            // check_ch1.subscribe{it -> it.copyTo("${out_path}/tsv")} // nothing if no file // not now
            uncheck_ch1.subscribe{it -> it.copyTo("${out_path}/tsv")}
            nb_wanted = nb_wanted.mix(check_ch1.countLines() - 1) // -1 for the header
            nb_unwanted = nb_unwanted.mix(uncheck_ch1.countLines() - 1) // -1 for the header
            check_check = check_ch1.ifEmpty {'NO_FILE'}   // marker \ue202turn0search2
            check_uncheck = uncheck_ch1.ifEmpty {'NO_FILE'}   // marker \ue202turn0search2
            warning_ch = warning_ch.mix(check_check.combine(check_uncheck).map{x,unx -> 
                if(x != 'NO_FILE'){
                    n = x.countLines() - 1   // here `x` is a Path, so it exists
                    if(n == -1){
                        throw new IllegalStateException("\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nWANTED FILE IS EMPTY, WHILE IT SHOULD HAVE AT LEAST A HEADER, FOLLOWING THE Igblast_chain_check PROCESS.\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n")
                    }
                    if(n == 0){
                        warn = "\n\nWARNING:\n0 WANTED SEQUENCE FOLLOWING THE Igblast_chain_check PROCESS.\nCHECK THAT THE\nigblast_organism\nigblast_loci\nigblast_B_heavy_chain\nigblast_B_lambda_chain\nigblast_B_kappa_chain\nigblast_T_alpha_chain\nigblast_T_beta_chain\nigblast_T_gamma_chain\nigblast_T_delta_chain ARE CORRECTLY SET IN THE nextflow.config FILE\nWORFLOW ENDED.\nWAIT FOR THE END OF THE Print_report PROCESS AND SEE THE PARTIAL RESULTS IN:\n${out_path}.\n\n"
                        print(warn)
                        return warn // accumulate
                    }
                }else{
                    return null
                }
                if(unx != 'NO_FILE'){
                    n = unx.countLines() - 1   // here `x` is a Path, so it exists
                    if(n == -1){
                        throw new IllegalStateException("\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nUNWANTED FILE IS EMPTY, WHILE IT SHOULD HAVE AT LEAST A HEADER, FOLLOWING THE igblast_B_heavy_chain PROCESS.\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n")
                    }else{
                        return null
                    }
                }else{
                    return null
                }
            })
            nb_productive.combine(nb_wanted).combine(nb_unwanted).subscribe{n,n1,n2 -> 
                if(n != -1 && n1 != -1 && n2 != -1){
                    if(n != n1 + n2){
                        throw new IllegalStateException("\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nTHE NUMBER OF LINES IN THE wanted_seq.tsv (${n1}) AND unwanted_seq.tsv (${n2}) IS NOT EQUAL TO THE NUMBER OF LINES IN productive_seq.tsv (${n})\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n")
                    }
                }
            }
            copyLogFile('igblast_chain_check.log', Igblast_chain_check.out.igblast_chain_check_log, out_path)


            // when: nb_wanted > 0

                TrimTranslate(
                    Igblast_chain_check.out.checked_tsv_ch.map{file -> nlines = file.text.readLines().size() ; tuple(file, nlines)}
                )
                //TrimTranslate.out.trimtranslate_ch.count().subscribe {n -> if ( n == 0 ){error "\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nO SEQUENCE RETURNED FOLLOWING THE TrimTranslate PROCESS\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n"}} // n is number of items in channel because of .count()
                trimtranslate_ch2 = TrimTranslate.out.trimtranslate_ch.collectFile(name: "trimtranslate.tsv", skip: 1, keepHeader: true) // wanted file with column sequence_alignment_aa added  // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files
                copyLogFile('trimtranslate.log', TrimTranslate.out.trimtranslate_log_ch, out_path)


                DistToNearest(
                    trimtranslate_ch2,
                    clone_model,
                    clone_normalize
                )
                nb_dist_ignored = nb_dist_ignored.mix(DistToNearest.out.dist_ignored_ch.countLines() - 1) // Minus 1 because 1st line = column names // either .tsv is empty (header alone) or has a header to remove to count the lines


                Distance_hist(
                    DistToNearest.out.distToNearest_ch,
                    cute_file, 
                    clone_model,
                    clone_normalize,
                    clone_distance
                )
                distance_hist_ch = distance_hist_ch.mix(Distance_hist.out.distance_hist_ch)


                Histogram_assembly(
                    Distance_hist.out.histogram_pdf_ch.collect()
                )


                Add_dotted_coord( // module
                    TrimTranslate.out.trimtranslate_ch
                )
                copyLogFile('add_dotted_coord.log', Add_dotted_coord.out.add_dotted_coord_log_ch, out_path)


                Translate_with_IMGT_gaps( // module
                    Add_dotted_coord.out.add_dotted_coord_ch
                )
                copyLogFile('Translate_with_IMGT_gaps.log', Translate_with_IMGT_gaps.out.translate_with_IMGT_gaps_ch_log, out_path)


                Seq_name_replacement(
                    Translate_with_IMGT_gaps.out.add_aa_imgt_ch,
                    meta_file,
                    meta_seq_names, 
                    meta_name_replacement,
                    meta_legend
                )
                seq_name_replacement_ch2 = Seq_name_replacement.out.seq_name_replacement_ch.collectFile(name: "replacement.tsv", skip: 1, keepHeader: true) // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files
                // seq_name_replacement_ch2.subscribe{it -> it.copyTo("${out_path}")}
                // tuple_seq_name_replacement = new Tuple("all", seq_name_replacement_ch2) # warning: this is not a channel but a variable now
                copyLogFile('seq_name_replacement.log', Seq_name_replacement.out.seq_name_replacement_log_ch, out_path)


                Data_assembly(
                    seq_name_replacement_ch2, 
                    DistToNearest.out.distToNearest_ch
                )
                warning_ch = warning_ch.mix(Data_assembly.out.data_assembly_warn_ch.filter{file(it).exists()}.map{file -> file.text}) //  file.text = contenu du fichier) // accumulate


                Metadata_check(
                    Data_assembly.out.wanted_ch,
                    meta_file, 
                    meta_seq_names, 
                    meta_name_replacement,
                    meta_legend
                )
                warning_ch = warning_ch.mix(Metadata_check.out.metadata_check_warn_ch.filter{file(it).exists()}.map{file -> file.text}) //  file.text = contenu du fichier) // accumulate


                Repertoire(
                    Data_assembly.out.wanted_ch,
                    Igblast_data_check.out.allele_names_tsv_all_ch,
                    cute_file
                )
                repertoire_png_ch = repertoire_png_ch.mix(Repertoire.out.repertoire_png_ch.collect())
                repertoire_constant_ch = repertoire_constant_ch.mix(Repertoire.out.repertoire_png_ch.flatten().filter {file -> 
                    file.name =~ /^.*(IG|TR).C_.*gene_non-zero\.png$/
                })
                repertoire_vj_ch =repertoire_vj_ch.mix( Repertoire.out.repertoire_png_ch.flatten().filter {file ->
                    file.name =~ /.*\/?(rep_gene_(IG|TR).V_.*non-zero\.png)$/
                })


                Clone_assignment(
                    Data_assembly.out.wanted_ch, 
                    clone_model,
                    clone_normalize,
                    clone_distance,
                    clone_strategy,
                    meta_file,
                    meta_legend
                )
                assign_ch1 = Clone_assignment.out.clone_ch // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files
                unassign_ch1 = Clone_assignment.out.failed_clone_ch // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files
                unassign_ch1.collectFile(name: "failed_clone_assigned_seq.tsv", skip: 1, keepHeader: true).subscribe{it -> it.copyTo("${out_path}/tsv")}
                nb_clone_assignment = nb_clone_assignment.mix(assign_ch1.countLines() - 1) // -1 for the header
                nb_clone_unassignment = nb_clone_unassignment.mix(unassign_ch1.countLines() - 1) // -1 for the header
                check_assign = assign_ch1.ifEmpty {'NO_FILE'}   // marker \ue202turn0search2
                check_unassign = unassign_ch1.ifEmpty {'NO_FILE'}   // marker \ue202turn0search2
                warning_ch = warning_ch.mix(check_assign.combine(check_unassign).map{x,unx -> 
                    if(x != 'NO_FILE'){
                        n = x.countLines() - 1   // here `x` is a Path, so it exists
                        if(n == -1){
                            throw new IllegalStateException("\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nASSIGNMENT FILE IS EMPTY, WHILE IT SHOULD HAVE AT LEAST A HEADER, FOLLOWING THE Clone_assignment PROCESS.\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n")
                        }
                        if(n == 0){
                            warn = "\n\nWARNING:\n0 ASSIGNED SEQUENCE FOLLOWING THE Clone_assignment PROCESS.\nGERMLINE SEQUENCE & TREE PARTS OF THE WORFLOW ENDED.\nWAIT FOR THE END OF THE Print_report PROCESS AND SEE THE PARTIAL RESULTS IN:\n${out_path}.\n\n"
                            print(warn)
                            return warn // accumulate
                        }
                    }else{
                        return null
                    }
                    if(unx != 'NO_FILE'){
                        n = unx.countLines() - 1   // here `x` is a Path, so it exists
                        if(n == -1){
                            throw new IllegalStateException("\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nUNASSIGNMENT FILE IS EMPTY, WHILE IT SHOULD HAVE AT LEAST A HEADER, FOLLOWING THE Clone_assignment PROCESS.\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n")
                        }else{
                            return null
                        }
                    }else{
                        return null
                    }
                })
                nb_wanted.combine(nb_clone_assignment).combine(nb_clone_unassignment).subscribe{n,n1,n2 -> 
                    if(n != -1 && n1 != -1 && n2 != -1){
                        if(n != n1 + n2){
                            throw new IllegalStateException("\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nTHE NUMBER OF LINES IN THE ASSIGNMENT FILE (${n1}) AND failed_clone_assigned_seq.tsv (${n2}) IS NOT EQUAL TO THE NUMBER OF LINES IN wanted_seq.tsv (${n})\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n")
                        }
                    }
                }


                // when: nb_clone_assignment > 0

                    Split_by_clones( // split by clone ID, not by single sequence
                        Clone_assignment.out.clone_ch.map{file -> nlines = file.text.readLines().size() ; tuple(file, nlines)}
                    )


                    Closest_germline(
                        Split_by_clones.out.clone_split_ch.flatten(), // flatten split the list into several objects (required for parallelization)
                        igblast_organism, 
                        igblast_variable_ref_files, 
                        clone_germline_kind,
                        meta_file,
                        meta_legend
                    )
                    closest_ch1 = Closest_germline.out.closest_ch.map { it[0] }.collectFile(name: "closest_seq.tsv", skip: 1, keepHeader: true) // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files // it[0] to take only the tsv, not the fasta
                    unclosest_ch1 = Closest_germline.out.failed_clonal_germline_ch.collectFile(name: "failed_clonal_germline_seq.tsv", skip: 1, keepHeader: true) // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files
                    closest_ch1.subscribe{it -> it.copyTo("${out_path}/tsv")} // nothing if no file // not now
                    unclosest_ch1.subscribe{it -> it.copyTo("${out_path}/tsv")}
                    nb_clone_germline = nb_clone_germline.mix(closest_ch1.countLines() - 1) // -1 for the header
                    nb_clone_ungermline = nb_clone_ungermline.mix(unclosest_ch1.countLines() - 1) // -1 for the header
                    check_closest = closest_ch1.ifEmpty{'NO_FILE'}   // marker \ue202turn0search2
                    check_unclosest = unclosest_ch1.ifEmpty{'NO_FILE'}   // marker \ue202turn0search2
                    warning_ch = warning_ch.mix(check_closest.combine(check_unclosest).map{x,unx -> 
                        if(x != 'NO_FILE'){
                            n = x.countLines() - 1   // here `x` is a Path, so it exists
                            if(n == -1){
                                throw new IllegalStateException("\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nCLOSEST GERMLINE FILE IS EMPTY, WHILE IT SHOULD HAVE AT LEAST A HEADER, FOLLOWING THE Closest_germline PROCESS.\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n")
                            }
                            if(n == 0){
                                warn = "\n\nWARNING: 0 CLOSEST GERMLINE SEQUENCE FOLLOWING THE Closest_germline PROCESS.\nGERMLINE SEQUENCE & TREE PARTS OF THE WORFLOW ENDED.\nWAIT FOR THE END OF THE Print_report PROCESS AND SEE THE PARTIAL RESULTS IN:\n${out_path}.\n\n"
                                print(warn)
                                return warn // accumulate
                            }
                        }else{
                            return null
                        }
                        if(unx != 'NO_FILE'){
                            n = unx.countLines() - 1   // here `x` is a Path, so it exists
                            if(n == -1){
                                throw new IllegalStateException("\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nUNCLOSEST FILE IS EMPTY, WHILE IT SHOULD HAVE AT LEAST A HEADER, FOLLOWING THE Closest_germline PROCESS.\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n")
                            }else{
                                return null
                            }
                        }else{
                            return null
                        }
                    })
                    nb_clone_assignment.combine(nb_clone_germline).combine(nb_clone_ungermline).subscribe{n,n1,n2 -> 
                        if(n != -1 && n1 != -1 && n2 != -1){
                            if(n != n1 + n2){
                                throw new IllegalStateException("\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nTHE NUMBER OF LINES IN THE CLOSEST GERMLINE FILE (${n1}) AND failed_clonal_germline_seq.tsv (${n2}) IS NOT EQUAL TO THE NUMBER OF LINES IN THE CLONE ASSIGNMENT FILE (${n})\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n")
                            }
                        }
                    }
                    copyLogFile('Closest_germline.log', Closest_germline.out.closest_log_ch, out_path)


                    // when: nb_clone_germline > 0

                        Igblast_germline_coords( // module igblast_germline.nf
                            Closest_germline.out.closest_ch.map{tsv, fasta -> tsv_nlines = tsv.text.readLines().size() ; tuple(tsv, fasta, tsv_nlines)}, // no ifEmpty{error} because can be empty for some of the parallelized processes
                            igblast_organism, 
                            igblast_loci
                        )
                        copyLogFile('Igblast_germline_coords_report.log', Igblast_germline_coords.out.log_ch, out_path)


                        AddGermlineSequences(
                            Igblast_germline_coords.out.germline_coords_ch,
                            igblast_organism, 
                            igblast_variable_ref_files
                        )
                        copyLogFile('AddGermlineSequences.log', AddGermlineSequences.out.add_germ_log_ch, out_path)



                        TranslateGermline(
                            AddGermlineSequences.out.add_germ_ch,
                            cute_file
                        )
                        TranslateGermline.out.translate_problem_ch.collectFile(name: "clonal_germline_translation_pb.tsv", skip: 1, keepHeader: true).subscribe{it -> it.copyTo("${out_path}/tsv")}
                        TranslateGermline.out.translate_germ_log_ch.collectFile(name: "translateGermline.log").subscribe{it -> it.copyTo("${out_path}/reports")}
                        warning_ch = warning_ch.mix(TranslateGermline.out.translate_germ_warn_ch.filter{file(it).exists()}.map{file -> file.text}) //  file.text = contenu du fichier) // accumulate



                        Mutation_load_germ_genes(
                            TranslateGermline.out.translate_germ_ch,
                            meta_file, 
                            meta_legend,
                            clone_mut_obs_seq,
                            clone_mut_germ_seq,
                            clone_mut_regionDefinition
                        )
                        clone_assigned_seq_ch = Mutation_load_germ_genes.out.mutation_load_ch.collectFile(name: "clone_assigned_seq.tsv", skip: 1, keepHeader: true) // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files
                        clone_assigned_seq_ch.subscribe{it -> it.copyTo("${out_path}/tsv")}
                        nb_clone_tot = nb_clone_tot.mix(clone_assigned_seq_ch.countLines() - 1) // -1 for the header
                        nb_unclone_tot = nb_unclone_tot.mix(nb_clone_unassignment.combine(nb_clone_ungermline).map{n1,n2 -> n1 + n2}) // -1 for the header
                        check_clone_tot = clone_assigned_seq_ch.ifEmpty {'NO_FILE'}   // marker \ue202turn0search2
                        warning_ch = warning_ch.mix(check_clone_tot.map {x -> 
                            if(x != 'NO_FILE'){
                                n = x.countLines() - 1   // here `x` is a Path, so it exists
                                if(n == -1){
                                    throw new IllegalStateException("\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nclone_assigned_seq.tsv FILE IS EMPTY, WHILE IT SHOULD HAVE AT LEAST A HEADER, FOLLOWING THE Mutation_load_germ_genes PROCESS.\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n")
                                }
                                if(n == 0){
                                    warn = "\n\nWARNING:\n0 CLONAL SEQUENCE FOLLOWING THE Mutation_load_germ_genes PROCESS.\nGERMLINE SEQUENCE & TREE PARTS OF THE WORFLOW ENDED.\nWAIT FOR THE END OF THE Print_report PROCESS AND SEE THE PARTIAL RESULTS IN:\n${out_path}.\n\n"
                                    print(warn)
                                    return warn // accumulate
                                }
                            }else{
                                return null
                            }
                        })
                        nb_wanted.combine(nb_clone_tot).combine(nb_unclone_tot).subscribe{n,n1,n2 -> 
                            if(n != -1 && n1 != -1 && n2 != -1){
                                if(n != n1 + n2){
                                    throw new IllegalStateException("\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nTHE NUMBER OF LINES IN THE clone_assigned_seq.tsv FILE (${n1}) AND failed_clone_assigned_seq.tsv + failed_clonal_germline_seq.tsv (${n2}) IS NOT EQUAL TO THE NUMBER OF LINES IN THE wanted_seq.tsv FILE (${n})\n\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\n\n========\n\n")
                                }
                            }
                        }
                        clone_assigned_seq_filtered_ch = Mutation_load_germ_genes.out.mutation_load_ch.filter{file -> file.countLines() > align_clone_nb.toInteger()} // Only keep clonal groups that have a number of sequences superior to align_clone_nb (variable defined in nextflow.config) 
                        copyLogFile('mutation_load_germ_genes.log', Mutation_load_germ_genes.out.mutation_load_log_ch, out_path)



                        Clone_id_count(
                            clone_assigned_seq_ch
                        )

                    // end when: nb_clone_germline > 0

                    /*

                    get_germ_tree(
                        Mutation_load_germ_genes.out.mutation_load_ch,
                        meta_file, // first() because get_germ_tree process is a parallele one and because meta_file is single
                        cute_file, 
                        align_clone_nb,
                        germ_tree_duplicate_seq,
                        igphylm_exe_path
                    )
                    
                    get_germ_tree.out.rdata_germ_tree_ch.count().subscribe {n -> if ( n == 0 ){
                        print("\n\nWARNING:\nEMPTY OUTPUT FOLLOWING THE get_germ_tree PROCESS -> NO TREE RETURNED\n\n")
                    }}
                    rdata_germ_tree_ch2 = get_germ_tree.out.rdata_germ_tree_ch.collect()
                    //rdata_germ_tree_ch2.view()

                    get_germ_tree.out.no_germ_tree_ch.count().subscribe {n -> if ( n == 0 ){
                        print("\n\nWARNING:\nALL SEQUENCES IN TREES FOLLOWING THE get_germ_tree PROCESS -> EMPTY germ_tree_dismissed_seq.tsv FILE RETURNED\n\n")
                    }}
                    no_germ_tree_ch2 = get_germ_tree.out.no_germ_tree_ch.collectFile(name: "germ_tree_dismissed_seq.tsv", skip: 1, keepHeader: true)
                    no_germ_tree_ch2.subscribe{it -> it.copyTo("${out_path}/tsv")}

                    get_germ_tree.out.germ_tree_ch.count().subscribe {n -> if ( n == 0 ){
                        print("\n\nWARNING:\nNO SEQUENCES IN TREES FOLLOWING THE get_germ_tree PROCESS -> EMPTY seq_for_germ_tree.tsv FILE RETURNED\n\n")
                    }}
                    germ_tree_ch2 = get_germ_tree.out.germ_tree_ch.collectFile(name: "seq_for_germ_tree.tsv", skip: 1, keepHeader: true)
                    // germ_tree_ch2.subscribe{it -> it.copyTo("${out_path}/tsv")}
                    germ_tree_ch3 = get_germ_tree.out.germ_tree_ch.flatten().filter{file -> file.countLines() > align_clone_nb.toInteger() }

                    get_germ_tree.out.no_cloneID_ch.count().subscribe {n -> if ( n == 0 ){
                        print("\n\nWARNING:\nALL SEQUENCES IN CLONAL GROUP FOLLOWING THE get_germ_tree PROCESS -> EMPTY germ_tree_dismissed_clone_id.tsv FILE RETURNED\n\n")
                    }}
                    no_cloneID_ch2 = get_germ_tree.out.no_cloneID_ch.collectFile(name: "germ_tree_dismissed_clone_id.tsv")
                    no_cloneID_ch2.subscribe{it -> it.copyTo("${out_path}/tsv")}

                    get_germ_tree.out.cloneID_ch.count().subscribe {n -> if ( n == 0 ){
                        print("\n\nWARNING:\nNO CLONAL GROUP FOLLOWING THE get_germ_tree PROCESS -> EMPTY germ_tree_clone_id.tsv and germ_tree.pdf FILES RETURNED\n\n")
                    }}
                    cloneID_ch2 = get_germ_tree.out.cloneID_ch.collectFile(name: "germ_tree_clone_id.tsv")
                    cloneID_ch2.subscribe{it -> it.copyTo("${out_path}/tsv")}

                    get_germ_tree.out.get_germ_tree_log_ch.collectFile(name: "get_germ_tree.log").subscribe{it -> it.copyTo("${out_path}/reports")} // 



                    germ_tree_vizu(
                        rdata_germ_tree_ch2,
                        germ_tree_kind,
                        align_clone_nb,
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
                        clone_assigned_seq_ch, // may be add .first()
                        meta_file, 
                        meta_legend,
                        cute_file
                    )
                    germ_tree_vizu.out.germ_tree_vizu_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE germ_tree_vizu PROCESS\n\n========\n\n"}
                    //germ_tree_vizu.out.germ_tree_vizu_ch.count().subscribe {n -> if ( n == 0 ){error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE germ_tree_vizu PROCESS\n\n========\n\n"}}
                    germ_tree_vizu.out.germ_tree_dup_seq_not_displayed_ch.ifEmpty{
                        print("\n\nWARNING:\n-> NO germ_tree_dup_seq_not_displayed.tsv FILE RETURNED\n\n")
                    }
                    germ_tree_dup_seq_not_displayed_ch2 = germ_tree_vizu.out.germ_tree_dup_seq_not_displayed_ch.flatten().collectFile(name: "germ_tree_dup_seq_not_displayed.tsv", skip: 1, keepHeader: true) // flatten split the list into several objects which is required by collectFile()
                    germ_tree_dup_seq_not_displayed_ch2.subscribe{it -> it.copyTo("${out_path}/tsv")}

                    */

                // end when: nb_clone_assignment > 0

                tempo1_ch = Channel.of("all", "annotated", "clone") // 1 channel with 3 values (not list)
                tempo2_ch = Data_assembly.out.wanted_ch.mix(Data_assembly.out.wanted_ch.mix(clone_assigned_seq_ch)) // 1 channel with 3 paths (do not use flatten() -> not list)
                tempo3_ch = tempo1_ch.merge(tempo2_ch) // 3 lists
                tempo4_ch = Channel.of("vj_allele", "c_allele", "vj_gene", "c_gene")
                tempo5_ch = tempo3_ch.combine(tempo4_ch) // 12 tuples
                Donut(
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
                    cute_file,
                    Igblast_data_check.out.allele_names_tsv_all_ch
                )
                donut_tsv_ch2 = Donut.out.donut_tsv_ch.collectFile(name: "donut_stats.tsv", skip: 1, keepHeader: true).subscribe{it -> it.copyTo("${out_path}/tsv")} // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files
                donuts_png_ch = donuts_png_ch.mix(Donut.out.donuts_png_ch.collect())

                Donut_assembly(
                    Donut.out.donut_pdf_ch.collect()
                )



                clone_labeled_ch = clone_assigned_seq_filtered_ch.map {file -> tuple(file, 'CLONE') }
                wanted_labeled_ch = Data_assembly.out.wanted_ch.map {file -> tuple(file, 'ALL') }
                imgt_labeled_ch = Data_assembly.out.wanted_ch.map {file -> tuple(file, 'IMGT') }
                all_files_ch = clone_labeled_ch.mix(wanted_labeled_ch, imgt_labeled_ch)
                Tsv2fasta(
                    all_files_ch,
                    align_seq, 
                    clone_germline_kind, 
                    align_clone_nb, 
                    cute_path
                )
                // fasta_align_imgt_ch = Tsv2fasta.out.fasta_align_ch.ifEmpty{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\nEMPTY OUTPUT FOLLOWING THE Tsv2fasta PROCESS\n\n========\n\n"}.filter {nuc, aa, kind -> nuc.name.endsWith("_imgt_nuc.fasta")}
                Tsv2fasta.out.fasta_align_ch.filter{nuc, aa, kind -> nuc.exists() && aa.exists() }.subscribe{nuc, aa, kind -> nuc.copyTo("${out_path}/fasta/for_alignment_nuc/${nuc.getName()}") ; aa.copyTo("${out_path}/fasta/for_alignment_aa/${aa.getName()}")} // copy the folder and content for_alignment_nuc/* into {out_path}/fasta
                Tsv2fasta.out.fasta_align_imgt_ch.filter{nuc, aa, kind -> nuc.exists() && aa.exists() }.subscribe{nuc, aa, kind -> nuc.copyTo("${out_path}/alignments/nuc/imgt/${nuc.getName()}") ; aa.copyTo("${out_path}/alignments/aa/imgt/${aa.getName()}")}
                fasta_align_imgt_aa_ch = Tsv2fasta.out.fasta_align_imgt_ch.map{nuc, aa, kind  -> [aa, nuc, kind] } // for aa printing into html
                // Tsv2fasta.out.fasta_align_imgt_ch.filter{nuc, aa, kind -> nuc.exists() }.subscribe{nuc, aa, kind -> nuc.copyTo("${out_path}/alignments/nuc/imgt/${nuc.getName()}")} // aa.copyTo("${out_path}/alignments/aa/${aa.getName()}") not used because no gaps in this aa sequences
                // fasta_align_ch2 = Tsv2fasta.out.fasta_align_ch.filter {nuc, aa, kind -> !nuc.name.endsWith("_imgt_nuc.fasta")}
                copyLogFile('Tsv2fasta.log', Tsv2fasta.out.tsv2fasta_log_ch, out_path)
                // Print warnings on the terminal:
                warning_ch = warning_ch.mix(Tsv2fasta.out.warning_ch.filter{file(it).exists()}.map{file -> file.text}) //  file.text = contenu du fichier
                has_fasta_align_imgt = Tsv2fasta.out.fasta_align_imgt_ch.map{true }.ifEmpty {[false] }.first() // has_fasta_align_imgt is false if Tsv2fasta.out.fasta_align_imgt_ch is empty



                if(align_soft == "abalign" && (align_seq == "query" || align_seq == "igblast_full" || align_seq == "trimmed" || align_seq == "fwr1" || align_seq == "fwr2" || align_seq == "fwr3" || align_seq == "fwr4" || align_seq == "cdr1" || align_seq == "cdr2" || align_seq == "cdr3" || align_seq == "junction" || align_seq == "d_sequence_alignment" || align_seq == "j_sequence_alignment" || align_seq == "c_sequence_alignment" || align_seq == "d_germline_alignment" || align_seq == "j_germline_alignment" || align_seq == "c_germline_alignment")){
                    align_soft = "mafft"
                    warn = "\n\nWARNING:\nalign_soft PARAMETER OF THE nextflow.config FILE RESET TO \"mafft\" SINCE align_soft PARAMETER WAS SET TO \"abalign\" BUT Abalign EITHER 1) REQUIRES AT LEAST A V DOMAIN IN THE SEQUENCES OR 2) TRUNKS THE C CONSTANT REGION,\nWHILE align_seq PARAMETER IS SET TO \"${align_seq}\"\n\n"
                    print(warn)
                    warning_ch = warning_ch.mix(Channel.value(warn)) // accumulate
                }
                if(align_soft == "mafft"){

                    // fasta_align_ch2 = Tsv2fasta.out.fasta_align_ch.map{nuc, aa, tag -> 
                    //     n = nuc.text.readLines().size()
                    //    return tuple(nuc, aa, tag, n) // accumulate
                    // }

                    Mafft_align(
                        Tsv2fasta.out.fasta_align_ch,
                        align_mafft_all_options,
                        align_mafft_clonal_options
                    )
                    align_aa_ch = Mafft_align.out.aligned_all_ch.map{nuc, aa, tag -> [aa, nuc, tag] }
                    align_nuc_ch = Mafft_align.out.aligned_all_ch.map{nuc, aa, tag -> [nuc, aa, tag] }
                    aligned_all_ch2 = Mafft_align.out.aligned_all_ch.map{nuc, aa, tag -> [nuc, aa] }
                    copyLogFile('mafft_align.log', Mafft_align.out.mafft_align_log_ch, out_path)


                }else if(align_soft == "abalign"){


                    Abalign_align_aa(
                        Tsv2fasta.out.fasta_align_ch,
                        igblast_organism,
                        igblast_B_heavy_chain,
                        align_abalign_options
                    )
                    copyLogFile('abalign_align_aa.log', Abalign_align_aa.out.abalign_align_aa_log_ch, out_path)



                    Abalign_rename(
                        Abalign_align_aa.out.aligned_aa_c
                    )
                    Abalign_rename.out.failed_abalign_align_ch.collectFile(name: "failed_abalign_align.tsv", skip: 1, keepHeader: true).subscribe{it -> it.copyTo("${out_path}/tsv")}
                    warning_ch = warning_ch.mix(Abalign_rename.out.abalign_rename_warn_ch.filter{file(it).exists()}.map{file -> file.text}) //  file.text = contenu du fichier) // accumulate
                    copyLogFile('abalign_rename.log', Abalign_rename.out.abalign_rename_log_ch, out_path)



                    Abalign_align_nuc(
                        Abalign_rename.out.renamed_aligned_aa_ch
                    )
                    align_nuc_ch = Abalign_align_nuc.out.aligned_all_ch.map{nuc, aa, tag -> [nuc, aa, tag] }
                    align_aa_ch = Abalign_align_nuc.out.aligned_all_ch.map{nuc, aa, tag -> [aa, nuc, tag] }
                    aligned_all_ch2 = Abalign_align_nuc.out.aligned_all_ch.map{nuc, aa, tag -> [nuc, aa] }
                    copyLogFile('abalign_align_nuc.log', Abalign_align_nuc.out.abalign_align_nuc_log_ch, out_path)



                }else{
                    error "\n\n========\n\nINTERNAL ERROR IN NEXTFLOW EXECUTION\n\nINVALID align_soft PARAMETER IN nextflow.config FILE:\n${align_soft}\n\n========\n\n"
                }
                align_nuc_ch.subscribe{nuc, aa, tag -> if(tag == "CLONE"){nuc.copyTo("${out_path}/alignments/nuc/clonal/${nuc.getName()}")}else if(tag == "ALL"){nuc.copyTo("${out_path}/alignments/nuc/all/${nuc.getName()}")}else{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\ntag CANNOT BE OTHER THAN ALL OR CLONE HERE.\n\n========\n\n"}}
                align_aa_ch.subscribe{aa, nuc, tag -> if(tag == "CLONE"){aa.copyTo("${out_path}/alignments/aa/clonal/${aa.getName()}")}else if(tag == "ALL"){aa.copyTo("${out_path}/alignments/aa/all/${aa.getName()}")}else{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\ntag CANNOT BE OTHER THAN ALL OR CLONE HERE.\n\n========\n\n"}}


                //Abalign_align_nuc.out.aligned_all_ch.map{nuc, aa, tag -> [y, z] }.view()
                //Mutation_load_germ_genes.out.mutation_load_ch.view()
                branches_nuc = align_nuc_ch.branch {
                    CLONE: it[2] == 'CLONE'
                    ALL  : it[2] == 'ALL'
                }
                clone_for_gff_nuc_ch = branches_nuc.CLONE.combine(clone_assigned_seq_filtered_ch.first())
                all_for_gff_nuc_ch  = branches_nuc.ALL.combine(Data_assembly.out.wanted_ch.first())
                for_gff_nuc_ch = clone_for_gff_nuc_ch.mix(all_for_gff_nuc_ch)
                branches_aa = align_aa_ch.branch {
                    CLONE: it[2] == 'CLONE'
                    ALL  : it[2] == 'ALL'
                }
                clone_for_gff_aa_ch = branches_aa.CLONE.combine(clone_assigned_seq_filtered_ch.first())
                all_for_gff_aa_ch  = branches_aa.ALL.combine(Data_assembly.out.wanted_ch.first())
                for_gff_aa_ch = clone_for_gff_aa_ch.mix(all_for_gff_aa_ch)



                GffNuc( // module gff.nf
                    for_gff_nuc_ch,
                    align_seq, 
                    align_clone_nb, 
                    cute_path
                )
                GffNuc.out.gff_ch.subscribe{gff_list, tag -> if(tag == "CLONE"){gff_list.each{gff -> gff.copyTo("${out_path}/alignments/nuc/clonal/${gff.getName()}")}}else if(tag == "ALL"){gff_list.each{gff -> gff.copyTo("${out_path}/alignments/nuc/all/${gff.getName()}")}}else{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\ntag CANNOT BE OTHER THAN ALL OR CLONE HERE.\n\n========\n\n"}} // .flatten() because gff several files. Otherwise, copyTo does not like it. Finally, .each() also solve the problem
                copyLogFile('gffNuc.log', GffNuc.out.gff_log_ch, out_path)
                warning_ch = warning_ch.mix(GffNuc.out.gff_warn_ch.filter{file(it).exists()}.map{file -> file.text}) //  file.text = contenu du fichier


                //for_gff_aa_ch.view()
                GffAa( // module gff.nf
                    for_gff_aa_ch,
                    align_seq, 
                    align_clone_nb, 
                    cute_path
                )
                GffAa.out.gff_ch.subscribe{gff_list, tag -> if(tag == "CLONE"){gff_list.each{gff -> gff.copyTo("${out_path}/alignments/aa/clonal/${gff.getName()}")}}else if(tag == "ALL"){gff_list.each{gff -> gff.copyTo("${out_path}/alignments/aa/all/${gff.getName()}")}}else{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\ntag CANNOT BE OTHER THAN ALL OR CLONE HERE.\n\n========\n\n"}}
                copyLogFile('gffAa.log', GffAa.out.gff_log_ch, out_path)
                warning_ch = warning_ch.mix(GffAa.out.gff_warn_ch.filter{file(it).exists()}.map{file -> file.text}) //  file.text = contenu du fichier
                GffAa.out.gff_approx_coord_ch.subscribe{tsv, tag -> if(tag == "CLONE"){tsv.copyTo("${out_path}/alignments/aa/clonal/${tsv.getName()}")}else if(tag == "ALL"){tsv.copyTo("${out_path}/alignments/aa/all/${tsv.getName()}")}else{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\ntag CANNOT BE OTHER THAN ALL OR CLONE HERE.\n\n========\n\n"}} // .flatten() because gff several files. Otherwise, copyTo does not like it. Finally, .each() also solve the problem
                //GffAa.out.gff_approx_coord_ch.collectFile(name: "gff_aa_approx_coord.tsv", skip: 1, keepHeader: true).subscribe{it -> it.copyTo("${out_path}/tsv")} // warning: skip: 1, keepHeader: true means that if the first file of the list is empty, then it is taken as reference to do not remove the header -> finally no header in the returned fusioned files


                PrintAlignmentNuc( // module print_alignment.nf
                    align_nuc_ch
                )
                PrintAlignmentNuc.out.alignment_html.subscribe{html, tag -> if(tag == "CLONE"){html.copyTo("${out_path}/alignments/nuc/clonal/${html.getName()}")}else if(tag == "ALL"){html.copyTo("${out_path}/alignments/nuc/all/${html.getName()}")}else{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\ntag CANNOT BE OTHER THAN ALL OR CLONE HERE.\n\n========\n\n"}}
                
                //.subscribe{it -> it.copyTo("${out_path}/alignments/nuc")}



                PrintAlignmentAa( // module print_alignment.nf
                    align_aa_ch
                )
                PrintAlignmentAa.out.alignment_html.subscribe{html, tag -> if(tag == "CLONE"){html.copyTo("${out_path}/alignments/aa/clonal/${html.getName()}")}else if(tag == "ALL"){html.copyTo("${out_path}/alignments/aa/all/${html.getName()}")}else{error "\n\n========\n\nERROR IN NEXTFLOW EXECUTION\n\ntag CANNOT BE OTHER THAN ALL OR CLONE HERE.\n\n========\n\n"}}



                if(has_fasta_align_imgt && (align_seq == "query" || align_seq == "igblast_full" || align_seq == "trimmed" || align_seq == "sequence_alignment")){


                    Gff_imgt( // module gff_imgt.nf
                        Tsv2fasta.out.fasta_align_imgt_ch, 
                        align_seq, 
                        Data_assembly.out.wanted_ch, 
                        cute_path
                    )
                    copyLogFile('gff_imgt.log', Gff_imgt.out.gff_log_ch, out_path)



                    PrintAlignmentIMGTnuc( // module print_alignment.nf
                        Tsv2fasta.out.fasta_align_imgt_ch
                    )
                    PrintAlignmentIMGTnuc.out.alignment_html.subscribe{html, tag -> html.copyTo("${out_path}/alignments/nuc/imgt")}



                    PrintAlignmentIMGTaa( // module print_alignment.nf
                        fasta_align_imgt_aa_ch
                    )
                    PrintAlignmentIMGTaa.out.alignment_html.subscribe{html, tag -> html.copyTo("${out_path}/alignments/aa/imgt")}


                }

                (aligned_all_nuc_ch, aligned_all_aa_ch) = aligned_all_ch2.multiMap { nuc, aa ->
                    nuc: nuc
                    aa: aa
                }

                Tree_nuc(
                    aligned_all_nuc_ch,
                    phylo_tree_model_file
                )

                Tree_aa(
                    aligned_all_aa_ch,
                    phylo_tree_model_file
                )


                if(phylo_tree_itol_subscription == "TRUE"){ // Warning: no run if tree does not exist (because tree_nuc_ch and tree_aa_ch are optional)

                    Meta2Itol (
                        meta_file,
                        meta_seq_names
                    )

                    tree_aa = Tree.out.tree_aa_ch
                    tree_nuc = Tree.out.tree_nuc_ch
                    //copyLogFile('tree.log', Tree.out.tree_log_ch, out_path)
                    tree = tree_aa.concat(tree_nuc)



                    // The ITOL process can only be executed if user has paid the subsription for automated visualization
                    Itol(
                        tree, // Warning: no run if tree does not exist (because tree_nuc_ch and tree_aa_ch are optional)
                        Meta2Itol.out.itol_out_ch,
                        phylo_tree_itolkey
                    )
                    itol_aa_ch = Itol.out.itol_aa_ch
                    itol_nuc_ch = Itol.out.itol_nuc_ch

                    itol_aa_ch.subscribe{it -> it.copyTo("${out_path}/phylo/aa")}
                    itol_nuc_ch.subscribe{it -> it.copyTo("${out_path}/phylo/nuc")}

                }


            // end when: nb_wanted > 0

        // end when: nb_productive > 0

    // end when: nb_igblast > 0
    Print_warnings(
        warning_ch.ifEmpty{''}.collectFile(name: "warnings_collect.txt")
    )


    Print_report(
        config_file, // from parameter
        template_rmd, // from parameter
        alignments_viz_rmd, // from parameter
        alignments_viz_html, // from parameter
        igblast_organism, // parameter
        igblast_loci, // parameter
        igblast_B_heavy_chain, // parameter
        igblast_B_lambda_chain, // parameter
        igblast_B_kappa_chain, // parameter
        igblast_T_alpha_chain, // parameter
        igblast_T_beta_chain, // parameter
        igblast_T_gamma_chain, // parameter
        igblast_T_delta_chain, // parameter
        clone_strategy, // parameter
        clone_model, // parameter
        clone_normalize, // parameter
        clone_distance, // parameter
        clone_germline_kind, // parameter
        clone_mut_obs_seq, // parameter
        clone_mut_germ_seq, // parameter
        clone_mut_regionDefinition, // parameter
        align_clone_nb, // parameter
        align_soft, // parameter
        align_seq, // parameter
        align_abalign_options, // parameter
        align_mafft_all_options, // parameter
        align_mafft_clonal_options, // parameter
        nb_input, // mandatory
        nb_igblast, // mandatory
        nb_unigblast, // mandatory
        nb_productive.ifEmpty{'EMPTY'}, 
        nb_unproductive.ifEmpty{'EMPTY'}, 
        nb_wanted.ifEmpty{'EMPTY'}, 
        nb_unwanted.ifEmpty{'EMPTY'}, 
        nb_dist_ignored.ifEmpty{'EMPTY'}, 
        nb_clone_tot.ifEmpty{'EMPTY'}, 
        nb_unclone_tot.ifEmpty{'EMPTY'}, 
        nb_clone_unassignment.ifEmpty{'EMPTY'}, 
        nb_clone_ungermline.ifEmpty{'EMPTY'}, 
        distance_hist_ch.ifEmpty{[empty_distance_hist]}, 
        Donut.out.donuts_png_ch.collect().ifEmpty{ [empty_donuts_png] }, 
        Repertoire.out.repertoire_png_ch.collect().ifEmpty{ [empty_repertoire_png] }, 
        repertoire_constant_ch.ifEmpty{'EMPTY'}, 
        repertoire_vj_ch.ifEmpty{'EMPTY'}, 
        phylo_tree_itol_subscription, // parameter
        Print_warnings.out.final_warning_ch // just so that print_report wait for all warnings // warning_ch.collect().map{it.join('\n\n')}.ifEmpty{'EMPTY'} // concatenate all warnings into a single string // finally, the gathered string is very loong. I prefer to use a file added in /reports/
    )


    Backup(
        config_file, 
        log_file
    )

}

    //////// end Main


//////// end Processes