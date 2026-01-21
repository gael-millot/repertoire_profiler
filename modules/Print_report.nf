// The render function only creates the html file with the rmardown
// Inputs:
//      template_rmd: rmardown template file used to create the html (file path is defined in config.nextflow)
//      nb_input: number of fasta sequences in initial input
//      nb_igblast: number of sequences that igblast could analyse
//      nb_unigblast: number of sequences that igblast could not alayse
//      donuts_png: provide links to the donut images needed in the rmd inside the work folder
//      repertoire_png: provide links to the repertoire images needed in the rmd inside the work folder
//      repertoire_constant_ch: names of the constant gene repertoire files to be displayed
//      repertoire_vj_ch: names of the variable gene repertoire files to be displayed
//      itol_subscription: nextflow.config parameter to know if user has paid the subscription to itol automated visualization of trees, process ITOL is only executed if TRUE
//      heavy_chain: to know if the analyzed data is VL or VH, because "Amino acid sequences phylogeny" section in html report is only displayed for VH
// Outputs:
//      "report.html": finalized html report for a specific run
//      "print_report.log": will contain any error or warning messages produced by rmardown::render
process Print_report{
    label 'r_ig_clustering'

    publishDir path: "${out_path}", mode: 'copy', pattern: "{report.html}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{alignments_viz.html}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{print_report.log}", overwrite: false
    cache 'false'

    input:
    file template_rmd
    file alignments_viz_rmd
    file alignments_viz_html
    val nb_input
    val nb_igblast
    val nb_unigblast
    val nb_productive
    val nb_unproductive
    val nb_wanted
    val nb_unwanted
    val nb_dist_ignored
    val nb_clone_assigned
    val nb_failed_clone
    val nb_failed_clone_assignment
    val nb_failed_clone_germline
    path distance_hist_ch
    path donuts_png
    path repertoire_png
    val repertoire_constant_ch
    val repertoire_vj_ch
    val clone_distance
    val align_soft
    val itol_subscription
    val warning

    output:
    file "report.html"
    file "alignments_viz.html"
    file "print_report.log"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    cp ${template_rmd} report_file.rmd
    cp ${alignments_viz_rmd} alignments_vizu.rmd
    cp ${alignments_viz_html} alignments_viz.html
    cp -r "${out_path}/phylo" .
    cp -r "${out_path}/tsv" .
    cp -r "${out_path}/pdf" .
    cp -r "${projectDir}/bin/doc_images" .
    Rscript -e '

        # Find the constant and vj repertoires to be displayed in the html file (file names differ depending on light/heavy chain)
        constant_files <- as.character("${repertoire_constant_ch}")
        vj_files <- as.character("${repertoire_vj_ch}")
        cleaned_repertoire_constant <- gsub("^\\\\[\\\\[|\\\\]\\\\]\$", "", constant_files)
        cleaned_repertoire_vj <- gsub("^\\\\[\\\\[|\\\\]\\\\]\$", "", vj_files)
        constant_paths <- strsplit(cleaned_repertoire_constant, ",\\\\s*")[[1]]
        vj_paths <- strsplit(cleaned_repertoire_vj, ",\\\\s*")[[1]]
        constant_names <- basename(constant_paths)
        vj_names <- basename(vj_paths)
        # Verification that the resulting file names are as expected
        constant_rep <- constant_names[grepl("^IG.C_.*gene_non-zero\\\\.png\$", constant_names)]
        vj_rep <- vj_names[grepl("^rep_gene_IG.V_.*non-zero\\\\.png\$", vj_names)]
        if(length(constant_rep) == 0 || length(vj_rep) == 0){
            stop(paste0("\\n\\n========\\n\\nERROR IN print_report PROCESS\\n\\nTHE REPERTOIRE PNG FILES TO BE DISPLAYED WERE NOT FOUND\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n"), call. = FALSE)
        }

        rmarkdown::render(
        input = "report_file.rmd",
        output_file = "report.html",
        # list of the variables waiting to be replaced in the rmd file:
        params = list(
            nb_input = ${nb_input},
            nb_igblast = ${nb_igblast}, 
            nb_unigblast = ${nb_unigblast},
            nb_productive = ${nb_productive},
            nb_unproductive = ${nb_unproductive},
            nb_wanted =  ${nb_wanted},
            nb_unwanted =  ${nb_unwanted},
            nb_dist_ignored = ${nb_dist_ignored},
            nb_clone_assigned = ${nb_clone_assigned},
            nb_failed_clone = ${nb_failed_clone},
            nb_failed_clone_assignment = ${nb_failed_clone_assignment}, 
            nb_failed_clone_germline = ${nb_failed_clone_germline}, 
            clone_distance = ${clone_distance},
            constant_rep = constant_rep,
            vj_rep = vj_rep,
            align_soft = "${align_soft}",
            itol_subscription = ${itol_subscription},
            warning_collect = ${warning}
        ),
        # output_dir = ".",
        # intermediates_dir = "./",
        # knit_root_dir = "./",
        run_pandoc = TRUE,
        quiet = TRUE,
        clean = TRUE
        )

        html_here_ok <- TRUE # set to FALSE to rerun the creation of the alignments_viz_html file
        if(html_here_ok == FALSE){
            # dir.create("reports", showWarnings = FALSE, recursive = TRUE)
            rmarkdown::render(input = "alignments_vizu.rmd", output_file = "alignments_viz.html", run_pandoc = TRUE, quiet = TRUE, clean = TRUE)
        }
    ' |& tee -a print_report.log
    """
}
