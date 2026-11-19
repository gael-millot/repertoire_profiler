#!/usr/bin/Rscript
#########################################################################
##                                                                     ##
##     donut.R                                                         ##
##                                                                     ##
##     Gael A. Millot                                                  ##
##     Chloe Taurel                                                    ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################




################################ Aim


################################ End Aim


################################ Introduction


################################ End Introduction


################################ Acknowlegments


################################ End Acknowlegments


################################ Initialization


# R version checking
if(version$version.string != "R version 4.4.2 (2024-10-31)"){
    stop(paste0("\n\n================\n\nERROR IN donut.R\n\n\n", version$version.string, " IS NOT THE 4.4.2 RECOMMANDED\n\n================\n\n"))
}
# other initializations
erase.objects = TRUE # write TRUE to erase all the existing objects in R before starting the algorithm and FALSE otherwise. Beginners should use TRUE
if(erase.objects == TRUE){
    rm(list = ls(all.names = TRUE))
    erase.objects = TRUE
}
erase.graphs = TRUE # write TRUE to erase all the graphic windows in R before starting the algorithm and FALSE otherwise
script <- "donut"


################################ End Initialization


################################ Parameters that need to be set by the user


################################ End Parameters that need to be set by the user


################################ Config import


tempo.cat <- "KIND OF RUN (SCRIPT, COPY-PASTE OR SOURCE): "
if(interactive() == FALSE){ # if(grepl(x = commandArgs(trailingOnly = FALSE), pattern = "R\.exe$|\/R$|Rcmd\.exe$|Rcmd$|Rgui\.exe$|Rgui$|Rscript\.exe$|Rscript$|Rterm\.exe$|Rterm$")){ # detection of script usage
    run.way <- "SCRIPT"
    cat(paste0("\n\n", tempo.cat, run.way, "\n"))
    command <- paste0(commandArgs(trailingOnly = FALSE), collapse = ",") # recover the full command
    args <- commandArgs(trailingOnly = TRUE) # recover arguments written after the call of the R script
    if(any(is.na(args))){
        stop(paste0("\n\n================\n\nERROR IN donut.R\n\n\nTHE args OBJECT HAS NA\n\n================\n\n"), call. = FALSE)
    }
    tempo.arg.names <- c(
        "file_name", 
        "kind", 
        "col",
        "donut_palette",
        "donut_hole_size",
        "donut_hole_text",
        "donut_hole_text_size",
        "donut_border_color",
        "donut_border_size",
        "donut_annotation_distance",
        "donut_annotation_size",
        "donut_annotation_force",
        "donut_annotation_force_pull",
        "donut_legend_width",
        "donut_legend_text_size", 
        "donut_legend_box_size", 
        "donut_legend_box_space", 
        "donut_legend_limit",
        "cute", 
        "repertoire_names_ch",
        "color_ref_files",
        "log"
    ) # objects names exactly in the same order as in the bash code and recovered in args. Here only one, because only the path of the config file to indicate after the donut.R script execution
    if(length(args) != length(tempo.arg.names)){
        stop(paste0("\n\n================\n\nERROR IN donut.R\n\n\nTHE NUMBER OF ELEMENTS IN args (", length(args),") IS DIFFERENT FROM THE NUMBER OF ELEMENTS IN tempo.arg.names (", length(tempo.arg.names),")\nargs:", paste0(args, collapse = ","), "\ntempo.arg.names:", paste0(tempo.arg.names, collapse = ","), "\n\n================\n\n"), call. = FALSE)
    }
    for(i1 in 1:length(tempo.arg.names)){
        assign(tempo.arg.names[i1], args[i1])
    }
    rm(args, i1)
}else if(sys.nframe() == 0L){ # detection of copy-paste/direct execution (for debugging). With script it is also 0, with source, it is 4
    run.way <- "COPY-PASTE"
    cat(paste0("\n\n", tempo.cat, run.way, "\n"))
}else{
    run.way <- "SOURCE" # using source(), sys.nframe() is 4
    cat(paste0("\n\n", tempo.cat, run.way, "\n"))
}
rm(tempo.cat)


################################ End Config import

################################ Test

# setwd("C:/Users/gmillot/Documents/Git_projects/repertoire_profiler/work/b7/b38b134c2cb42833aa1998e96051e8")
# file_name = "./wanted_seq.tsv"
# kind = "annotated"
# col = "c_allele"
# donut_palette = "NULL" 
# donut_hole_size = "0.5" 
# donut_hole_text = "TRUE" 
# donut_hole_text_size = "14" 
# donut_border_color = "gray50" 
# donut_border_size = "0.1" 
# donut_annotation_distance = "0" 
# donut_annotation_size = "3" 
# donut_annotation_force = "1" 
# donut_annotation_force_pull = "100" 
# donut_legend_width = "0.25"
# donut_legend_text_size = "10" 
# donut_legend_box_size = "5" 
# donut_legend_box_space = "2" 
# donut_legend_limit = "0.1"
# cute = "C:/Users/gmillot/Documents/Git_projects/repertoire_profiler/bin/cute_little_R_functions_v12.8.R"
# repertoire_names_ch = "imgt_human_IGHC2.tsv imgt_human_IGHD2.tsv imgt_human_IGHJ2.tsv imgt_human_IGHV2.tsv"
# color_ref_files = "color_ref/wanted_seq.tsv color_ref/clone_assigned_seq.tsv"
# log = "all_donut.log"



################################ end Test

################################ Recording of the initial parameters


param.list <- c(
    "erase.objects", 
    "erase.graphs", 
    "script", 
    "run.way",
    "tempo.arg.names", 
    if(run.way == "SCRIPT"){"command"}, 
    "file_name", 
    "kind", 
    "col",
    "donut_palette",
    "donut_hole_size",
    "donut_hole_text",
    "donut_hole_text_size",
    "donut_border_color",
    "donut_border_size",
    "donut_annotation_distance",
    "donut_annotation_size",
    "donut_annotation_force",
    "donut_annotation_force_pull",
    "donut_legend_width",
    "donut_legend_text_size", 
    "donut_legend_box_size", 
    "donut_legend_box_space", 
    "donut_legend_limit",
    "cute", 
    "repertoire_names_ch",
    "color_ref_files",
    "log"
)
if(any(duplicated(param.list))){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 1 IN donut.R\n\nTHE param.list OBJECT CONTAINS DUPLICATED ELEMENTS:\n", paste(param.list[duplicated(param.list)], collapse = " "), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE) # message for developers
}
if(erase.objects == TRUE){
    created.object.control <- ls()[ ! ls() %in% "param.list"]
    if( ! (all(created.object.control %in% param.list) & all(param.list %in% created.object.control))){
        stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 2 IN donut.R\n\nINCONSISTENCIES BETWEEN THE ARGUMENTS USED AND THE PARAMETERS REQUIRED IN THE EXECUTABLE CODE FILE\nTHE ARGUMENTS NOT PRESENT IN THE EXECUTABLE FILE (donut.R) ARE:\n", paste(created.object.control[ ! created.object.control %in% param.list], collapse = " "), "\nTHE PARAMETERS OF THE EXECUTABLE FILE (donut.R) NOT PRESENT IN THE ARGUMENTS ARE:\n", paste(param.list[ ! param.list %in% created.object.control], collapse = " "), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE) # message for developers
    }
}
char.length <- nchar(param.list)
space.add <- max(char.length) - char.length + 5
param.ini.settings <- character(length = length(param.list))
for(i in 1:length(param.list)){
    param.ini.settings[i] <- paste0("\n", param.list[i], paste0(rep(" ", space.add[i]), collapse = ""), paste0(get(param.list[i]), collapse = ",")) # no env = sys.nframe(), inherit = FALSE in get() because look for function in the classical scope
}


################################ End Recording of the initial parameters


################################ Functions


# Functions are built such that they should have no direct use of Global objects (going through the R scope), and only use function arguments
# 1) Cute little function is sourced for the moment into the .GlobalEnv environment, but may be interesting to put it into a new environement just above .GlobalEnv environment. See https://stackoverflow.com/questions/9002544/how-to-add-functions-in-an-existing-environment
# 2) Argument names of each function must not be a name of Global objects (error message otherwise)
# 3) Argument name of each function ends with "_fun" in the first function, "_2fun" in the second, etc. This prevent conflicts with the argument partial names when using these functions, notably when they are imbricated

################ import functions from cute little functions toolbox


if(length(cute) != 1){
    stop(paste0("\n\n============\n\nERROR IN donut.R\n\ncute PARAMETER MUST BE LENGTH 1: ", paste(cute, collapse = " "), "\n\n============\n\n"), call. = FALSE)
}else if(grepl(x = cute, pattern = "^http")){
    tempo.try <- try(suppressWarnings(suppressMessages(source(cute, local = .GlobalEnv))), silent = TRUE)
    if(any(grepl(x = tempo.try, pattern = "^[Ee]rror"))){
        stop(paste0("\n\n============\n\nERROR IN donut.R\n\nHTTP INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else if( ! grepl(x = cute, pattern = "^http")){
    if( ! file.exists(cute)){
        stop(paste0("\n\n============\n\nERROR IN donut.R\n\nFILE INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else{
    tempo.cat <- paste0("\n\n================\n\nINTERNAL CODE ERROR 3 IN donut.R: CODE HAS TO BE MODIFIED\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n============\n\n")
    stop(tempo.cat, call. = FALSE)
}


# required cute function checking
req.function <- c(
    "fun_check",
    "fun_pack", 
    "fun_report"
)
tempo <- NULL
for(i1 in req.function){
    if(length(find(i1, mode = "function")) == 0L){
        tempo <- c(tempo, i1)
    }
}
if( ! is.null(tempo)){
    tempo.cat <- paste0("ERROR IN donut.R\n\nREQUIRED cute FUNCTION", ifelse(length(tempo) > 1, "S ARE", " IS"), " MISSING IN THE R ENVIRONMENT:\n", paste0(tempo, collapse = "()\n"))
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop(), not in tempo.cat, to be able to add several messages between ==
}
# end required function checking


################ local function: package import


# R Packages required
req.package.list <- c(
    "lubridate", 
    "ggplot2",
    "ggrepel"
)
for(i in 1:length(req.package.list)){suppressMessages(library(req.package.list[i], character.only = TRUE))}
# fun_pack(req.package = req.package.list, load = TRUE, lib_path = NULL) # packages are imported even if inside functions are written as package.name::function() in the present code

################ local function


# These functions extract the loci from a igblast_ref_file
# Example of input : "imgt_human_IGLV.fasta imgt_human_IGKV.fasta"
# The output would be ["IGL", "IGK"]
extract_locus <- function(x) {
    parts <- strsplit(x, "_")[[1]]
    ig_part <- parts[grepl("^IG[A-Z]{2}", parts)]
    return(substr(ig_part, 1, 3))
}
extract_chain <- function(x) {
    files_list <- unlist(strsplit(x, " "))
    chain <- sapply(files_list, extract_locus)
    return(chain)
}


# Class table of a donut, i.e., number of sequences per allele or gene (V and J pair for vj_allele and vj_gene)
# Used for the plotted donut and for the donuts of the other kinds, so that a class has the same color in all the donuts
# obs_fun: data frame of the wanted_seq.tsv file (kinds "all" and "annotated") or of the clone_assigned_seq.tsv file (kind "clone"), with at least one row
# kind_fun: kind of donut, "all", "annotated" or "clone"
# col_fun: classes of the donut, "vj_allele", "c_allele", "vj_gene" or "c_gene"
# Returns a list:
# $obs2: data frame of the classes of all the sequences (first column named col_fun, then Freq, Prop, kind, and labels if sequence names have been replaced), in decreasing order of Prop
# $obs3: obs2 restricted to the classes with at least one replaced sequence name if kind_fun is "annotated", obs2 otherwise
# $all.annotation.log: TRUE if no sequence name has been replaced (no annotation)
# $nb_na: number of NA among the selected sequences
# $nb_total: number of selected sequences
donut_classes <- function(obs_fun, kind_fun, col_fun){
    if(kind_fun == "clone"){
        if (col_fun == "vj_allele") {
            tempo.primary <- obs_fun$germline_v_call
            tempo.secondary <- obs_fun$germline_j_call
        } else if (col_fun == "vj_gene") {
            tempo.primary <- obs_fun$clonal_germline_v_gene
            tempo.secondary <- obs_fun$clonal_germline_j_gene
        } else if (col_fun == "c_allele") {
            tempo.primary <- obs_fun$c_call
            tempo.secondary <- NULL
        } else if (col_fun == "c_gene") {
            tempo.primary <- obs_fun$c_gene
            tempo.secondary <- NULL
        }
    }else{
        if (col_fun == "vj_allele") {
            tempo.primary <- obs_fun$v_call
            tempo.secondary <- obs_fun$j_call
        } else if (col_fun == "vj_gene") {
            tempo.primary <- obs_fun$v_gene
            tempo.secondary <- obs_fun$j_gene
        } else if (col_fun == "c_allele") {
            tempo.primary <- obs_fun$c_call
            tempo.secondary <- NULL
        } else if (col_fun == "c_gene") {
            tempo.primary <- obs_fun$c_gene
            tempo.secondary <- NULL
        }
    }

    # keep only the first allele or gene if igblast hesitated between several
    tempo.primary <- unlist(lapply(tempo.primary, function(x) {
        if (is.na(x)) return(NA) else return(strsplit(x, ",")[[1]][1])
    }))
    tempo.secondary <- unlist(lapply(tempo.secondary, function(x) {
        if (is.na(x)) return(NA) else return(strsplit(x, ",")[[1]][1])
    }))


    # # warning: I can use all.annotation.log because I have duplicated the first column of the dataset in the second column, in order to change the name in the first column with metadata. If all(obs_fun[ , 1] == obs_fun[ , 2]) == TRUE, it means no annotation added
    if(grepl(x = names(obs_fun)[2], pattern = "^initial_")){ # means that fonctional annotations are present
        annotation.log <- obs_fun[ , 1] == obs_fun[ , 2]
        all.annotation.log <- all(annotation.log) # if one difference between obs_fun[ , 1] == obs_fun[ , 2], then all.annotation.log is FALSE
    }else{
        all.annotation.log <- TRUE
    }
    if(all.annotation.log == FALSE){
        tempo.labels <- obs_fun[ , 1]
        tempo.labels[annotation.log] <- NA
    }

    #inactivated because now chain can be both IGL and IGK
    # if(length(chain) != 1){
    # stop(paste0("\n\n============\n\nINTERNAL CODE ERROR 4 IN donut.R for kind ", kind_fun, ": chain MUST BE A SINGLE VALUE.\nHERE IT IS: ", paste(chain, collapse = " "), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
    # }else if(chain != unique(substr(tempo.v, 1, 3))){
    # stop(paste0("\n\n============\n\nINTERNAL CODE ERROR 5 IN donut.R for kind ", kind_fun, ": CHAIN OF J DIFFERENT FROM CHAIN OF V.\nCHAIN OF V: ", paste(chain, collapse = " "), "\nCHAIN OF J: ", paste(unique(substr(tempo.v, 1, 3)), collapse = " "), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
    # }

    if (!is.null(tempo.secondary) && !all(is.na(tempo.secondary)) && !all(is.na(tempo.primary))) { # if there is a secondary column (meaning we are looking at VJ alleles/genes) check that it is the same chain
        check1 <- substr(tempo.primary, 1, 3)
        check2 <- substr(tempo.secondary, 1, 3)
        if( ! all(check1 == check2)){
            stop(paste0("\n\n============\n\nINTERNAL CODE ERROR 5 IN donut.R for kind ", kind_fun, ": CHAIN OF J DIFFERENT FROM CHAIN OF V.\nCHAIN OF V: ", paste(tempo.primary[check1 != check2], collapse = " "), "\nCHAIN OF J: ", paste(tempo.secondary[check1 != check2], collapse = " "), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
        }
    }

    # remove "IGH" if heavy chain ; remove "IG" if light chain to keep the info about K or L locus (if all values are NA it is correctly handled)
    third_character <- substr(tempo.primary, 3, 3)
    third_character_non_na <- third_character[!is.na(third_character)]
    if(length(third_character_non_na) == 0 || !all(third_character_non_na == "H")){
        tronc <- 3
    } else {
        tronc <- 4
    }
    tempo.primary <- ifelse(is.na(tempo.primary), NA, substring(tempo.primary, tronc))
    if (!is.null(tempo.secondary)) {
        tempo.secondary <- ifelse(is.na(tempo.secondary), NA, substring(tempo.secondary, tronc))
        clone.name <- paste0(tempo.primary, "_", tempo.secondary)
    } else {
        clone.name <- tempo.primary
    }

    # nb of na
    tempo_primary_na <- is.na(tempo.primary)
    if(is.null(tempo.secondary)){
        nb_na_log <- tempo_primary_na
    }else{
        tempo_secondary_na <- is.na(tempo.secondary)
        nb_na_log <- tempo_primary_na | tempo_secondary_na
    }
    if(kind_fun == "annotated" & all.annotation.log != TRUE){
        tempo_log <- clone.name %in% clone.name[ ! is.na(tempo.labels)]
        nb_na <- sum(nb_na_log[tempo_log])
        nb_total <- sum(tempo_log, na.rm = TRUE)
    }else{
        nb_na <- sum(nb_na_log)
        nb_total <- length(tempo.primary)
    }
    # end nb of na

    obs2 <- data.frame(table(clone.name))
    if (nrow(obs2) > 0) {
        names(obs2)[1] <- col_fun
        obs2 <- data.frame(obs2, Prop = obs2$Freq / sum(obs2$Freq), kind = kind_fun)
        obs2[[col_fun]] <- factor(obs2[[col_fun]], levels = obs2[[col_fun]][order(obs2$Prop, decreasing = TRUE)]) # reorder so that the donut is according to decreasing proportion starting at the top in a clockwise direction
        obs2 <- obs2[order(as.numeric(obs2[[col_fun]]), decreasing = FALSE), ] # obs2 with rows in decreasing order, according to Prop

        if(all.annotation.log == FALSE){
            tempo.data1 <- aggregate(tempo.labels ~ clone.name, FUN = function(x){paste(x, collapse = ",")})
            obs2 <- data.frame(obs2, labels = tempo.data1$tempo.labels[match(obs2[[col_fun]], tempo.data1$clone.name)])
        }

        if(kind_fun == "annotated" & all.annotation.log != TRUE){
            obs3 <- obs2[ ! is.na(obs2$labels) ,]
        }else{
            obs3 <- obs2
        }

    } else { # Handle the situation in which all values are NA in the original column
        # structure vide avec les bonnes colonnes pour éviter les erreurs en aval
        obs2 <- data.frame(
            temp_col = character(0),
            Freq = integer(0),
            Prop = numeric(0),
            kind = character(0)
        )
        names(obs2)[1] <- col_fun
        obs2[[col_fun]] <- factor(obs2[[col_fun]], levels = character(0))  # nécessaire pour éviter erreur au reorder
        obs3 <- obs2
    }
    return(list(obs2 = obs2, obs3 = obs3, all.annotation.log = all.annotation.log, nb_na = nb_na, nb_total = nb_total))
}


# Classes of a donut displayed in the legend: same selection as the legend.limit argument of gg_donut() (cute_little_R_functions file)
# freq_fun: frequencies of the classes (Freq column of obs3)
# class_fun: names of the classes (first column of obs3)
# limit_fun: donut_legend_limit, i.e., NULL (ten biggest classes), "all" (all the classes) or a proportion (classes with a proportion above it, or the ten biggest classes if none)
# Returns the names of the classes in the legend
legend_classes <- function(freq_fun, class_fun, limit_fun){
    tempo.prop <- freq_fun / sum(freq_fun)
    if(is.character(limit_fun)){ # "all"
        limit_fun <- 0
    }else if(is.numeric(limit_fun)){
        if(sum(tempo.prop >= limit_fun) == 0){ # limit higher than the biggest proportion
            limit_fun <- NULL
        }
    }
    if(is.null(limit_fun)){
        if(length(tempo.prop) >= 10){
            limit_fun <- sort(tempo.prop, decreasing = TRUE)[10] # ten biggest classes, and those with the same proportion as the tenth
        }else{
            limit_fun <- 0
        }
    }
    return(as.character(class_fun[tempo.prop >= limit_fun]))
}


# Colors of the classes that are in the legend of at least one donut
# n_fun: number of classes
# palette_fun: donut_palette, i.e., NULL (default ggplot2 colors) or the name of a ColorBrewer palette
# Returns n_fun colors, the first ones being for the biggest classes. No grey, used for the classes not in the legend
class_colors <- function(n_fun, palette_fun){
    if(n_fun == 0){
        return(character(0))
    }
    if(is.null(palette_fun)){
        tempo.colors <- scales::hue_pal()(n_fun) # default ggplot2 colors
    }else if(RColorBrewer::brewer.pal.info[palette_fun, "category"] == "qual"){
        tempo.colors <- RColorBrewer::brewer.pal(n = RColorBrewer::brewer.pal.info[palette_fun, "maxcolors"], name = palette_fun)
        tempo.rgb <- grDevices::col2rgb(tempo.colors)
        tempo.colors <- tempo.colors[apply(tempo.rgb, 2, max) > apply(tempo.rgb, 2, min)] # grey of the palette removed
        if(n_fun <= length(tempo.colors)){
            return(tempo.colors[1:n_fun]) # colors already made to be distinct, in this order
        }
        # more classes than colors: colors interpolated between the successive colors of the palette in the HCL space, along the shortest way around the hue circle (grDevices::colorRampPalette() gives greyish mixes of opposite colors, which look like the grey of the classes not in the legend)
        tempo.hcl <- farver::decode_colour(tempo.colors, to = "hcl") # columns h, c and l
        tempo.pos <- seq(1, nrow(tempo.hcl), length.out = n_fun) # position of each color among the colors of the palette
        tempo.i <- pmin(floor(tempo.pos), nrow(tempo.hcl) - 1) # previous color of the palette
        tempo.f <- tempo.pos - tempo.i # fraction of the way to the next color of the palette
        tempo.dh <- ((tempo.hcl[tempo.i + 1, "h"] - tempo.hcl[tempo.i, "h"] + 180) %% 360) - 180 # shortest hue difference
        tempo.colors <- farver::encode_colour(cbind(
            h = (tempo.hcl[tempo.i, "h"] + tempo.f * tempo.dh) %% 360,
            c = tempo.hcl[tempo.i, "c"] + tempo.f * (tempo.hcl[tempo.i + 1, "c"] - tempo.hcl[tempo.i, "c"]),
            l = tempo.hcl[tempo.i, "l"] + tempo.f * (tempo.hcl[tempo.i + 1, "l"] - tempo.hcl[tempo.i, "l"])
        ), from = "hcl")
    }else if(n_fun <= RColorBrewer::brewer.pal.info[palette_fun, "maxcolors"]){ # sequential or diverging palette
        tempo.colors <- scales::brewer_pal(palette = palette_fun)(n_fun) # as ggplot2::scale_fill_brewer()
    }else{
        tempo.colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = RColorBrewer::brewer.pal.info[palette_fun, "maxcolors"], name = palette_fun))(n_fun)
    }
    # colors taken along the scale by steps of the golden ratio: two classes of successive sizes get distant colors, whatever the number of classes
    tempo.order <- rank((0:(n_fun - 1) * (sqrt(5) - 1) / 2) %% 1, ties.method = "first")
    return(tempo.colors[tempo.order])
}


################################ End Functions


################################ Pre-ignition checking


# reserved words
# end reserved words
# argument primary checking
arg.check <- NULL #
text.check <- NULL #
checked.arg.names <- NULL # for function debbuging: used by r_debugging_tools
ee <- expression(arg.check <- c(arg.check, tempo$problem) , text.check <- c(text.check, tempo$text) , checked.arg.names <- c(checked.arg.names, tempo$object.name))
for(i0 in tempo.arg.names){
    tempo <- fun_check(data = get(i0), class = "vector", typeof = "character", length = 1) ; eval(ee)
}
if(any(arg.check) == TRUE){ # normally no NA
    stop(paste0("\n\n================\n\n", paste(text.check[arg.check], collapse = "\n"), "\n\n================\n\n"), call. = FALSE) # == in stop(), not in tempo.cat, to be able to add several messages between == #
}
# end argument primary checking
# second round of checking and data preparation
# management of NA arguments
# end management of NA arguments
# management of NULL arguments, WARNING: only for donut.R because NULL is "NULL" in the nextflow.config file
tempo.arg <-c(
    "file_name", 
    "kind", 
    "col",
    # "donut_palette", # can be NULL
    "donut_hole_size",
    "donut_hole_text",
    "donut_hole_text_size",
    "donut_border_color",
    "donut_border_size",
    "donut_annotation_distance",
    "donut_annotation_size",
    "donut_annotation_force",
    "donut_annotation_force_pull",
    "donut_legend_width",
    "donut_legend_text_size", 
    "donut_legend_box_size", 
    "donut_legend_box_space", 
    # "donut_legend_limit", # can be NULL
    "cute",
    "repertoire_names_ch",
    "color_ref_files",
    "log"
)
tempo.log <- sapply(lapply(tempo.arg, FUN = get, env = sys.nframe(), inherit = FALSE), FUN = function(x){x == "NULL"}) # WARNING: only for donut.R
if(any(tempo.log) == TRUE){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN donut.R:\n", ifelse(sum(tempo.log, na.rm = TRUE) > 1, "THESE ARGUMENTS\n", "THIS ARGUMENT\n"), paste0(tempo.arg[tempo.log], collapse = "\n"),"\nCANNOT BE NULL")
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop(), not in tempo.cat, to be able to add several messages between ==
}
# end management of NULL arguments, WARNING: only for donut.R because NULL is "NULL" in the nextflow.config file
# seed
set.seed(1)
# end seed
# warning initiation
ini.warning.length <- options()$warning.length
options(warning.length = 8170)
# warn.count <- 0 # not required
# end warning initiation
# other checkings (not full checked because already checked in the .nf file)
if( ! file.exists(file_name)){
    tempo.cat <- paste0("ERROR IN donut.R:\nTHE file_name PARAMETER MUST BE A VALID PATH OF A FILE IS NOT \"NULL\"\nHERE IT IS: \n", paste0(file_name, collapse = " "))
    stop(paste0("\n\n================\n\n", tempo$text, "\n\n================\n\n"), call. = FALSE)
}

tempo <- fun_check(data = kind, options = c("all", "clone", "annotated"), length = 1) ; eval(ee)
if(tempo$problem == TRUE){
    stop(paste0("\n\n================\n\n", tempo$text, "\n\n================\n\n"), call. = FALSE)
}
tempo <- fun_check(data = col, options = c("vj_allele", "c_allele", "vj_gene", "c_gene"), length = 1) ; eval(ee)
if(tempo$problem == TRUE){
    stop(paste0("\n\n================\n\n", tempo$text, "\n\n================\n\n"), call. = FALSE)
}

# following parameter are those of the gg_donut() function and are checked by this one
if(donut_palette == "NULL"){
    donut_palette <- NULL
}else if( ! donut_palette %in% rownames(RColorBrewer::brewer.pal.info)){ # also checked by gg_donut(), but used before by class_colors()
    tempo.cat <- paste0("ERROR IN donut.R:\nTHE donut_palette PARAMETER MUST BE \"NULL\" OR THE NAME OF A COLORBREWER PALETTE (SEE https://ggplot2.tidyverse.org/reference/scale_brewer.html#palettes)\nHERE IT IS: \n", paste0(donut_palette, collapse = " "))
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE)
}
donut_hole_size <- as.numeric(donut_hole_size) # numeric string already checked by nextflow
if( ! (length(donut_hole_text) == 1 & any(donut_hole_text %in% c("TRUE", "FALSE")))){ # positive numeric
    tempo.cat <- paste0("ERROR IN donut.R:\nTHE donut_hole_text PARAMETER MUST BE \"TRUE\" OR \"FALSE\"\nHERE IT IS: \n", paste0(donut_hole_text, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}else if(donut_hole_text == "TRUE"){
    donut_hole_text <- TRUE
}else{
    donut_hole_text <- FALSE
}
donut_hole_text_size <- as.numeric(donut_hole_text_size) # numeric string already checked by nextflow
# nothing to check for donut_border_color
donut_border_size <- as.numeric(donut_border_size) # numeric string already checked by nextflow
donut_annotation_distance <- as.numeric(donut_annotation_distance) # numeric string already checked by nextflow
donut_annotation_size <- as.numeric(donut_annotation_size) # numeric string already checked by nextflow
donut_annotation_force <- as.numeric(donut_annotation_force) # numeric string already checked by nextflow
donut_annotation_force_pull <- as.numeric(donut_annotation_force_pull) # numeric string already checked by nextflow
donut_legend_width <- as.numeric(donut_legend_width) # numeric string already checked by nextflow
donut_legend_text_size <- as.numeric(donut_legend_text_size) # numeric string already checked by nextflow
donut_legend_box_size <- as.numeric(donut_legend_box_size) # numeric string already checked by nextflow
donut_legend_box_space <- as.numeric(donut_legend_box_space) # numeric string already checked by nextflow
if(donut_legend_limit == "NULL"){
    donut_legend_limit <- NULL
}else if(donut_legend_limit == "all"){
    donut_legend_limit <- "all" # kept as a character string: gg_donut() then puts all the sectors in the legend
}else{
    donut_legend_limit <- as.numeric(donut_legend_limit) # numeric string already checked by nextflow
}
# end following parameter are those of the gg_donut() function and are checked by this one

# other checkings (not full checked because already checked in the .nf file)
# reserved word checking
# end reserved word checking
# end second round of checking and data preparation
# package checking
# end package checking


################################ End pre-ignition checking


################################ Main code


################ Ignition


fun_report(data = paste0("\n\n################################################################ donut PROCESS\n\n"), output = log, path = "./", overwrite = FALSE)
ini.date <- Sys.time()
ini.time <- as.numeric(ini.date) # time of process begin, converted into seconds
fun_report(data = paste0("\n\n################################ RUNNING DATE AND STARTING TIME\n\n"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0(ini.date, "\n\n"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0("\n\n################################ RUNNING\n\n"), output = log, path = "./", overwrite = FALSE)


################ End ignition


################ Graphical parameter initialization


################ End graphical parameter initialization





################ Data import

obs <- read.table(file_name, sep = "\t", header = TRUE, comment.char = "")

tempo.cat <- paste0("\n\nREPERTOIRE NAMES : ", repertoire_names_ch, "\n\n")
fun_report(data = tempo.cat, output = log, path = "./", overwrite = FALSE)
cat(tempo.cat)

################ End Data import

################ data modification

type_text <- switch(
    col,
    "vj_allele" = "V and J alleles",
    "c_allele" = "C alleles",
    "vj_gene" = "V and J genes",
    "c_gene" = "C genes",
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR: unexpected allele_type value: ", col, "\n\n================\n\n"))
)

tempo.title <- paste0(
    ifelse(
        kind == "all", 
        paste0(
            "Donut plot of the wanted sequences (see the corresponding wanted_seq.tsv) grouped by same ", type_text, "\n\n",
            "Kind of sequences: ",
            kind
        ), 
        ifelse(
            kind == "clone", 
            paste0(
                "Donut plot of the sequences in clones (see the corresponding clone_assigned_seq.tsv), grouped by same ", type_text, "\n\n",
                "Kind of sequences: ",
                kind
            ), 
            ifelse(
                kind == "annotated", 
                paste0(
                    "Donut plot of the wanted sequences (see the corresponding wanted_seq.tsv) grouped by same ", type_text, ",\nfor which at least one name replacement is present (according to the meta_name_replacement parameter of the nextflow.config file)\n\n",
                    "Kind of sequences: ",
                    kind
                ), 
                stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 4 IN donut.R for kind THAT CAN ONLY BE \"all\", \"annotated\" OR \"clone\".\nHER IT IS: ", kind, "\n\n================\n\n"))
            )
        )
    )
)

# Chain is determined with igblast_ref_files because if we studied IGL and IGK but only one was found, both still need to be in the title
loci <- sapply(repertoire_names_ch, extract_chain)
chain <- unique(loci[!is.na(loci)]) # extract the IGH or IGK name (ignoring NA values)

if(nrow(obs) > 0){
    tempo.classes <- donut_classes(obs_fun = obs, kind_fun = kind, col_fun = col) # same function as for the color reference below
    obs2 <- tempo.classes$obs2 # classes of all the sequences
    obs3 <- tempo.classes$obs3 # classes of the donut
    all.annotation.log <- tempo.classes$all.annotation.log
    nb_na <- tempo.classes$nb_na
    nb_total <- tempo.classes$nb_total

    # classes in the legend of the donut (the other sectors are in grey). NULL if no donut plotted
    if(nrow(obs3) > 0 & ! (kind == "annotated" & all.annotation.log == TRUE)){
        legend.class <- legend_classes(freq_fun = obs3$Freq, class_fun = obs3[[col]], limit_fun = donut_legend_limit)
    }else{
        legend.class <- NULL
    }

    # colors of the classes, the same in all the donuts of the same col (kinds "all", "annotated" and "clone"): an allele or gene keeps its color from one donut to the other
    # each class in the legend of at least one of these donuts gets a color, in the order of the legends of the "all", then "annotated", then "clone" donuts
    # the legends of all the kinds are computed from the color_ref_files files (wanted_seq.tsv for "all" and "annotated", clone_assigned_seq.tsv for "clone" if present), with the code used for the plotted donut
    legend.classes <- list(all = NULL, annotated = NULL, clone = NULL) # classes in the legend of each kind of donut
    for(i1 in strsplit(color_ref_files, split = " ")[[1]]){
        if(basename(i1) == "wanted_seq.tsv"){
            tempo.kinds <- c("all", "annotated")
        }else if(basename(i1) == "clone_assigned_seq.tsv"){
            tempo.kinds <- "clone"
        }else{
            stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 6 IN donut.R: THE color_ref_files PARAMETER CAN ONLY CONTAIN wanted_seq.tsv AND clone_assigned_seq.tsv FILES.\nHERE IT IS: ", color_ref_files, "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
        }
        tempo.ref <- read.table(i1, sep = "\t", header = TRUE, comment.char = "")
        if(nrow(tempo.ref) > 0){
            for(i2 in tempo.kinds){
                tempo.ref.classes <- donut_classes(obs_fun = tempo.ref, kind_fun = i2, col_fun = col)
                if(nrow(tempo.ref.classes$obs3) > 0 & ! (i2 == "annotated" & tempo.ref.classes$all.annotation.log == TRUE)){ # donut plotted
                    legend.classes[[i2]] <- legend_classes(freq_fun = tempo.ref.classes$obs3$Freq, class_fun = tempo.ref.classes$obs3[[col]], limit_fun = donut_legend_limit)
                }
            }
        }
    }
    color.classes <- unique(c(legend.classes$all, legend.classes$annotated, legend.classes$clone))
    class.colors <- class_colors(n_fun = length(color.classes), palette_fun = donut_palette)
    names(class.colors) <- color.classes
    if( ! all(legend.class %in% color.classes)){
        stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 7 IN donut.R for kind ", kind, ": CLASSES OF THE LEGEND WITHOUT COLOR (ABSENT FROM THE color_ref_files FILES):\n", paste(legend.class[ ! legend.class %in% color.classes], collapse = " "), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
    }
    # colors of the sectors of the donut, in the order of the rows of obs3 (decreasing frequencies), as required by the fill.color argument of gg_donut()
    sector.colors <- rep("grey80", nrow(obs3)) # grey: classes not in the legend
    tempo.log <- as.character(obs3[[col]]) %in% legend.class
    sector.colors[tempo.log] <- class.colors[as.character(obs3[[col]])[tempo.log]]
    tempo.cat <- paste0("\n\nCOLORS OF THE CLASSES IN THE LEGEND OF AT LEAST ONE ", col, " DONUT (SAME COLOR IN ALL THESE DONUTS):\n", paste0(color.classes, ": ", class.colors, collapse = "\n"), "\n\n")
    fun_report(data = tempo.cat, output = log, path = "./", overwrite = FALSE)


    ################ End data modification
    
    ################ Plotting
    
    tempo.title <- paste0(
        tempo.title,
        "\nChain: ", paste(chain, collapse = " "), "\n",
        "Proportions of sequences with the indicated ", type_text, " are indicated in the legend, after the labels.\n", 
        "The total number of non NA ", type_text, " detected is indicated in the center of the donut.\n",
        "See the donut_stats.tsv file for the complete stats.\n", 
        "Number of sequences selected (depending on \"all\", \"annotated\" or \"clone\"): ",
        nb_total,
        ".\nNumber of NA removed from the ", type_text, " among the selected number of sequences: ",
        nb_na,
        "."
    )
    if( ! is.null(legend.class)){ # donut plotted
        if(length(legend.class) < nrow(obs3)){ # some classes not in the legend (never with donut_legend_limit "all")
            if(is.null(donut_legend_limit)){
                tempo.limit <- "NULL: ten biggest classes"
            }else if(any(obs3$Freq / sum(obs3$Freq) >= donut_legend_limit)){
                tempo.limit <- round(donut_legend_limit, 2)
            }else{
                tempo.limit <- paste0(round(donut_legend_limit, 2), ", above all the proportions: ten biggest classes")
            }
            tempo.title <- paste0(tempo.title, "\nClasses not in the legend are in grey (donut_legend_limit = ", tempo.limit, ").")
        }
    }
    
    if(length(obs3$Freq) != 0){
        if( ! (all(obs3$Freq == 0) | all(is.na(obs3$Freq)))){
            pdf(NULL) # used because I need plot = TRUE for return.gtable
            tempo.plot <- gg_donut(
                data1 = obs3, # select only the classes with annotations if kind == "annotated" & all.annotation.log != TRUE
                freq = "Freq", 
                categ = col, 
                fill.palette = donut_palette, # not used: overridden by fill.color
                fill.color = sector.colors, # same color for a class in all the donuts, grey for the classes not in the legend
                hole.size = donut_hole_size, 
                hole.text = donut_hole_text, 
                hole.text.size = donut_hole_text_size, 
                border.color = donut_border_color, 
                border.size = donut_border_size, 
                title = tempo.title, 
                title.text.size = 5, 
                annotation = if(all.annotation.log == FALSE){"labels"}else{NULL},
                annotation.distance = donut_annotation_distance,
                annotation.size = donut_annotation_size,
                annotation.force = donut_annotation_force,
                annotation.force.pull = donut_annotation_force_pull,
                legend.show = TRUE, 
                legend.width = donut_legend_width, 
                legend.name = type_text, 
                legend.text.size = donut_legend_text_size, 
                legend.box.size = donut_legend_box_size, 
                legend.box.space = donut_legend_box_space, 
                legend.limit = donut_legend_limit, 
                legend.add.prop = TRUE,
                add = NULL, 
                return = TRUE, 
                return.ggplot = FALSE,
                return.gtable = TRUE,
                plot = FALSE, 
                warn.print = TRUE, 
                lib_path = NULL
            )
            final.plot <- suppressMessages(suppressWarnings(gridExtra::grid.arrange(tempo.plot$gtable))) # , left = " ", right = " " : trick to add margins in the plot. padding =  unit(0.5, "inch") is for top margin above the title
            dev.off()
        }else{
            # no need to use pdf(NULL) with fun_gg_empty_graph()
            final.plot <- fun_gg_empty_graph(text = "NO DONUT PLOTTED\nNO FREQUENCY DETECTED", text.size = 3, title = tempo.title, title.size = 4)
        }
    }else{
        # no need to use pdf(NULL) with fun_gg_empty_graph()
        final.plot <- fun_gg_empty_graph(text = "NO DONUT PLOTTED\nNO FREQUENCY DETECTED", text.size = 3, title = tempo.title, title.size = 4)
    }
    if(kind == "annotated" & all.annotation.log == TRUE){ # overwrite final.plot
        # no need to use pdf(NULL) with fun_gg_empty_graph()
        final.plot <- fun_gg_empty_graph(text = "NO DONUT PLOTTED\nNO SEQUENCE NAME REPLACEMENT DETECTED\n(meta_name_replacement PARAMETER OF THE nextflow.config FILE)", text.size = 3, title = tempo.title, title.size = 4)
    }
}else{
    final.plot <- fun_gg_empty_graph(text = "NO DONUT PLOTTED\nNO SEQUENCE DETECTED", text.size = 3, title = tempo.title, title.size = 4)
    obs3 <- obs
}




################ end Plotting

################ saving plots

ggplot2::ggsave(filename = paste0(col, "_", kind, "_donutchart.png"), plot = final.plot, device = "png", path = ".", width = 120, height = 120, units = "mm", dpi = 300, bg = "white")
ggplot2::ggsave(filename = paste0(col, "_", kind, "_donutchart.svg"), plot = final.plot, device = "svg", path = ".", width = 120, height = 120, units = "mm", dpi = 300, bg = "white")
ggplot2::ggsave(filename = paste0(col, "_", kind, "_donutchart.pdf"), plot = final.plot, device = "pdf", path = ".", width = 120, height = 120, units = "mm", dpi = 300, bg = "white")


################ end saving plots

################ data saving

write.table(obs3, file = paste0("./", kind, "_donut_stats.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)

################ End data saving

################ Pdf window closing


graphics.off()


################ end Pdf window closing


################ Seeding inactivation


set.seed(NULL)


################ end Seeding inactivation


################ Environment saving


save(list = ls(), file = "./all_objects.RData")
fun_report(data = paste0("\n\n################################ RUNNING END"), output = log, path = "./", overwrite = FALSE)
end.date <- Sys.time()
end.time <- as.numeric(end.date)
total.lapse <- round(lubridate::seconds_to_period(end.time - ini.time))
fun_report(data = paste0("\n\nEND TIME: ", end.date), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0("\n\nTOTAL TIME LAPSE: ", total.lapse), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0("\n\nALL DATA SAVED IN all_objects.RData"), output = log, path = "./", overwrite = FALSE)


################ end Environment saving


################ Warning messages


on.exit(exp = options(warning.length = ini.warning.length), add = TRUE)


################ end Warning messages


################ Parameter printing


fun_report(data = paste0("\n\n################################ INITIAL SETTINGS OF PARAMETERS"), output = log, path = "./", overwrite = FALSE)
fun_report(data = param.ini.settings, output = log, path = "./", overwrite = FALSE, , vector.cat = TRUE)
fun_report(data = paste0("\n\n################################ R SYSTEM AND PACKAGES"), output = log, path = "./", overwrite = FALSE)
tempo <- sessionInfo()
tempo$otherPkgs <- tempo$otherPkgs[order(names(tempo$otherPkgs))] # sort the packages
tempo$loadedOnly <- tempo$loadedOnly[order(names(tempo$loadedOnly))] # sort the packages
fun_report(data = tempo, output = log, path = "./", overwrite = FALSE, , vector.cat = TRUE)
fun_report(data = paste0("\n\n################################ JOB END\n\nTIME: ", end.date, "\n\nTOTAL TIME LAPSE: ", total.lapse, "\n"), output = log, path = "./", overwrite = FALSE)


################ end Parameter printing


################################ End Main code
