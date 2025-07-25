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
if(version$version.string != "R version 4.1.2 (2021-11-01)"){
    stop(paste0("\n\n================\n\nERROR IN donut.R\n\n\n", version$version.string, " IS NOT THE 4.1.2 RECOMMANDED\n\n================\n\n"))
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

# setwd("C:/Users/gael/Documents/Git_projects/19532_marbouty/dataset/test")
# file_name = "./caca.tsv"
# kind = "all"
# col = "vj_allele"
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
# cute = "https://gitlab.pasteur.fr/gmillot/cute_little_R_functions/-/raw/v12.2.0/cute_little_R_functions.R"
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
# fun_pack(req.package = req.package.list, load = TRUE, lib.path = NULL) # packages are imported even if inside functions are written as package.name::function() in the present code

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
    "donut_legend_width",
    "donut_legend_text_size", 
    "donut_legend_box_size", 
    # "donut_legend_limit", # can be NULL
    "cute",
    "repertoire_names_ch",
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
warn <- NULL
# warn.count <- 0 # not required
# end warning initiation
# other checkings (not full checked because already checked in the .nf file)
if( ! file.exists(file_name)){
    tempo.cat <- paste0("ERROR IN donut.R:\nTHE file_name PARAMETER MUST BE A VALID PATH OF A FILE IS NOT \"NULL\"\nHERE IT IS: \n", paste0(file_name, collapse = " "))
    stop(paste0("\n\n================\n\n", tempo$text, "\n\n================\n\n"), call. = FALSE)
}

tempo <- fun_check(data = kind, options = c("all", "tree", "annotated"), length = 1) ; eval(ee)
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
}
donut_hole_size <- as.numeric(donut_hole_size) # numeric string already checked by nextflow
if( ! (length(donut_hole_text) == 1 & any(donut_hole_text %in% c("TRUE", "FALSE")))){ # positive numeric
    tempo.cat <- paste0("ERROR IN donut.R:\nTHE tree_label_size PARAMETER MUST BE \"TRUE\" OR \"FALSE\"\nHERE IT IS: \n", paste0(donut_hole_text, collapse = " "))
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

tempo.cat <- paste0("\n\nREPERTOIRENAMES : ", repertoire_names_ch, "\n\n")
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
            "Donut plot of the productive sequences (see the corresponding productive_seq.tsv) grouped by same ", type_text, "\n\n",
            "Kind of sequences: ",
            kind
        ), 
        ifelse(
            kind == "tree", 
            paste0(
                "Donut plot of the clone-assigned sequences (see the corresponding germ_tree_seq.tsv), grouped by same ", type_text, "\n\n",
                "Kind of sequences: ",
                kind
            ), 
            ifelse(
                kind == "annotated", 
                paste0(
                    "Donut plot of the productive sequences (see the corresponding productive_seq.tsv) grouped by same ", type_text, ",\nfor which at least one name replacement is present (according to the meta_name_replacement parameter of the nextflow.config file)\n\n",
                    "Kind of sequences: ",
                    kind
                ), 
                stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 4 IN donut.R for kind THAT CAN ONLY BE \"all\", \"annotated\" OR \"tree\".\nHER IT IS: ", kind, "\n\n================\n\n"))
            )
        )
    )
)

# Chain is determined with igblast_ref_files because if we studied IGL and IGK but only one was found, both still need to be in the title
loci <- sapply(repertoire_names_ch, extract_chain)
chain <- unique(loci[!is.na(loci)]) # extract the IGH or IGK name (ignoring NA values)

if(nrow(obs) > 0){
    
    if(kind == "tree"){
        if (col == "vj_allele") {
            tempo.primary <- obs$germline_v_call
            tempo.secondary <- obs$germline_j_call
        } else if (col == "vj_gene") {
            tempo.primary <- obs$germline_v_gene
            tempo.secondary <- obs$germline_j_gene
        } else if (col == "c_allele") {
            tempo.primary <- obs$c_call
            tempo.secondary <- NULL
        } else if (col == "c_gene") {
            tempo.primary <- obs$c_gene
            tempo.secondary <- NULL
        }
    }else{
        if (col == "vj_allele") {
            tempo.primary <- obs$v_call
            tempo.secondary <- obs$j_call
        } else if (col == "vj_gene") {
            tempo.primary <- obs$v_gene
            tempo.secondary <- obs$j_gene
        } else if (col == "c_allele") {
            tempo.primary <- obs$c_call
            tempo.secondary <- NULL
        } else if (col == "c_gene") {
            tempo.primary <- obs$c_gene
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
    
    #inactivated because now chain can be both IGL and IGK
    # if(length(chain) != 1){
    # stop(paste0("\n\n============\n\nINTERNAL CODE ERROR 4 IN donut.R for kind ", kind, ": chain MUST BE A SINGLE VALUE.\nHERE IT IS: ", paste(chain, collapse = " "), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE) 
    # }else if(chain != unique(substr(tempo.v, 1, 3))){
    # stop(paste0("\n\n============\n\nINTERNAL CODE ERROR 5 IN donut.R for kind ", kind, ": CHAIN OF J DIFFERENT FROM CHAIN OF V.\nCHAIN OF V: ", paste(chain, collapse = " "), "\nCHAIN OF J: ", paste(unique(substr(tempo.v, 1, 3)), collapse = " "), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE) 
    # }
    
    if (!is.null(tempo.secondary) && !all(is.na(tempo.secondary)) && !all(is.na(tempo.primary))) { # if there is a secondary column (meaning we are looking at VJ alleles/genes) check that it is the same chain
        check1 <- substr(tempo.primary, 1, 3)
        check2 <- substr(tempo.secondary, 1, 3)
        if( ! all(check1 == check2)){
            stop(paste0("\n\n============\n\nINTERNAL CODE ERROR 5 IN donut.R for kind ", kind, ": CHAIN OF J DIFFERENT FROM CHAIN OF V.\nCHAIN OF V: ", paste(tempo.v[check1 != check2], collapse = " "), "\nCHAIN OF J: ", paste(tempo.j[check1 != check2], collapse = " "), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE) 
        }
    }
    
    # remove "IGH" if heavy chain ; remove "IG" if light chain to keep the info about K or L locus (if all values are NA it is correctly handled)
    third_character <- substr(tempo.primary, 3, 3)
    third_character_non_na <- third_character[!is.na(third_character)]
    if(length(third_character_non_na) == 0 || !all(third_character_non_na == "H")) {
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
    
    
    obs2 <- data.frame(table(clone.name))
    if (nrow(obs2) > 0) {
        names(obs2)[1] <- col
        obs2 <- data.frame(obs2, Prop = obs2$Freq / sum(obs2$Freq), kind = kind)
        obs2[[col]] <- factor(obs2[[col]], levels = obs2[[col]][order(obs2$Prop, decreasing = TRUE)]) # reorder so that the donut is according to decreasing proportion starting at the top in a clockwise direction
        obs2 <- obs2[order(as.numeric(obs2[[col]]), decreasing = FALSE), ] # obs2 with rows in decreasing order, according to Prop
    } else { # Handle the situation in which all values are NA in the original column
        # structure vide avec les bonnes colonnes pour éviter les erreurs en aval
        obs2 <- data.frame(
            temp_col = character(0),
            Freq = integer(0),
            Prop = numeric(0),
            kind = character(0)
        )
        names(obs2)[1] <- col
        obs2[[col]] <- factor(obs2[[col]], levels = character(0))  # nécessaire pour éviter erreur au reorder
    }
    # # warning: I can use all.annotation.log because I have duplicated the first column of the dataset in the second column, in order to change the name in the first column with metadata. If all(obs[ , 1] == obs[ , 2]) == TRUE, it means no annotation added
    if(grepl(x = names(obs)[2], pattern = "^initial_")){ # means that fonctional annotations are present
        annotation.log <- obs[ , 1] == obs[ , 2]
        all.annotation.log <- all(annotation.log) # if one difference between obs[ , 1] == obs[ , 2], then all.annotation.log is FALSE
    }else{
        all.annotation.log <- TRUE
    }
    if(all.annotation.log == FALSE){
        tempo.labels <- obs[ , 1]
        tempo.labels[annotation.log] <- NA
        tempo.data1 <- aggregate(tempo.labels ~ clone.name, FUN = function(x){paste(x, collapse = ",")})
        obs2 <- data.frame(obs2, labels = tempo.data1$tempo.labels[match(obs2[[col]], tempo.data1$clone.name)])
    }
    
    
    ################ End data modification
    
    ################ Plotting
    
    tempo.title <- paste0(
        tempo.title,
        "\nChain: ", paste(chain, collapse = " "), "\n",
        "Proportions of sequences with the indicated ", type_text, " are indicated in the legend, after the labels.\n", 
        "The total number of sequences of the donut is indicated in the center.\nSee the donut_stats.tsv file for the complete stats."
    )
    if( ! is.null(donut_legend_limit)){
        if(sum(obs2$Prop >= donut_legend_limit) < length(obs2$Prop)){
            tempo.title <- paste0(tempo.title, "\nSome classes have been removed from the legend (donut_legend_limit = ", round(donut_legend_limit, 2), "). See the donut_stats.tsv file")
        }
    }
    
    if(kind == "annotated" & all.annotation.log != TRUE){
        obs3 <- obs2[ ! is.na(obs2$labels) ,]
    }else{
        obs3 <- obs2
    }
    
    if(length(obs3$Freq) != 0){
        if( ! (all(obs3$Freq == 0) | all(is.na(obs3$Freq)))){
            pdf(NULL) # used because I need plot = TRUE for return.gtable
            tempo.plot <- fun_gg_donut(
                data1 = obs3, # select only the classes with annotations if kind == "annotated" & all.annotation.log != TRUE
                freq = "Freq", 
                categ = col, 
                fill.palette = donut_palette,
                fill.color = NULL, 
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
                lib.path = NULL
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




################ end Plotting tree

################ saving plots

ggplot2::ggsave(filename = paste0(col, "_", kind, "_donutchart.png"), plot = final.plot, device = "png", path = ".", width = 120, height = 120, units = "mm", dpi = 300)
ggplot2::ggsave(filename = paste0(col, "_", kind, "_donutchart.svg"), plot = final.plot, device = "svg", path = ".", width = 120, height = 120, units = "mm", dpi = 300)
ggplot2::ggsave(filename = paste0(col, "_", kind, "_donutchart.pdf"), plot = final.plot, device = "pdf", path = ".", width = 120, height = 120, units = "mm", dpi = 300)


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


fun_report(data = paste0("\n\n################################ RECAPITULATION OF WARNING MESSAGES"), output = log, path = "./", overwrite = FALSE)
if( ! is.null(warn)){
    tempo.cat <- paste0("IN donut.R OF THE NEXFLOW EXECUTION:\n\n", warn)
    fun_report(data = tempo.cat, output = log, path = "./", overwrite = FALSE)
    cat(tempo.cat)
}else{
    fun_report(data = paste0("\n\nNO WARNING MESSAGE TO REPORT"), output = log, path = "./", overwrite = FALSE)
}
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
