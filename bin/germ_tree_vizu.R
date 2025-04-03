#!/usr/bin/Rscript
#########################################################################
##                                                                     ##
##     germ_tree_vizu.R                                                ##
##                                                                     ##
##     Gael A. Millot                                                  ##
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
    stop(paste0("\n\n================\n\nERROR IN germ_tree_vizu.R\n", version$version.string, " IS NOT THE 4.1.2 RECOMMANDED\n\n================\n\n"))
}
# other initializations
erase.objects = TRUE # write TRUE to erase all the existing objects in R before starting the algorithm and FALSE otherwise. Beginners should use TRUE
if(erase.objects == TRUE){
    rm(list = ls(all.names = TRUE))
    erase.objects = TRUE
}
erase.graphs = TRUE # write TRUE to erase all the graphic windows in R before starting the algorithm and FALSE otherwise
script <- "germ_tree_vizu"


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
        stop(paste0("\n\n================\n\nERROR IN germ_tree_vizu.R\nTHE args OBJECT HAS NA\n\n================\n\n"), call. = FALSE)
    }
    tempo.arg.names <- c(
        "germ_tree_kind", 
        "clone_nb_seq", 
        "germ_tree_duplicate_seq", 
        "germ_tree_leaf_color", 
        "germ_tree_leaf_shape", 
        "germ_tree_leaf_size", 
        "germ_tree_label_size", 
        "germ_tree_label_hjust", 
        "germ_tree_label_rigth", 
        "germ_tree_label_outside", 
        "germ_tree_right_margin", 
        "germ_tree_legend", 
        "germ_tree_data_assembly_path", 
        "germ_tree_meta_path", 
        "germ_tree_meta_legend", 
        "cute", 
        "log"
    ) # objects names exactly in the same order as in the bash code and recovered in args. Here only one, because only the path of the config file to indicate after the germ_tree_vizu.R script execution
    if(length(args) != length(tempo.arg.names)){
        stop(paste0("\n\n================\n\nERROR IN germ_tree_vizu.R\nTHE NUMBER OF ELEMENTS IN args (", length(args),") IS DIFFERENT FROM THE NUMBER OF ELEMENTS IN tempo.arg.names (", length(tempo.arg.names),")\nargs:", paste0(args, collapse = ","), "\ntempo.arg.names:", paste0(tempo.arg.names, collapse = ","), "\n\n================\n\n"), call. = FALSE)
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

# setwd("C:/Users/gmillot/Documents/Git_projects/repertoire_profiler/dev/test")
# setwd("C:/Users/gmillot/Documents/Git_projects/repertoire_profiler/work/e6/9a4aa4f90243b1a9c9308eb9952edb")
# setwd("Z:/Alice/sort1_VH/repertoire_profiler-v8.1/work/6a/5804f66785412c808f6f7d9ca46cca")
# setwd("Z:/Alice/human_Patient_B9_886_VH/ig_clustering-v8.8/work/ea/d1287c99e6a36d2790910c98fffcc0")
# setwd("C:/Users/gmillot/Documents/Git_projects/repertoire_profiler/work/e6/9a4aa4f90243b1a9c9308eb9952edb")
# germ_tree_kind = "rectangular"
# clone_nb_seq = "3"
# germ_tree_duplicate_seq = "FALSE" 
# germ_tree_leaf_color = "NULL" 
# germ_tree_leaf_shape = "23" 
# germ_tree_leaf_size = "3" 
# germ_tree_label_size = "2" 
# germ_tree_label_hjust = "-0.25" 
# germ_tree_label_rigth = "FALSE" 
# germ_tree_label_outside = "TRUE"
# germ_tree_right_margin = "1.5" 
# germ_tree_legend = "TRUE" 
# germ_tree_data_assembly_path = "caca"
# germ_tree_meta_path = "./metadata3.tsv" 
# germ_tree_meta_legend = "KD"
# cute = "https://gitlab.pasteur.fr/gmillot/cute_little_R_functions/-/raw/v12.8/cute_little_R_functions.R"
# log = "germ_tree_vizu.log"
# file.remove(c("./all_objects.RData", "./all_germ_trees.RData", "./germ_tree.pdf", "./germ_tree_vizu.log"))





################################ end Test

################################ Recording of the initial parameters


param.list <- c(
    "erase.objects", 
    "erase.graphs", 
    "script", 
    "run.way",
    "tempo.arg.names", 
    if(run.way == "SCRIPT"){"command"}, 
    "germ_tree_kind", 
    "clone_nb_seq", 
    "germ_tree_duplicate_seq", 
    "germ_tree_leaf_color", 
    "germ_tree_leaf_shape", 
    "germ_tree_leaf_size", 
    "germ_tree_label_size", 
    "germ_tree_label_hjust", 
    "germ_tree_label_rigth", 
    "germ_tree_label_outside", 
    "germ_tree_right_margin", 
    "germ_tree_legend", 
    "germ_tree_data_assembly_path", 
    "germ_tree_meta_path", 
    "germ_tree_meta_legend", 
    "cute", 
    "log"
)
if(any(duplicated(param.list))){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 1 IN germ_tree_vizu.R\nTHE param.list OBJECT CONTAINS DUPLICATED ELEMENTS:\n", paste(param.list[duplicated(param.list)], collapse = " "), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE) # message for developers
}
if(erase.objects == TRUE){
    created.object.control <- ls()[ ! ls() %in% "param.list"]
    if( ! (all(created.object.control %in% param.list) & all(param.list %in% created.object.control))){
        stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 2 IN germ_tree_vizu.R\nINCONSISTENCIES BETWEEN THE ARGUMENTS USED AND THE PARAMETERS REQUIRED IN THE EXECUTABLE CODE FILE\nTHE ARGUMENTS NOT PRESENT IN THE EXECUTABLE FILE (germ_tree_vizu.R) ARE:\n", paste(created.object.control[ ! created.object.control %in% param.list], collapse = " "), "\nTHE PARAMETERS OF THE EXECUTABLE FILE (germ_tree_vizu.R) NOT PRESENT IN THE ARGUMENTS ARE:\n", paste(param.list[ ! param.list %in% created.object.control], collapse = " "), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE) # message for developers
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
    stop(paste0("\n\n============\n\nERROR IN germ_tree_vizu.R\ncute PARAMETER MUST BE LENGTH 1: ", paste(cute, collapse = " "), "\n\n============\n\n"), call. = FALSE)
}else if(grepl(x = cute, pattern = "^http")){
    tempo.try <- try(suppressWarnings(suppressMessages(source(cute, local = .GlobalEnv))), silent = TRUE)
    if(any(grepl(x = tempo.try, pattern = "^[Ee]rror"))){
        stop(paste0("\n\n============\n\nERROR IN germ_tree_vizu.R\nHTTP INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else if( ! grepl(x = cute, pattern = "^http")){
    if( ! file.exists(cute)){
        stop(paste0("\n\n============\n\nERROR IN germ_tree_vizu.R\nFILE INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else{
    tempo.cat <- paste0("\n\n================\n\nINTERNAL CODE ERROR 3 IN germ_tree_vizu.R: CODE HAS TO BE MODIFIED\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n============\n\n")
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
    tempo.cat <- paste0("ERROR IN germ_tree_vizu.R\nREQUIRED cute FUNCTION", ifelse(length(tempo) > 1, "S ARE", " IS"), " MISSING IN THE R ENVIRONMENT:\n", paste0(tempo, collapse = "()\n"))
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end required function checking


################ local function: package import


# R Packages required
req.package.list <- c(
    "lubridate", 
    "tibble", 
    "svglite", 
    "ggplot2", 
    "ggtree",
    "dplyr",
    "dowser", 
    "qpdf"


)
for(i in 1:length(req.package.list)){suppressMessages(library(req.package.list[i], character.only = TRUE))}
# fun_pack(req.package = req.package.list, load = TRUE, lib.path = NULL) # packages are imported even if inside functions are written as package.name::function() in the present code


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
# management of NULL arguments, WARNING: only for germ_tree_vizu.R because NULL is "NULL" in the nextflow.config file
tempo.arg <-c(
    "germ_tree_kind", 
    "clone_nb_seq", 
    "germ_tree_duplicate_seq", 
    "germ_tree_leaf_color", 
    "germ_tree_leaf_shape", 
    "germ_tree_leaf_size", 
    "germ_tree_label_size", 
    "germ_tree_label_hjust", 
    "germ_tree_label_rigth", 
    "germ_tree_label_outside", 
    "germ_tree_right_margin", 
    "germ_tree_legend", 
    "germ_tree_data_assembly_path",
    "germ_tree_meta_path", 
    "germ_tree_meta_legend", 
    "log"
)
tempo.log <- sapply(lapply(tempo.arg, FUN = get, env = sys.nframe(), inherit = FALSE), FUN = is.null)
if(any(tempo.log) == TRUE){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN germ_tree_vizu.R:\n", ifelse(sum(tempo.log, na.rm = TRUE) > 1, "THESE ARGUMENTS\n", "THIS ARGUMENT\n"), paste0(tempo.arg[tempo.log], collapse = "\n"),"\nCANNOT BE NULL")
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end management of NULL arguments, WARNING: only for germ_tree_vizu.R because NULL is "NULL" in the nextflow.config file
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

arg.check2 <- NULL #
text.check2 <- NULL #
checked.arg.names2 <- NULL # for function debbuging: used by r_debugging_tools
ee <- expression(arg.check2 <- c(arg.check2, tempo$problem) , text.check2 <- c(text.check2, tempo$text) , checked.arg.names2 <- c(checked.arg.names2, tempo$object.name))
tempo <- fun_check(data = germ_tree_kind, options = c("rectangular", "roundrect", "slanted", "ellipse", "circular", "fan", "equal_angle", "daylight"), length = 1) ; eval(ee)

if(length(clone_nb_seq) != 1 & any(grepl(clone_nb_seq, pattern = "\\D"))){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN germ_tree_vizu.R:\nTHE clone_nb_seq PARAMETER MUST BE A SINGLE INTEGER\nHERE IT IS: \n", paste0(clone_nb_seq, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}else{
    clone_nb_seq <- as.integer(clone_nb_seq)
    tempo <- fun_check(data = clone_nb_seq, class = "vector", typeof = "integer", length = 1) ; eval(ee)
}

if( ! (length(germ_tree_duplicate_seq) == 1 & any(germ_tree_duplicate_seq %in% c("TRUE", "FALSE")))){ # positive numeric
    tempo.cat <- paste0("ERROR IN germ_tree_vizu.R:\nTHE germ_tree_label_size PARAMETER MUST BE \"TRUE\" OR \"FALSE\"\nHERE IT IS: \n", paste0(germ_tree_duplicate_seq, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}

if(germ_tree_leaf_color == "NULL"){
    germ_tree_leaf_color <- NULL
}else{
    tempo1 <- fun_check(data = germ_tree_leaf_color, class = "vector", mode = "character", length = 1)
    tempo2 <- fun_check(data = germ_tree_leaf_color, class = "vector", typeof = "integer", double.as.integer.allowed = TRUE, length = 1)
    checked.arg.names <- c(checked.arg.names, tempo2$object.name)
    if(tempo1$problem == TRUE & tempo2$problem == TRUE){
        tempo.cat <- paste0("ERROR IN germ_tree_vizu.R\ngerm_tree_leaf_color ARGUMENT MUST BE (1) A HEXADECIMAL COLOR STRING STARTING BY #, OR (2) A COLOR NAME GIVEN BY colors(), OR (3) AN INTEGER VALUE")
        text.check2 <- c(text.check2, tempo.cat)
        arg.check2 <- c(arg.check2, TRUE)
    }else if(tempo1$problem == FALSE & tempo2$problem == TRUE){
        if( ! all(germ_tree_leaf_color %in% colors() | grepl(pattern = "^#", germ_tree_leaf_color), na.rm = TRUE)){
            tempo.cat <- paste0("ERROR IN germ_tree_vizu.R\ngerm_tree_leaf_color ARGUMENT MUST BE (1) A HEXADECIMAL COLOR STRING STARTING BY #, OR (2) A COLOR NAME GIVEN BY colors(), OR (3) AN INTEGER VALUE")
            text.check2 <- c(text.check2, tempo.cat)
            arg.check2 <- c(arg.check2, TRUE)
        }
    }else{
        # no fun_check test here, it is just for checked.arg.names
        tempo <- fun_check(data = germ_tree_leaf_color, class = "vector")
        checked.arg.names <- c(checked.arg.names, tempo$object.name)
    }
}


if(length(germ_tree_leaf_shape) != 1 & any(grepl(germ_tree_leaf_shape, pattern = "\\D"))){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN germ_tree_vizu.R:\nTHE germ_tree_leaf_shape PARAMETER MUST BE A SINGLE INTEGER BETWEEN 0 AND 25\nHERE IT IS: \n", paste0(germ_tree_leaf_shape, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}else{
    germ_tree_leaf_shape <- as.integer(germ_tree_leaf_shape)
    tempo <- fun_check(data = germ_tree_leaf_shape, class = "vector", typeof = "integer", length = 1) ; eval(ee)
    if( ! germ_tree_leaf_shape %in% 0:25){
        tempo.cat <- paste0("ERROR IN germ_tree_vizu.R:\nTHE germ_tree_leaf_shape PARAMETER MUST BE A SINGLE INTEGER BETWEEN 0 AND 25\nHERE IT IS:\n", paste0(germ_tree_leaf_shape, collapse = " "))
        text.check2 <- c(text.check2, tempo.cat)
        arg.check2 <- c(arg.check2, TRUE)
    }
}

if(length(germ_tree_leaf_size) != 1 & any(grepl(germ_tree_leaf_size, pattern = "^[0-9]+\\.{0,1}[0-9]*$"))){ # positive numeric
    tempo.cat <- paste0("ERROR IN germ_tree_vizu.R:\nTHE germ_tree_leaf_size PARAMETER MUST BE A SINGLE POSITIVE NUMBER\nHERE IT IS: \n", paste0(germ_tree_leaf_size, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}else{
    germ_tree_leaf_size <- as.numeric(germ_tree_leaf_size)
    tempo <- fun_check(data = germ_tree_leaf_size, class = "vector", mode = "numeric", neg.values = FALSE, length = 1) ; eval(ee)
}

if(length(germ_tree_label_size) != 1 & any(grepl(germ_tree_label_size, pattern = "^[0-9]+\\.{0,1}[0-9]*$"))){ # positive numeric
    tempo.cat <- paste0("ERROR IN germ_tree_vizu.R:\nTHE germ_tree_label_size PARAMETER MUST BE A SINGLE POSITIVE NUMBER\nHERE IT IS: \n", paste0(germ_tree_label_size, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}else{
    germ_tree_label_size <- as.numeric(germ_tree_label_size)
    tempo <- fun_check(data = germ_tree_label_size, class = "vector", mode = "numeric", neg.values = FALSE, length = 1) ; eval(ee)
}

if(length(germ_tree_label_hjust) != 1 & any(grepl(germ_tree_label_hjust, pattern = "^\\-{0,1}[0-9]+\\.{0,1}[0-9]*$"))){ # positive numeric
    tempo.cat <- paste0("ERROR IN germ_tree_vizu.R:\nTHE germ_tree_label_hjust PARAMETER MUST BE A SINGLE NUMBER\nHERE IT IS: \n", paste0(germ_tree_label_hjust, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}else{
    germ_tree_label_hjust <- as.numeric(germ_tree_label_hjust)
    tempo <- fun_check(data = germ_tree_label_hjust, class = "vector", mode = "numeric", length = 1) ; eval(ee)
}

if( ! (length(germ_tree_label_rigth) == 1 & any(germ_tree_label_rigth %in% c("TRUE", "FALSE")))){ # positive numeric
    tempo.cat <- paste0("ERROR IN germ_tree_vizu.R:\nTHE germ_tree_label_size PARAMETER MUST BE \"TRUE\" OR \"FALSE\"\nHERE IT IS: \n", paste0(germ_tree_label_rigth, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}

if( ! (length(germ_tree_label_outside) == 1 & any(germ_tree_label_outside %in% c("TRUE", "FALSE")))){ # positive numeric
    tempo.cat <- paste0("ERROR IN germ_tree_vizu.R:\nTHE germ_tree_label_outside PARAMETER MUST BE \"TRUE\" OR \"FALSE\"\nHERE IT IS: \n", paste0(germ_tree_label_outside, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}

if(length(germ_tree_right_margin) != 1 & any(grepl(germ_tree_right_margin, pattern = "^[0-9]+\\.{0,1}[0-9]*$"))){ # positive numeric
    tempo.cat <- paste0("ERROR IN germ_tree_vizu.R:\nTHE germ_tree_right_margin PARAMETER MUST BE A SINGLE POSITIVE NUMBER\nHERE IT IS: \n", paste0(germ_tree_right_margin, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}else{
    germ_tree_right_margin <- as.numeric(germ_tree_right_margin)
    tempo <- fun_check(data = germ_tree_right_margin, class = "vector", mode = "numeric", neg.values = FALSE, length = 1) ; eval(ee)
}

if( ! (length(germ_tree_legend) == 1 & any(germ_tree_legend %in% c("TRUE", "FALSE")))){ # positive numeric
    tempo.cat <- paste0("ERROR IN germ_tree_vizu.R:\nTHE germ_tree_legend PARAMETER MUST BE \"TRUE\" OR \"FALSE\"\nHERE IT IS: \n", paste0(germ_tree_legend, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}

if(germ_tree_meta_path == "NULL"){
    germ_tree_meta_path <- NULL
}else if( ! file.exists(germ_tree_meta_path)){
    tempo.cat <- paste0("ERROR IN germ_tree_vizu.R:\nTHE meta_path PARAMETER MUST BE A VALID PATH OF A FILE IF NOT \"NULL\"\nHERE IT IS: \n", paste0(germ_tree_meta_path, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}

if(germ_tree_meta_legend == "NULL"){
    germ_tree_meta_legend <- NULL
}

if(any(arg.check2) == TRUE){ # normally no NA
    stop(paste0("\n\n================\n\n", paste(text.check2[arg.check2], collapse = "\n"), "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between == #
}


# other checkings (not full checked because already checked in the .nf file)
# reserved word checking
# end reserved word checking
# end second round of checking and data preparation
# package checking
# end package checking


################################ End pre-ignition checking


################################ Main code


################ Ignition


fun_report(data = paste0("\n\n################################################################ germ_tree_vizu PROCESS\n\n"), output = log, path = "./", overwrite = FALSE)
ini.date <- Sys.time()
ini.time <- as.numeric(ini.date) # time of process begin, converted into seconds
fun_report(data = paste0("\n\n################################ RUNNING DATE AND STARTING TIME\n\n"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0(ini.date, "\n\n"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0("\n\n################################ RUNNING\n\n"), output = log, path = "./", overwrite = FALSE)


################ End ignition


################ Graphical parameter initialization


################ End graphical parameter initialization


################ Control

empty.tsv <- TRUE
tempo.list <- list.files(path = ".", pattern = "RData$") # gather all the .RData file names into tempo.list
if(length(tempo.list) == 0){
    tempo.cat <- "\\n\\nNO SEQUENCE DATA COMPUTED (NO .RData FILE DETECTED) -> NO GRAPH PLOTTED"
    cat(tempo.cat)
    fun_report(data = tempo.cat, output = log, path = "./", overwrite = FALSE)
    # no need to use pdf(NULL) with fun_gg_empty_graph()
    final.plot <- fun_gg_empty_graph(text = "NO GRAPH PLOTTED\nNO .RData FILE DETECTED", text.size = 3)
    ggplot2::ggsave(filename = paste0("germ_tree.png"), plot = final.plot, device = "png", path = ".", width = 5, height = 5, units = "in", dpi = 300)
    ggplot2::ggsave(filename = paste0("germ_tree.svg"), plot = final.plot, device = "svg", path = ".", width = 5, height = 5, units = "in", dpi = 300)
    ggplot2::ggsave(filename = paste0("germ_tree.pdf"), plot = final.plot, device = "pdf", path = ".", width = 5, height = 5, units = "in", dpi = 300)
    ggplot2::ggsave(filename = paste0("germ_no_tree.png"), plot = final.plot, device = "png", path = ".", width = 5, height = 5, units = "in", dpi = 300)
    ggplot2::ggsave(filename = paste0("germ_no_tree.svg"), plot = final.plot, device = "svg", path = ".", width = 5, height = 5, units = "in", dpi = 300)
    ggplot2::ggsave(filename = paste0("germ_no_tree.pdf"), plot = final.plot, device = "pdf", path = ".", width = 5, height = 5, units = "in", dpi = 300)
    tempo.df <- matrix(c("sequence_id", "clone_id", "clone_name", "chain", "identical_to"), nrow = 1)
    write.table(tempo.df, file = paste0("./germ_tree_dup_seq_not_displayed.tsv"), row.names = FALSE, col.name = FALSE, sep = "\t", quote = FALSE)
}else{
    clone.id <- sapply(tempo.list, FUN = function(x){strsplit(x, split = "_")[[1]][1]})

################ end Control

################ concatenation of all tibble trees into a same tibble and plot


    # gathering trees and db from each .RData in tempo.list and trees -> tempo and db -> tempo.db[[i3]]
    tempo.names <- c("clone_id", "data", "locus", "seqs", "parameters", "trees")
    tempo.names <- setNames(tempo.names, nm = tempo.names) # dplyr::bind_rows use named vectors
    tempo.trees <- dplyr::bind_rows(tempo.names)[0, ]
    tempo.trees <- dplyr::"%>%"(tempo.trees, dplyr::mutate(dplyr::across(clone_id, as.double)))
    tempo.trees <- dplyr::"%>%"(tempo.trees, dplyr::mutate(dplyr::across(data, as.list)))
    tempo.trees <- dplyr::"%>%"(tempo.trees, dplyr::mutate(dplyr::across(locus, as.character)))
    tempo.trees <- dplyr::"%>%"(tempo.trees, dplyr::mutate(dplyr::across(seqs, as.integer)))
    tempo.trees <- dplyr::"%>%"(tempo.trees, dplyr::mutate(dplyr::across(parameters, as.list)))
    tempo.trees <- dplyr::"%>%"(tempo.trees, dplyr::mutate(dplyr::across(trees, as.list)))
    tempo.trees.ini <- tempo.trees
    tempo.clones <- tibble()
    tempo.db <- vector("list", length(tempo.list))
    for(i1 in 1:length(tempo.list)){
        suppressWarnings(rm(trees))
        suppressWarnings(rm(clones))
        suppressWarnings(rm(db))
        load(tempo.list[i1])
        if( ! exists("clones", where = ".GlobalEnv", inherit = FALSE)){
            stop("\n\n========\n\nINTERNAL CODE ERROR 4 IN germ_tree_vizu.R:\nclones OBJECT CANNOT BE ABSENT FROM ", tempo.list[i1], "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n========\n\n")
        }
        # Warning : here [[1]] is ok
        if( ! exists("trees", where = ".GlobalEnv", inherit = FALSE)){ # creation of an empty tibble but with the clone id added
            if( ! is.null(fun_get_message("clones$data[[1]]@clone"))){
                stop("\n\n========\n\nINTERNAL CODE ERROR 5 IN germ_tree_vizu.R:\nclones$data[[1]]@clone OBJECT CANNOT BE ABSENT FROM ", tempo.list[i1], "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n========\n\n")
            }
            trees <- dplyr::bind_rows(tempo.trees.ini, data.frame(clone_id = as.double(clones$data[[1]]@clone), data = NA, locus = as.character(NA), seqs = as.integer(NA), parameters = NA, trees = NA))
        }
        tempo.trees <- dplyr::bind_rows(tempo.trees, trees)
        tempo.clones <- dplyr::bind_rows(tempo.clones, clones)
        tempo.db[[i1]] <- db
    }
    suppressWarnings(rm(trees))
    suppressWarnings(rm(clones))
    suppressWarnings(rm(db))
    # end gathering trees and db from each .RData in tempo.list and trees -> tempo and db -> tempo.db[[i3]]
    trees <- tempo.trees # creation of the trees tibble object, made of all the trees
    if(all(sapply(trees$trees, FUN = is.null))){
        plots <- vector(mode = "list", length = length(tempo.db))
    }else{
        plots <- suppressMessages(dowser::plotTrees( # creation of the plots list object, using the trees tibble object
            # option in https://dowser.readthedocs.io/en/latest/topics/plotTrees/
            trees[ ! sapply(trees$trees, FUN = is.null),], 
            layout = germ_tree_kind, 
            title = TRUE
        ))
    }
    clones <- tempo.clones
    db.list <- tempo.db # creation of the db.list tibble object, made of all the db
    save(list = c("trees", "plots", "clones", "db.list", "clone.id"), file = "./all_trees.RData")


################ end concatenation of all tibble trees into a same tibble and plot


################ Data import


    if( ! is.null(germ_tree_meta_path)){
        new_meta_df <- read.table(germ_tree_data_assembly_path, sep = "\t", header = TRUE, comment.char = "") # new because it is not the meta file anymore but the generated table used here has meta file
    }


################ End Data import

################ Plotting tree

    count.final.plot <- 0
    count.final.no.plot <- 0
    for(i3 in 1:length(db.list)){
        final.plot <- NULL
        final.no.plot <- NULL
        empty.tsv <- TRUE
        tempo.title <- paste0(
            "Clonal Group: ", ifelse(is.null(fun_get_message("clones$data[[i3]]@clone")), clones$data[[i3]]@clone, ""), "\n",
            "Chain: ", ifelse(is.null(fun_get_message("clones$locus[i3]")), clones$locus[i3], ""), "\n", 
            "Clonal Group full name: ", ifelse(is.null(fun_get_message("clones$data[[i3]]@v_gene")) & is.null(fun_get_message("clones$data[[i3]]@j_gene")), paste0(clones$data[[i3]]@v_gene, "_", clones$data[[i3]]@j_gene), ""), "\n",
            "Clone ID: ", ifelse(is.null(fun_get_message("clones$clone_id[i3]")), clones$clone_id[i3], ""), "\n",
            "CDR3 junction length: ", ifelse(is.null(fun_get_message("clones$data[[i3]]@junc_len")), clones$data[[i3]]@junc_len, ""), "\n",
            "Number of leafs (different sequences): ", ifelse(is.null(fun_get_message("clones$data[[i3]]@data")), nrow(clones$data[[i3]]@data), ""), "\n",
            "Number of sequences in the clonal group: ", ifelse(is.null(fun_get_message("nrow(db.list[[i3]])")), nrow(db.list[[i3]]), ""), "\n",
            "Minimal number of different sequences for tree (set by the clone_nb_seq parameter): ", clone_nb_seq
        )
        add.text <- NULL
        pdf(NULL)
        legend.final <- grid::grid.rect(gp = grid::gpar(col = "white"))
        dev.off()
        if( ! (is.null(trees$trees[[i3]]) | length(trees$trees[[i3]]$parameters$nseq) == 0)){
            if(trees$trees[[i3]]$parameters$nseq != 0 | is.na(trees$trees[[i3]]$parameters$nseq)){
                if(any(trees$data[[i3]]@clone != clones$data[[i3]]@clone)){
                    stop("\n\n========\n\nINTERNAL CODE ERROR 6 IN germ_tree_vizu.R FOR CLONE ID ", clone.id[[i3]], ":\ntrees$data[[i3]]@clone AND clones$data[[i3]]@clone SHOULD BE iDENTICAL\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n========\n\n")
                }
                if(any(trees$data[[i3]]@v_gene != clones$data[[i3]]@v_gene)){
                    stop("\n\n========\n\nINTERNAL CODE ERROR 7 IN germ_tree_vizu.R FOR CLONE ID ", clone.id[[i3]], ":\ntrees$data[[i3]]@v_gene AND clones$data[[i3]]@v_gene SHOULD BE iDENTICAL\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n========\n\n")
                }
                if(any(trees$data[[i3]]@j_gene != clones$data[[i3]]@j_gene)){
                    stop("\n\n========\n\nINTERNAL CODE ERROR 8 IN germ_tree_vizu.R FOR CLONE ID ", clone.id[[i3]], ":\ntrees$data[[i3]]@j_gene AND clones$data[[i3]]@j_gene SHOULD BE iDENTICAL\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n========\n\n")
                }
                if(any(trees$data[[i3]]@junc_len != clones$data[[i3]]@junc_len)){
                    stop("\n\n========\n\nINTERNAL CODE ERROR 9 IN germ_tree_vizu.R FOR CLONE ID ", clone.id[[i3]], ":\ntrees$data[[i3]]@junc_len AND clones$data[[i3]]@junc_len SHOULD BE iDENTICAL\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n========\n\n")
                }
                if(nrow(trees$data[[i3]]@data) != nrow(clones$data[[i3]]@data)){
                    stop("\n\n========\n\nINTERNAL CODE ERROR 10 IN germ_tree_vizu.R FOR CLONE ID ", clone.id[[i3]], ":\nnrow(trees$data[[i3]]@data) AND nrow(clones$data[[i3]]@data) SHOULD BE iDENTICAL\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n========\n\n")
                }
                if(any(clone.id[[i3]] != trees$data[[i3]]@clone)){
                    stop("\n\n========\n\nINTERNAL CODE ERROR 11 IN germ_tree_vizu.R FOR CLONE ID ", clone.id[[i3]], ":\nclone.id[[i3]] AND trees$data[[i3]]@clone SHOULD BE iDENTICAL\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n========\n\n")
                }
                # color of leafs and germline
                tempo.graph.info <- ggplot2::ggplot_build(ggtree::ggtree(trees$trees[[i3]], layout = if(germ_tree_kind %in% c("roundrect", "ellipse")){"rectangular"}else{germ_tree_kind}) +  ggtree::geom_tippoint() + ggtree::geom_tiplab()) # plots[[i3]] is equivalent to ggtree::ggtree(trees$trees[[i3]])
                leaf.nodes <- NULL
                tempo.log <- NULL
                for(i4 in 1:length(tempo.graph.info$data)){
                    if(any(names(tempo.graph.info$data[[i4]]) == "label")){
                        leaf.nodes <- c(leaf.nodes, tempo.graph.info$data[[i4]]$node)
                        tempo.log <- c(tempo.log, tempo.graph.info$data[[i4]]$label == "Germline")
                    }
                }
                if(is.null(leaf.nodes)){
                    stop("\n\n========\n\nINTERNAL CODE ERROR 12 IN germ_tree_vizu.R FOR CLONE ID ", paste(unique(db.list[[i3]]$clone_id)), ":\nEMPTY leaf.nodes OBJECT GENERATED\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n========\n\n")
                }
                leaf.node.germline <- leaf.nodes[tempo.log]
                leaf.node.not.germline <- leaf.nodes[ ! tempo.log]
                nodes <- tempo.graph.info$data[[1]]$node # here [[1]] is ok
                tempo.log.germline <- nodes == leaf.node.germline
                tempo.log.leaf <- nodes %in% leaf.nodes
                tip.kind <- rep("Seq", length(nodes)) # kind of tip: either germline or seq
                tip.kind[tempo.log.germline] <- "Germline"
                tip.kind[ ! (tempo.log.leaf | tempo.log.germline)] <- NA
                # end color of leafs and germline
                # title prep and leaf labeling modification
                tempo.v <- trees$data[[i3]]@v_gene
                tempo.j <- trees$data[[i3]]@j_gene
                chain <- substr(tempo.j, 1, 3) # extract the IGH or IGK name
                if(chain != substr(tempo.v, 1, 3)){
                    stop(paste0("\n\n============\n\nERROR IN germ_tree_vizu.R\nTHE CHAIN OF THE clone_id ", trees$data[[i3]]@clone, " IS NOT THE SAME BETWEEN V (", tempo.v, ") AND J (", tempo.j, ")\n\n============\n\n"), call. = FALSE)
                }
                tempo.v <- substring(tempo.v, 4)
                tempo.j <- substring(tempo.j, 4)
                clone.name <- paste0(tempo.v, "_", tempo.j)
                removed.seq <- NULL
                trees$trees[[i3]]$old.tip.label <- trees$trees[[i3]]$tip.label # backup the initial tip.label into a new 13th compartment of the trees$trees[[i3]] list
                trees$trees[[i3]]$new.tip.label <- trees$trees[[i3]]$tip.label # make a new 14th compartment of the trees$trees[[i3]] list for the tip labels with seq removed added
                if(germ_tree_duplicate_seq == "TRUE" & nrow(trees$data[[i3]]@data) != nrow(db.list[[i3]])){
                    stop(paste0("\n\n============\n\nINTERNAL CODE ERROR 13 IN germ_tree_vizu.R for clone ID ", clone.id[[i3]], "\nTHE germ_tree_duplicate_seq PARAMETER IS SET TO \"TRUE\"\nBUT THE NUMBER OF ROWS IN trees$data[[i3]]@data (n=", nrow(trees$data[[i3]]@data), ")\nIS DIFFERENT FROM THE NUMBER OF ROWS IN db (n=", nrow(db.list[[i3]]), ")\nAS IF SOME SEQUENCES WHERE REMOVED\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n============\n\n"), call. = FALSE)
                }else if(germ_tree_duplicate_seq == "TRUE"){
                    tempo.cat <- "All sequences of the tree are displayed"
                    add.text <- paste0(ifelse(is.null(add.text), tempo.cat, paste0(add.text, "\n", tempo.cat)))
                }else if(germ_tree_duplicate_seq == "FALSE" & nrow(trees$data[[i3]]@data) == nrow(db.list[[i3]])){
                    tempo.cat <- "All sequences of the tree are displayed"
                    add.text <- paste0(ifelse(is.null(add.text), tempo.cat, paste0(add.text, "\n", tempo.cat)))
                }else if(germ_tree_duplicate_seq == "FALSE" & nrow(trees$data[[i3]]@data) != nrow(db.list[[i3]])){
                    # get removed sequences info
                    no.removed.seq.log <- db.list[[i3]][[1]] %in% trees$data[[i3]]@data[[1]] # of note, the trees$data[[i3]]@data$collapsed report also the names in the first column trees$data[[i3]]@data[[1]] = trees$data[[i3]]@data$sequence_id. It is not the sequences that are not in the tree anymore !
                    removed.seq.log <- ! no.removed.seq.log
                    if( ! any(removed.seq.log)){
                        stop(paste0("\n\n============\n\nINTERNAL CODE ERROR 14 IN germ_tree_vizu.R for clone ID ", clone.id[[i3]], "\nTHE germ_tree_duplicate_seq PARAMETER IS SET TO \"FALSE\"\nBUT NO SEQ NAMES COLLAPSED IN THE TREE IN trees$data[[", i3, "]]@data[[1]] COMPARED TO db.list[[", i3, "]][[1]]\ntrees$data[[i3]]@data[[1]]: ", paste(trees$data[[i3]]@data[[1]], collapse = " "), "\ndb.list[[", i3, "]][[1]]: ", paste(db.list[[i3]][[1]], collapse = " "), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n============\n\n"), call. = FALSE)
                    }else{
                        not.removed.seq <- db.list[[i3]][[1]][no.removed.seq.log]
                        removed.seq <- db.list[[i3]][[1]][removed.seq.log]
                        tempo.warn <- paste0(length(removed.seq), " sequences removed from the display (Parameter germ_tree_duplicate_seq == \"FALSE\". See germ_tree_dup_seq_not_displayed.tsv)")
                        tempo.warn2 <- paste0("FOR CLONE ID ", paste(unique(db.list[[i3]]$clone_id)), "\n", tempo.warn)
                        warn <- paste0(ifelse(is.null(warn), tempo.warn2, paste0(warn, "\n\n", tempo.warn2)))
                        tempo.cat <- paste0("Warning: ", tempo.warn, "\nNumber of total sequences are indicated between brackets in leaf labelings (removed = total - 1)")
                        add.text <- paste0(ifelse(is.null(add.text), tempo.cat, paste0(add.text, "\n", tempo.cat)))
                        identical.seq <- NULL
                        for(i4 in 1:length(removed.seq)){
                            tempo.log <- grepl(trees$data[[i3]]@data$collapsed, pattern = removed.seq[i4]) # collapsed names are separated by comma during dowser::formatClones()
                            identical.seq <- c(identical.seq, trees$data[[i3]]@data[[1]][tempo.log])
                        }
                        if(any( ! removed.seq %in% unlist(strsplit(trees$data[[i3]]@data$collapsed, split = ",")))){
                            stop(paste0("\n\n============\n\nINTERNAL CODE ERROR 15 IN germ_tree_vizu.R for clone ID ", clone.id[[i3]], "\nTHESE SEQUENCES NAMES, PRESENT IN db.list[[", i3, "]][[1]]: ", paste(db.list[[i3]][[1]], collapse = " "), "\nDOES NOT APPEAR ANYMORE IN removed.seq ", paste(removed.seq, collapse = " "), "\nTHE PROBLEM COMES FROM THE COLLAPSE IN clones$data[[", i3, "]]@data$collapsed: ", paste(clones$data[[", i3, "]]@data$collapsed, collapse = " "), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n============\n\n"), call. = FALSE)
                        }
                        if(length(removed.seq) != length(identical.seq)){
                            stop(paste0("\n\n============\n\nINTERNAL CODE ERROR 16 IN germ_tree_vizu.R for clone ID ", clone.id[[i3]], "\nidentical.seq SHOULD HAVE ", length(removed.seq), " SEQUENCES NAMES\nLENGTH OF removed.seq: ", paste(removed.seq, collapse = " "), "\nSHOULD BE IDENTICAL TO LENGTH OF identical.seq: ", paste(identical.seq, collapse = " "), "\ntrees$data[[i3]]@data$collapsed: ", paste(trees$data[[i3]]@data$collapsed, collapse = " "), "\nclones$data[[", i3, "]]@data$collapsed: ", paste(clones$data[[i3]]@data$collapsed, collapse = " "), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n============\n\n"), call. = FALSE)
                        }
                        suppressWarnings(rm(tempo.df))
                        tempo.df <- data.frame(sequence_id = removed.seq, clone_id = clone.id[[i3]], clone_name = clone.name, chain = chain, identical_to = identical.seq)
                        write.table(tempo.df, file = paste0("./", clone.id[[i3]], "_germ_tree_dup_seq_not_displayed.tsv"), row.names = FALSE, col.name = TRUE, sep = "\t", quote = FALSE)
                        empty.tsv <- FALSE
                        # end get removed sequences info
                        # modif of the tree tip labeling
                        if(length(trees$trees[[i3]]$new.tip.label) - 1 != length(trees$data[[i3]]@data$sequence_id)){ # - 1 because "Germline" label removed
                            stop(paste0("\n\n============\n\nINTERNAL CODE ERROR 17 IN germ_tree_vizu.R for clone ID ", clone.id[[i3]], "\nLENGTH OF trees$trees[[i3]]$new.tip.label SHOULD BE -1 OF LENGTH OF trees$data[[i3]]@data$sequence_id. HERE IT IS:\n\nLENGTH - 1 OF trees$trees[[i3]]$new.tip.label: ", length(trees$trees[[i3]]$new.tip.label) - 1, " \nLENGTH OF trees$data[[i3]]@data$sequence_id: ", length(trees$data[[i3]]@data$sequence_id), "\n\ntrees$trees[[i3]]$new.tip.label:\n", paste(trees$trees[[i3]]$new.tip.label, collapse = "\n"), "\n\ntrees$data[[i3]]@data$sequence_id:\n", paste(trees$data[[i3]]@data$sequence_id, collapse = "\n"), "\n\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n============\n\n"), call. = FALSE)
                        }else{
                            tempo.pos <- match(trees$trees[[i3]]$new.tip.label[-length(trees$trees[[i3]]$new.tip.label)], trees$data[[i3]]@data$sequence_id) # [-length(trees$trees[[i3]]$new.tip.label)] to remove the "Germline" label from the elements
                            if( ! all(trees$data[[i3]]@data$sequence_id[tempo.pos] == trees$trees[[i3]]$new.tip.label[-length(trees$trees[[i3]]$new.tip.label)])){
                                stop(paste0("\n\n============\n\nINTERNAL CODE ERROR 18 IN germ_tree_vizu.R for clone ID ", clone.id[[i3]], "\nELEMENTS SHOULD BE THE SAME. HERE IT IS:\n\ntrees$data[[i3]]@data$sequence_id[tempo.pos]:\n", paste(trees$data[[i3]]@data$sequence_id[tempo.pos], collapse = "\n"), "\n\ntrees$trees[[i3]]$new.tip.label[-length(trees$trees[[i3]]$new.tip.label)]:\n", paste(trees$trees[[i3]]$new.tip.label[-length(trees$trees[[i3]]$new.tip.label)], collapse = "\n"), "\n\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n============\n\n"), call. = FALSE)
                            }else{
                                trees$trees[[i3]]$new.tip.label[-length(trees$trees[[i3]]$new.tip.label)] <- paste0("(", trees$data[[i3]]@data$collapse_count[tempo.pos], ") ", trees$data[[i3]]@data$sequence_id[tempo.pos]) # -1 because it is the number of "removed" sequences, not "total" number of sequences
                            }
                        }
                        # end modif of the tree tip labeling
                    }
                }
                # end title prep and leaf labeling modification
                # ggplot building
                tempo.gg.name <- "gg.indiv.plot."
                tempo.gg.count <- 0
                # tree with metadata
                if( ( ! is.null(germ_tree_meta_path)) & ! is.null(germ_tree_meta_legend)){
                    # merge of the meta data into the ggtree object. See https://yulab-smu.top/treedata-book/chapter7.html#attach-operator
                    if( ! germ_tree_meta_legend %in% names(new_meta_df)){
                        stop(paste0("\n\n============\n\nERROR IN germ_tree_vizu.R for clone ID ", paste(unique(db.list[[i3]]$clone_id)), "\nIF NOT \"NULL\", THE meta_legend PARAMETER MUST BE A COLUMN NAME OF THE meta_path PARAMETER. HERE IT IS:\ngerm_tree_meta_legend: ", germ_tree_meta_legend, "\nCOLUMN NAMES OF meta_path: ", paste(names(new_meta_df), collapse = " "), "\n\n============\n\n"), call. = FALSE)
                    }
                    if( ! any(trees$trees[[i3]]$tip.label %in% new_meta_df[ , 1])){
                        tempo.cat <- paste0("Warning: meta_legend parameter indicated in the nextflow.config file but no metadata present in this tree (check the first column of the clone_assigned_seq.tsv file used here as meta file?).")
                        paste0(ifelse(is.null(add.text), tempo.cat, paste0(add.text, "\n", tempo.cat)))
                    }
                    tempo.added.trees <- ggtree::"%<+%"( # it seems that this command uses the tip.label compartment to merge new_meta_df into ggtree::ggtree(trees$trees[[i3]], layout = germ_tree_kind)
                        ggtree::ggtree(trees$trees[[i3]], layout = germ_tree_kind),
                        new_meta_df
                    )
                    # adding new tip labels into tempo.added.trees (because they are lost with ggtree::"%<+%")
                    tempo.log <- tempo.added.trees$data$label %in% trees$trees[[i3]]$tip.label
                    if(all(tempo.added.trees$data$label[tempo.log] == trees$trees[[i3]]$tip.label)){
                        tempo.added.trees$data$label[tempo.log] <- trees$trees[[i3]]$new.tip.label
                    }else{
                        stop(paste0("\n\n============\n\nINTERNAL CODE ERROR 19 IN germ_tree_vizu.R for clone ID ", clone.id[[i3]], "\nELEMENTS SHOULD BE THE SAME. HERE IT IS:\n\ntempo.added.trees$data$label[tempo.log]:\n", paste(tempo.added.trees$data$label[tempo.log], collapse = "\n"), "\n\ntrees$trees[[i3]]$tip.label:\n", paste(trees$trees[[i3]]$tip.label, collapse = "\n"), "\n\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n============\n\n"), call. = FALSE)
                    }
                    # end adding new tip labels into tempo.added.trees (because they are lost with ggtree::"%<+%")
                    # 
                    tempo.added.trees$data <- tibble::add_column(tempo.added.trees$data, tip.kind)
                    Annotated <- character(length = nrow(tempo.added.trees$data))
                    Annotated[] <- NA
                    Annotated[tempo.added.trees$data$tip.kind == "Seq"] <- "No"
                    if(any( ! is.na(tempo.added.trees$data$Name))){
                        Annotated[ ! is.na(tempo.added.trees$data$Name)] <- "Yes"
                    }
                    tempo.added.trees$data <- tibble::add_column(tempo.added.trees$data, Annotated)
                    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), tempo.added.trees)
                    tempo.col.values <- germ_tree_meta_legend # name of the legend scale for the tippoints
                    assign(tempo.col.values, tempo.added.trees$data[[germ_tree_meta_legend]])
                    if(is.numeric(get(tempo.col.values))){
                        assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::geom_tippoint(
                            ggplot2::aes_string(
                                fill = "tip.kind", 
                                size = if( ! is.null(get(tempo.col.values))){get(tempo.col.values)}else{NA},
                            ),
                            pch = germ_tree_leaf_shape
                        ))
                        tempo.scale <- seq(
                            from = min(new_meta_df[ , germ_tree_meta_legend], na.rm = TRUE), 
                            to = max(new_meta_df[ , germ_tree_meta_legend], na.rm = TRUE), 
                            length.out = 5
                        )
                        assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::scale_size(
                            name = germ_tree_meta_legend,
                            breaks = tempo.scale,
                            labels = formatC(tempo.scale),
                            limits = range(c(new_meta_df[ , germ_tree_meta_legend]), na.rm = TRUE),
                            range = c(0, 5), 
                            guide = ggplot2::guide_legend(
                                override.aes = list(
                                    fill = ifelse(is.null(germ_tree_leaf_color), "tomato", germ_tree_leaf_color)
                                )
                            )
                        ))
                        assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::scale_discrete_manual(
                            aesthetics = "fill", 
                            name = NULL, 
                            values = c("black", ifelse(is.null(germ_tree_leaf_color), "tomato", germ_tree_leaf_color)),
                            labels = c("Germline", "Seq")
                        ))
                        assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::guides(fill = "none")) # no legend for fill in this context
                    }else{
                        if( ! is.null(get(tempo.col.values))){
                            if( ! all(is.na(get(tempo.col.values)))){
                                assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::geom_tippoint(
                                    ggplot2::aes_string(fill = tempo.col.values),
                                    pch = germ_tree_leaf_shape, 
                                    size = germ_tree_leaf_size
                                ))
                            }else{
                                assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::geom_tippoint(
                                    fill = "grey50",
                                    pch = germ_tree_leaf_shape, 
                                    size = germ_tree_leaf_size
                                ))
                            }
                        }else{
                            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::geom_tippoint(
                                fill = "grey50",
                                pch = germ_tree_leaf_shape, 
                                size = germ_tree_leaf_size
                            ))
                        }
                    }
                    if(any(germ_tree_kind %in% c("rectangular", "roundrect", "slanted", "ellipse"))){
                        if(any(names(tempo.added.trees$data) == "Annotated")){
                            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::geom_tiplab(
                                ggplot2::aes(color = Annotated), 
                                hjust = germ_tree_label_hjust,
                                size = germ_tree_label_size,
                                as_ylab = ifelse(germ_tree_label_rigth == "TRUE" & germ_tree_kind == "rectangular", TRUE, FALSE)
                            ))
                            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::scale_color_manual(
                                name = "Annotated", 
                                values = c("black", ifelse(is.null(germ_tree_leaf_color), "tomato", germ_tree_leaf_color)),
                                labels = c("No", "Yes")
                            ))
                            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::guides(
                                color = "none"
                            ))
                        }else{
                            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::geom_tiplab(
                                hjust = germ_tree_label_hjust,
                                size = germ_tree_label_size,
                                as_ylab = ifelse(germ_tree_label_rigth == "TRUE" & germ_tree_kind == "rectangular", TRUE, FALSE)
                            ))
                        }
                    }else if(any(germ_tree_kind %in% c("circular", "fan", "equal_angle", "daylight"))){
                        if(any(names(tempo.added.trees$data) == "Annotated")){
                            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::geom_tiplab(
                                ggplot2::aes(angle = angle, color = Annotated), 
                                hjust = germ_tree_label_hjust,
                                size = germ_tree_label_size
                            ))
                            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::scale_color_manual(
                                name = "Annotated", 
                                values = c("black", ifelse(is.null(germ_tree_leaf_color), "tomato", germ_tree_leaf_color)),
                                labels = c("No", "Yes")
                            ))
                            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::guides(
                                color = "none"
                            ))
                        }else{
                            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::geom_tiplab(
                                ggplot2::aes(angle = angle), 
                                hjust = germ_tree_label_hjust,
                                size = germ_tree_label_size
                            ))
                        }
                    }else{
                        if(any(names(tempo.added.trees$data) == "Annotated")){
                            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::geom_tiplab(
                                ggplot2::aes(color = Annotated), 
                                size = germ_tree_label_size
                            ))
                            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::scale_color_manual(
                                name = "Annotated", 
                                values = c("black", ifelse(is.null(germ_tree_leaf_color), "tomato", germ_tree_leaf_color)),
                                labels = c("No", "Yes")
                            ))
                            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::guides(
                                color = "none"
                            ))
                        }else{
                            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::geom_tiplab(
                                size = germ_tree_label_size
                            ))
                        }
                    }
                # end tree with metadata
                # tree with no metadata
                }else{ # tree with no metadata
                    if(is.null(germ_tree_meta_path) & ! is.null(germ_tree_meta_legend)){
                        tempo.warn <- paste0("FOR CLONE ID ", paste(unique(db.list[[i3]]$clone_id)), "\nTHE meta_legend PARAMETER IS NOT \"NULL\" BUT THE meta_path PARAMETER IS \"NULL\"")
                        warn <- paste0(ifelse(is.null(warn), tempo.warn, paste0(warn, "\n\n", tempo.warn)))
                    }
                    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::ggtree(trees$trees[[i3]], layout = germ_tree_kind))
                    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::geom_tippoint(
                        ggplot2::aes(fill = tip.kind),
                        pch = germ_tree_leaf_shape, 
                        size = germ_tree_leaf_size
                    ))
                    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::scale_discrete_manual(
                        aesthetics = "fill", 
                        name = NULL, 
                        values = c("black", ifelse(is.null(germ_tree_leaf_color), "tomato", germ_tree_leaf_color)),
                        labels = c("Germline", "Seq")
                    ))
                    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::guides(fill = "none")) # no legend for fill in this context
                    if(any(germ_tree_kind %in% c("rectangular", "roundrect", "slanted", "ellipse"))){
                        assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::geom_tiplab(
                            hjust = germ_tree_label_hjust,
                            size = germ_tree_label_size,
                            as_ylab = ifelse(germ_tree_label_rigth == "TRUE" & germ_tree_kind == "rectangular", TRUE, FALSE)
                        ))
                    }else if(any(germ_tree_kind %in% c("circular", "fan", "equal_angle", "daylight"))){
                        assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::geom_tiplab(
                            ggplot2::aes(angle = angle), 
                            hjust = germ_tree_label_hjust,
                            size = germ_tree_label_size
                        ))
                    }else{
                        assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::geom_tiplab(
                            size = germ_tree_label_size
                        ))
                    }
                }
                if( ! any(germ_tree_kind %in% c("circular", "fan"))){
                    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::coord_cartesian( # to extend the text display outside of the plot region
                        clip = ifelse(germ_tree_label_outside == "TRUE", "off", "on")
                    ))
                }
                if( ! any(germ_tree_kind %in% c("circular", "fan", "equal_angle", "daylight"))){
                    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::geom_treescale(width = 0.01, offset = 0.05, fontsize = 1.5, linesize = 0.25))
                }
                assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::theme(plot.margin = ggplot2::margin(t = 0.25, l = 0.1, b = 0.1, r = germ_tree_right_margin, unit = "in")))
                # legend
                if( ( ! is.null(germ_tree_meta_path)) & ! is.null(germ_tree_meta_legend)){
                    bef.final.plot <- ggplot2::ggplot_build(eval(parse(text = paste(paste0(tempo.gg.name, 1:tempo.gg.count), collapse = " + "))))
                    legend.final <- fun_gg_get_legend(ggplot_built = bef.final.plot) # get legend
                    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::guides(fill = "none", color = "none", size = "none")) # inactivate the initial legend
                    if(germ_tree_legend == "FALSE"){ # even if any(unlist(legend.disp)) is TRUE
                        legend.final <- ggplot2::ggplot()+ggplot2::theme_void() # empty graph instead of legend
                    }
                }
                # end legend
                count.final.plot <- count.final.plot + 1
                final.plot <- suppressMessages(suppressWarnings(eval(parse(text = paste(paste0(tempo.gg.name, 1:tempo.gg.count), collapse = " + ")))))
            }else{
                # no need to use pdf(NULL) with fun_gg_empty_graph()
                count.final.no.plot <- count.final.no.plot + 1
                final.no.plot <- fun_gg_empty_graph(text = paste0("NO GRAPH PLOTTED FOR CLONE ID ", paste(unique(db.list[[i3]]$clone_id)), "\nNOT ENOUGH SEQUENCES DETECTED"), text.size = 3)
            }
        }else{
            # no need to use pdf(NULL) with fun_gg_empty_graph()
            count.final.no.plot <- count.final.no.plot + 1
            final.no.plot <- fun_gg_empty_graph(text = paste0("NO GRAPH PLOTTED FOR CLONE ID ", paste(unique(db.list[[i3]]$clone_id)), "\nNOT ENOUGH SEQUENCES DETECTED"), text.size = 3)
        }
        if( ( ! is.null(germ_tree_meta_path)) & ! is.null(germ_tree_meta_legend)){
            tempo.log <- clones$data[[i3]]@data$sequence_id %in% new_meta_df$Name
            if(any(tempo.log)){
                tempo.cat <- paste0("Annotated sequences in this clonal group: ", paste(clones$data[[i3]]@data$sequence_id[tempo.log], collapse = ", "))
                add.text <- paste0(ifelse(is.null(add.text), tempo.cat, paste0(add.text, "\n", tempo.cat)))
            }
        }
        if(is.null(germ_tree_meta_path)){
            legend.width = 0
        }else{
            legend.width = 0.15 # single proportion (between 0 and 1) indicating the relative width of the legend sector (on the right of the plot) relative to the width of the plot. Value 1 means that the window device width is split in 2, half for the plot and half for the legend. Value 0 means no room for the legend, which will overlay the plot region. Write NULL to inactivate the legend sector. In such case, ggplot2 will manage the room required for the legend display, meaning that the width of the plotting region can vary between graphs, depending on the text in the legend
        }
        # end legend
        #title
        tempo.title <- paste0(
            tempo.title, 
            ifelse(is.null(add.text), "", paste0("\n", add.text))
        )
        title.grob <- grid::textGrob(
            label = tempo.title,
            x = grid::unit(0, "lines"), 
            y = grid::unit(0, "lines"),
            hjust = 0,
            vjust = 0,
            gp = grid::gpar(fontsize = 5)
        )
        if( ! is.null(final.plot)){
            pdf(NULL)
            final.plot <- suppressMessages(suppressWarnings(gridExtra::arrangeGrob(grobs = list(final.plot, legend.final), ncol=2, widths=c(1, legend.width), top = title.grob, left = " ", right = " "))) # , left = " ", right = " " : trick to add margins in the plot. padding =  unit(0.5, "inch") is for top margin above the title
            dev.off()
            # saving plots
            plots[[i3]] <- final.plot # just for all_objects.RData
            ggplot2::ggsave(filename = paste0("germ_tree_cloneID_", trees$clone_id[i3], ".png"), plot = final.plot, device = "png", path = ".", width = 5, height = 5, units = "in", dpi = 300)
            ggplot2::ggsave(filename = paste0("germ_tree_cloneID_", trees$clone_id[i3], ".svg"), plot = final.plot, device = "svg", path = ".", width = 5, height = 5, units = "in", dpi = 300)
            ggplot2::ggsave(filename = paste0("germ_tree_cloneID_", trees$clone_id[i3], ".pdf"), plot = final.plot, device = "pdf", path = ".", width = 5, height = 5, units = "in", dpi = 300)
        }else if( ! is.null(final.no.plot)){
            pdf(NULL)
            final.no.plot <- suppressMessages(suppressWarnings(gridExtra::arrangeGrob(grobs = list(final.no.plot, legend.final), ncol=2, widths=c(1, legend.width), top = title.grob, left = " ", right = " "))) # , left = " ", right = " " : trick to add margins in the plot. padding =  unit(0.5, "inch") is for top margin above the title
            dev.off()
            # saving plots
            ggplot2::ggsave(filename = paste0("germ_no_tree_cloneID_", trees$clone_id[i3], ".png"), plot = final.no.plot, device = "png", path = ".", width = 5, height = 5, units = "in", dpi = 300)
            ggplot2::ggsave(filename = paste0("germ_no_tree_cloneID_", trees$clone_id[i3], ".svg"), plot = final.no.plot, device = "svg", path = ".", width = 5, height = 5, units = "in", dpi = 300)
            ggplot2::ggsave(filename = paste0("germ_no_tree_cloneID_", trees$clone_id[i3], ".pdf"), plot = final.no.plot, device = "pdf", path = ".", width = 5, height = 5, units = "in", dpi = 300)
        }else{
            stop("\n\n========\n\nINTERNAL CODE ERROR 20 IN germ_tree_vizu.R FOR CLONE ID ", clone.id[[i3]], ":\n\nfinal.plot AND final.no.plot CANNOT BE BOTH NULL.\n\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n========\n\n")
        }

################ end Plotting tree

################ save empty tsv

        if(empty.tsv == TRUE){
            tempo.df <- matrix(c("sequence_id", "clone_id", "clone_name", "chain", "identical_to"), nrow = 1)
            write.table(tempo.df, file = paste0("./", clone.id[[i3]], "_germ_tree_dup_seq_not_displayed.tsv"), row.names = FALSE, col.name = FALSE, sep = "\t", quote = FALSE)
        }

################ end save empty tsv

################ saving plots

    }
    if(length(list.files(path = ".", pattern = "^germ_tree_cloneID_.*pdf$")) > 0){
        # dowser::treesToPDF(plots, file = "trees.pdf", nrow=2, ncol=2)
        tempo <- qpdf::pdf_combine(input = list.files(path = ".", pattern = "^germ_tree_cloneID_.*pdf$"), output = "./germ_tree.pdf") # assignation to prevent a returned element
    }else{
        final.plot <- fun_gg_empty_graph(text = "NO GRAPH PLOTTED\nNO GRONAL GROUP WITH AT LEAST 3 SEQUENCES", text.size = 3)
        ggplot2::ggsave(filename = paste0("germ_tree.png"), plot = final.plot, device = "png", path = ".", width = 5, height = 5, units = "in", dpi = 300)
        ggplot2::ggsave(filename = paste0("germ_tree.svg"), plot = final.plot, device = "svg", path = ".", width = 5, height = 5, units = "in", dpi = 300)
        ggplot2::ggsave(filename = paste0("germ_tree.pdf"), plot = final.plot, device = "pdf", path = ".", width = 5, height = 5, units = "in", dpi = 300)
    }
    if(length(list.files(path = ".", pattern = "^germ_no_tree_cloneID_.*pdf$")) > 0){
        # dowser::treesToPDF(plots, file = "trees.pdf", nrow=2, ncol=2)
        tempo <- qpdf::pdf_combine(input = list.files(path = ".", pattern = "^germ_no_tree_cloneID_.*pdf$"), output = "./germ_no_tree.pdf") # assignation to prevent a returned element
    }else{
        final.no.plot <- fun_gg_empty_graph(text = "NO GRAPH PLOTTED\nONLY GRONAL GROUPS WITH AT LEAST 3 SEQUENCES", text.size = 3)
        ggplot2::ggsave(filename = paste0("germ_no_tree.png"), plot = final.no.plot, device = "png", path = ".", width = 5, height = 5, units = "in", dpi = 300)
        ggplot2::ggsave(filename = paste0("germ_no_tree.svg"), plot = final.no.plot, device = "svg", path = ".", width = 5, height = 5, units = "in", dpi = 300)
        ggplot2::ggsave(filename = paste0("germ_no_tree.pdf"), plot = final.no.plot, device = "pdf", path = ".", width = 5, height = 5, units = "in", dpi = 300)
    }
}

################ end saving plots


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
    tempo.cat <- paste0("IN germ_tree_vizu.R OF THE NEXFLOW EXECUTION:\n\n", warn)
    fun_report(data = tempo.cat, output = log, path = "./", overwrite = FALSE)
    cat(tempo.cat)
}else{
    fun_report(data = paste0("\n\nNO WARNING MESSAGE TO REPORT"), output = log, path = "./", overwrite = FALSE)
}
on.exit(exp = options(warning.length = ini.warning.length), add = TRUE)


################ end Warning messages


################ Parameter printing


fun_report(data = paste0("\n\n################################ INITIAL SETTINGS OF PARAMETERS"), output = log, path = "./", overwrite = FALSE)
fun_report(data = param.ini.settings, output = log, path = "./", overwrite = FALSE, vector.cat = TRUE)
fun_report(data = paste0("\n\n################################ R SYSTEM AND PACKAGES"), output = log, path = "./", overwrite = FALSE)
tempo <- sessionInfo()
tempo$otherPkgs <- tempo$otherPkgs[order(names(tempo$otherPkgs))] # sort the packages
tempo$loadedOnly <- tempo$loadedOnly[order(names(tempo$loadedOnly))] # sort the packages
fun_report(data = tempo, output = log, path = "./", overwrite = FALSE, vector.cat = TRUE)
fun_report(data = paste0("\n\n################################ JOB END\n\nTIME: ", end.date, "\n\nTOTAL TIME LAPSE: ", total.lapse, "\n"), output = log, path = "./", overwrite = FALSE)


################ end Parameter printing


################################ End Main code
