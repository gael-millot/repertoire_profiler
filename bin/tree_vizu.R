#!/usr/bin/Rscript
#########################################################################
##                                                                     ##
##     tree_vizu.R                                                     ##
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
    stop(paste0("\n\n================\n\nERROR IN tree_vizu.R\n", version$version.string, " IS NOT THE 4.1.2 RECOMMANDED\n\n================\n\n"))
}
# other initializations
erase.objects = TRUE # write TRUE to erase all the existing objects in R before starting the algorithm and FALSE otherwise. Beginners should use TRUE
if(erase.objects == TRUE){
    rm(list = ls(all.names = TRUE))
    erase.objects = TRUE
}
erase.graphs = TRUE # write TRUE to erase all the graphic windows in R before starting the algorithm and FALSE otherwise
script <- "tree_vizu"


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
        stop(paste0("\n\n================\n\nERROR IN tree_vizu.R\nTHE args OBJECT HAS NA\n\n================\n\n"), call. = FALSE)
    }
    tempo.arg.names <- c(
        "tree_kind", 
        "tree_duplicate_seq", 
        "tree_leaf_color", 
        "tree_leaf_shape", 
        "tree_leaf_size", 
        "tree_label_size", 
        "tree_label_hjust", 
        "tree_label_rigth", 
        "tree_label_outside", 
        "tree_right_margin", 
        "tree_legend", 
        "meta_file", 
        "tree_meta_path_names", 
        "cute", 
        "log"
    ) # objects names exactly in the same order as in the bash code and recovered in args. Here only one, because only the path of the config file to indicate after the tree_vizu.R script execution
    if(length(args) != length(tempo.arg.names)){
        stop(paste0("\n\n================\n\nERROR IN tree_vizu.R\nTHE NUMBER OF ELEMENTS IN args (", length(args),") IS DIFFERENT FROM THE NUMBER OF ELEMENTS IN tempo.arg.names (", length(tempo.arg.names),")\nargs:", paste0(args, collapse = ","), "\ntempo.arg.names:", paste0(tempo.arg.names, collapse = ","), "\n\n================\n\n"), call. = FALSE)
    }
    for(i1 in 1:length(tempo.arg.names)){
        assign(tempo.arg.names[i1], args[i1])
    }
    rm(tempo.arg.names, args, i1)
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

# setwd("C:/Users/gael/Documents/Git_projects/ig_clustering/dev/test")
# tree_kind = "rectangular"
# tree_duplicate_seq = "TRUE" 
# tree_leaf_color = "NULL" 
# tree_leaf_shape = "23" 
# tree_leaf_size = "3" 
# tree_label_size = "2" 
# tree_label_hjust = "-0.25" 
# tree_label_rigth = "FALSE" 
# tree_label_outside = "TRUE"
# tree_right_margin = "1.5" 
# tree_legend = "TRUE" 
# meta_file = "./metadata.tsv"
# tree_meta_path_names = "KD"
# cute = "https://gitlab.pasteur.fr/gmillot/cute_little_R_functions/-/raw/v11.8.0/cute_little_R_functions.R"
# log = "tree_vizu.log"
# file.remove(c("./all_objects.RData", "./all_trees.RData", "./trees.pdf", "./tree_vizu.log"))


################################ end Test

################################ Recording of the initial parameters


param.list <- c(
    "erase.objects", 
    "erase.graphs", 
    "script", 
    "run.way",
    if(run.way == "SCRIPT"){"command"}, 
    "tree_kind", 
    "tree_duplicate_seq", 
    "tree_leaf_color", 
    "tree_leaf_shape", 
    "tree_leaf_size", 
    "tree_label_size", 
    "tree_label_hjust", 
    "tree_label_rigth", 
    "tree_label_outside", 
    "tree_right_margin", 
    "tree_legend", 
    "meta_file", 
    "tree_meta_path_names", 
    "cute", 
    "log"
)
if(any(duplicated(param.list))){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 1 IN tree_vizu.R\nTHE param.list OBJECT CONTAINS DUPLICATED ELEMENTS:\n", paste(param.list[duplicated(param.list)], collapse = " "), "\n\n================\n\n"), call. = FALSE) # message for developers
}
if(erase.objects == TRUE){
    created.object.control <- ls()[ ! ls() %in% "param.list"]
    if( ! (all(created.object.control %in% param.list) & all(param.list %in% created.object.control))){
        stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 2 IN tree_vizu.R\nINCONSISTENCIES BETWEEN THE ARGUMENTS USED AND THE PARAMETERS REQUIRED IN THE EXECUTABLE CODE FILE\nTHE ARGUMENTS NOT PRESENT IN THE EXECUTABLE FILE (tree_vizu.R) ARE:\n", paste(created.object.control[ ! created.object.control %in% param.list], collapse = " "), "\nTHE PARAMETERS OF THE EXECUTABLE FILE (tree_vizu.R) NOT PRESENT IN THE ARGUMENTS ARE:\n", paste(param.list[ ! param.list %in% created.object.control], collapse = " "), "\n\n================\n\n"), call. = FALSE) # message for developers
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
    stop(paste0("\n\n============\n\nERROR IN tree_vizu.R\ncute PARAMETER MUST BE LENGTH 1: ", paste(cute, collapse = " "), "\n\n============\n\n"), call. = FALSE)
}else if(grepl(x = cute, pattern = "^http")){
    tempo.try <- try(suppressWarnings(suppressMessages(source(cute, local = .GlobalEnv))), silent = TRUE)
    if(any(grepl(x = tempo.try, pattern = "^[Ee]rror"))){
        stop(paste0("\n\n============\n\nERROR IN tree_vizu.R\nHTTP INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else if( ! grepl(x = cute, pattern = "^http")){
    if( ! file.exists(cute)){
        stop(paste0("\n\n============\n\nERROR IN tree_vizu.R\nFILE INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else{
    tempo.cat <- paste0("\n\n================\n\nINTERNAL CODE ERROR 3 IN tree_vizu.R: CODE HAS TO BE MODIFIED\n\n============\n\n")
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
    tempo.cat <- paste0("ERROR IN tree_vizu.R\nREQUIRED cute FUNCTION", ifelse(length(tempo) > 1, "S ARE", " IS"), " MISSING IN THE R ENVIRONMENT:\n", paste0(tempo, collapse = "()\n"))
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end required function checking


################ local function: package import


# R Packages required
req.package.list <- c(
    "lubridate", 
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
tempo <- fun_check(data = tree_kind, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = tree_duplicate_seq, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = tree_leaf_color, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = tree_leaf_shape, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = tree_leaf_size, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = tree_label_size, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = tree_label_hjust, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = tree_label_rigth, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = tree_label_outside, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = tree_right_margin, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = tree_legend, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = meta_file, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = tree_meta_path_names, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = log, class = "vector", typeof = "character", length = 1) ; eval(ee)
if(any(arg.check) == TRUE){ # normally no NA
    stop(paste0("\n\n================\n\n", paste(text.check[arg.check], collapse = "\n"), "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between == #
}
# end argument primary checking
# second round of checking and data preparation
# management of NA arguments
# end management of NA arguments
# management of NULL arguments
tempo.arg <-c(
    "tree_kind", 
    "tree_duplicate_seq", 
    "tree_leaf_color", 
    "tree_leaf_shape", 
    "tree_leaf_size", 
    "tree_label_size", 
    "tree_label_hjust", 
    "tree_label_rigth", 
    "tree_label_outside", 
    "tree_right_margin", 
    "tree_legend", 
    "meta_file", 
    "tree_meta_path_names", 
    "log"
)
tempo.log <- sapply(lapply(tempo.arg, FUN = get, env = sys.nframe(), inherit = FALSE), FUN = is.null)
if(any(tempo.log) == TRUE){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN tree_vizu.R:\n", ifelse(sum(tempo.log, na.rm = TRUE) > 1, "THESE ARGUMENTS\n", "THIS ARGUMENT\n"), paste0(tempo.arg[tempo.log], collapse = "\n"),"\nCANNOT BE NULL")
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end management of NULL arguments
# seed
set.seed(1)
# end seed
# warning initiation
ini.warning.length <- options()$warning.length
options(warning.length = 8170)
warn <- NULL
# warn.count <- 0 # not required
# end warning initiation
# other checkings

arg.check2 <- NULL #
text.check2 <- NULL #
checked.arg.names2 <- NULL # for function debbuging: used by r_debugging_tools
ee <- expression(arg.check2 <- c(arg.check2, tempo$problem) , text.check2 <- c(text.check2, tempo$text) , checked.arg.names2 <- c(checked.arg.names2, tempo$object.name))
tempo <- fun_check(data = tree_kind, options = c("rectangular", "roundrect", "slanted", "ellipse", "circular", "fan", "equal_angle", "daylight"), length = 1) ; eval(ee)

if( ! (length(tree_duplicate_seq) == 1 & any(tree_duplicate_seq %in% c("TRUE", "FALSE")))){ # positive numeric
    tempo.cat <- paste0("ERROR IN tree_vizu.R:\nTHE tree_label_size PARAMETER MUST BE \"TRUE\" OR \"FALSE\"\nHERE IT IS: \n", paste0(tree_duplicate_seq, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}

if(tree_leaf_color == "NULL"){
    tree_leaf_color <- NULL
}else{
    tempo1 <- fun_check(data = tree_leaf_color, class = "vector", mode = "character", length = 1)
    tempo2 <- fun_check(data = tree_leaf_color, class = "vector", typeof = "integer", double.as.integer.allowed = TRUE, length = 1)
    checked.arg.names <- c(checked.arg.names, tempo2$object.name)
    if(tempo1$problem == TRUE & tempo2$problem == TRUE){
        tempo.cat <- paste0("ERROR IN tree_vizu.R\ntree_leaf_color ARGUMENT MUST BE (1) A HEXADECIMAL COLOR STRING STARTING BY #, OR (2) A COLOR NAME GIVEN BY colors(), OR (3) AN INTEGER VALUE")
        text.check2 <- c(text.check2, tempo.cat)
        arg.check2 <- c(arg.check2, TRUE)
    }else if(tempo1$problem == FALSE & tempo2$problem == TRUE){
        if( ! all(tree_leaf_color %in% colors() | grepl(pattern = "^#", tree_leaf_color), na.rm = TRUE)){
            tempo.cat <- paste0("ERROR IN tree_vizu.R\ntree_leaf_color ARGUMENT MUST BE (1) A HEXADECIMAL COLOR STRING STARTING BY #, OR (2) A COLOR NAME GIVEN BY colors(), OR (3) AN INTEGER VALUE")
            text.check2 <- c(text.check2, tempo.cat)
            arg.check2 <- c(arg.check2, TRUE)
        }
    }else{
        # no fun_check test here, it is just for checked.arg.names
        tempo <- fun_check(data = tree_leaf_color, class = "vector")
        checked.arg.names <- c(checked.arg.names, tempo$object.name)
    }
}


if(length(tree_leaf_shape) != 1 & any(grepl(tree_leaf_shape, pattern = "\\D"))){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN tree_vizu.R:\nTHE tree_leaf_shape PARAMETER MUST BE A SINGLE INTEGER BETWEEN 0 AND 25\nHERE IT IS: \n", paste0(tree_leaf_shape, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}else{
    tree_leaf_shape <- as.integer(tree_leaf_shape)
    tempo <- fun_check(data = tree_leaf_shape, class = "vector", typeof = "integer", length = 1) ; eval(ee)
    if( ! tree_leaf_shape %in% 0:25){
        tempo.cat <- paste0("ERROR IN tree_vizu.R:\nTHE tree_leaf_shape PARAMETER MUST BE A SINGLE INTEGER BETWEEN 0 AND 25\nHERE IT IS:\n", paste0(tree_leaf_shape, collapse = " "))
        text.check2 <- c(text.check2, tempo.cat)
        arg.check2 <- c(arg.check2, TRUE)
    }
}

if(length(tree_leaf_size) != 1 & any(grepl(tree_leaf_size, pattern = "^[0-9]+\\.{0,1}[0-9]*$"))){ # positive numeric
    tempo.cat <- paste0("ERROR IN tree_vizu.R:\nTHE tree_leaf_size PARAMETER MUST BE A SINGLE POSITIVE NUMBER\nHERE IT IS: \n", paste0(tree_leaf_size, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}else{
    tree_leaf_size <- as.numeric(tree_leaf_size)
    tempo <- fun_check(data = tree_leaf_size, class = "vector", mode = "numeric", neg.values = FALSE, length = 1) ; eval(ee)
}

if(length(tree_label_size) != 1 & any(grepl(tree_label_size, pattern = "^[0-9]+\\.{0,1}[0-9]*$"))){ # positive numeric
    tempo.cat <- paste0("ERROR IN tree_vizu.R:\nTHE tree_label_size PARAMETER MUST BE A SINGLE POSITIVE NUMBER\nHERE IT IS: \n", paste0(tree_label_size, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}else{
    tree_label_size <- as.numeric(tree_label_size)
    tempo <- fun_check(data = tree_label_size, class = "vector", mode = "numeric", neg.values = FALSE, length = 1) ; eval(ee)
}

if(length(tree_label_hjust) != 1 & any(grepl(tree_label_hjust, pattern = "^\\-{0,1}[0-9]+\\.{0,1}[0-9]*$"))){ # positive numeric
    tempo.cat <- paste0("ERROR IN tree_vizu.R:\nTHE tree_label_hjust PARAMETER MUST BE A SINGLE NUMBER\nHERE IT IS: \n", paste0(tree_label_hjust, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}else{
    tree_label_hjust <- as.numeric(tree_label_hjust)
    tempo <- fun_check(data = tree_label_hjust, class = "vector", mode = "numeric", length = 1) ; eval(ee)
}

if( ! (length(tree_label_rigth) == 1 & any(tree_label_rigth %in% c("TRUE", "FALSE")))){ # positive numeric
    tempo.cat <- paste0("ERROR IN tree_vizu.R:\nTHE tree_label_size PARAMETER MUST BE \"TRUE\" OR \"FALSE\"\nHERE IT IS: \n", paste0(tree_label_rigth, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}

if( ! (length(tree_label_outside) == 1 & any(tree_label_outside %in% c("TRUE", "FALSE")))){ # positive numeric
    tempo.cat <- paste0("ERROR IN tree_vizu.R:\nTHE tree_label_outside PARAMETER MUST BE \"TRUE\" OR \"FALSE\"\nHERE IT IS: \n", paste0(tree_label_outside, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}

if(length(tree_right_margin) != 1 & any(grepl(tree_right_margin, pattern = "^[0-9]+\\.{0,1}[0-9]*$"))){ # positive numeric
    tempo.cat <- paste0("ERROR IN tree_vizu.R:\nTHE tree_right_margin PARAMETER MUST BE A SINGLE POSITIVE NUMBER\nHERE IT IS: \n", paste0(tree_right_margin, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}else{
    tree_right_margin <- as.numeric(tree_right_margin)
    tempo <- fun_check(data = tree_right_margin, class = "vector", mode = "numeric", neg.values = FALSE, length = 1) ; eval(ee)
}

if( ! (length(tree_legend) == 1 & any(tree_legend %in% c("TRUE", "FALSE")))){ # positive numeric
    tempo.cat <- paste0("ERROR IN tree_vizu.R:\nTHE tree_legend PARAMETER MUST BE \"TRUE\" OR \"FALSE\"\nHERE IT IS: \n", paste0(tree_legend, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}

if(meta_file == "NULL"){
    meta_file <- NULL
}else if( ! file.exists(meta_file)){
    tempo.cat <- paste0("ERROR IN tree_vizu.R:\nTHE meta_file PARAMETER MUST BE A VALID PATH OF A FILE IS NOT \"NULL\"\nHERE IT IS: \n", paste0(meta_file, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}

if(tree_meta_path_names == "NULL"){
    tree_meta_path_names <- NULL
}

if(any(arg.check2) == TRUE){ # normally no NA
    stop(paste0("\n\n================\n\n", paste(text.check2[arg.check2], collapse = "\n"), "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between == #
}


# end other checkings
# reserved word checking
# end reserved word checking
# end second round of checking and data preparation
# package checking
# end package checking


################################ End pre-ignition checking


################################ Main code


################ Ignition


fun_report(data = paste0("\n\n################################################################ tree_vizu PROCESS\n\n"), output = log, path = "./", overwrite = FALSE)
ini.date <- Sys.time()
ini.time <- as.numeric(ini.date) # time of process begin, converted into seconds
fun_report(data = paste0("\n\n################################ RUNNING DATE AND STARTING TIME\n\n"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0(ini.date, "\n\n"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0("\n\n################################ RUNNING\n\n"), output = log, path = "./", overwrite = FALSE)


################ End ignition


################ Graphical parameter initialization


################ End graphical parameter initialization


################ Control


tempo.list <- list.files(path = ".", pattern = "RData")
if(length(tempo.list) == 0){
    cat("\\n\\nNO TREE DATA COMPUTED (no .RData detected) -> NO GRAPH PLOTTED")
}else{


################ end Control

################ concatenation of all tibble trees into a same tibble and plot

    tempo <- NULL
    tempo.db <- vector("list", length(tempo.list))
    for(i3 in 1:length(tempo.list)){
        suppressWarnings(rm(trees))
        suppressWarnings(rm(db))
        load(tempo.list[i3])
        tempo <- dplyr::bind_rows(tempo, trees)
        tempo.db[[i3]] <- db
    }
    suppressWarnings(rm(trees))
    trees <- tempo
    plots <- suppressMessages(dowser::plotTrees(
        # option in https://dowser.readthedocs.io/en/latest/topics/plotTrees/
        trees, 
        layout = tree_kind, 
        title = TRUE
    ))
    db.list <- tempo.db
    save(list = c("trees", "plots", "db.list"), file = "./all_trees.RData")

################ end concatenation of all tibble trees into a same tibble and plot


################ Data import

    if( ! is.null(meta_file)){
        meta.df <- read.table(meta_file, sep = "\t", header = TRUE)
    }


################ End Data import

################ Plotting tree

    for(i3 in 1:length(plots)){
        # color of leafs and germline
        tempo.graph.info <- ggplot2::ggplot_build(ggtree::ggtree(trees$trees[[i3]], layout = if(tree_kind %in% c("roundrect", "ellipse")){"rectangular"}else{tree_kind}) +  ggtree::geom_tippoint() + ggtree::geom_tiplab()) # plots[[i3]] is equivalent to ggtree::ggtree(trees$trees[[i3]])
        leaf.nodes <- NULL
        tempo.log <- NULL
        for(i4 in 1:length(tempo.graph.info$data)){
            if(any(names(tempo.graph.info$data[[i4]]) == "label")){
                leaf.nodes <- c(leaf.nodes, tempo.graph.info$data[[i4]]$node)
                tempo.log <- c(tempo.log, tempo.graph.info$data[[i4]]$label == "Germline")
            }
        }
        if(is.null(leaf.nodes)){
            stop("\n\n========\n\nINTERNAL CODE ERROR 4 IN tree_vizu.R for clone ID ", paste(unique(db.list[[i3]]$clone_id)), ":\nEMPTY leaf.nodes OBJECT GENERATED\n\n========\n\n")
        }
        leaf.node.germline <- leaf.nodes[tempo.log]
        leaf.node.not.germline <- leaf.nodes[ ! tempo.log]
        nodes <- tempo.graph.info$data[[1]]$node
        tempo.log.germline <- nodes == leaf.node.germline
        tempo.log.leaf <- nodes %in% leaf.nodes
        tempo.names <- rep("Seq", length(nodes))
        tempo.names[tempo.log.germline] <- "Germline"
        tempo.names[ ! (tempo.log.leaf | tempo.log.germline)] <- NA
        # end color of leafs and germline
        # ggplot building
        tempo.gg.name <- "gg.indiv.plot."
        tempo.gg.count <- 0
        if(is.null(meta_file)){
            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::ggtree(trees$trees[[i3]], layout = tree_kind))
            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::geom_tippoint(
                ggplot2::aes(fill = tempo.names),
                pch = tree_leaf_shape, 
                size = tree_leaf_size
            ))
            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::scale_discrete_manual(
                aesthetics = "fill", 
                name = NULL, 
                values = c("black", ifelse(is.null(tree_leaf_color), "tomato", tree_leaf_color)),
                labels = c("Germline", "Seq")
            ))
            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::guides(fill = "none")) # never legend for fill in this context
        }else{
            # merge of the meta data into the ggtree object. See https://yulab-smu.top/treedata-book/chapter7.html#attach-operator
            if( ! tree_meta_path_names %in% names(meta.df)){
                stop(paste0("\n\n============\n\nERROR IN tree_vizu.R for clone ID ", paste(unique(db.list[[i3]]$clone_id)), "\nIF NOT \"NULL\", THE tree_meta_path_names PARAMETER MUST BE A COLUMN NAME OF THE tree_meta_path PARAMETER: ", tree_meta_path_names, "\n\n============\n\n"), call. = FALSE)
            }
            tempo.added.trees <- ggtree::"%<+%"(
                ggtree::ggtree(trees$trees[[i3]], layout = tree_kind),
                meta.df
            )
            tempo.added.trees$data <- tibble::add_column(tempo.added.trees$data, tempo.names)
            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), tempo.added.trees)
            tempo.col.values <- tree_meta_path_names
            assign(tempo.col.values, tempo.added.trees$data[[tree_meta_path_names]])
            if(is.numeric(get(tempo.col.values))){
                assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::geom_tippoint(
                    ggplot2::aes_string(
                        fill = "tempo.names", 
                        size = if( ! is.null(get(tempo.col.values))){get(tempo.col.values)}else{NA},
                    ),
                    pch = tree_leaf_shape
                ))
                tempo.scale <- seq(
                    from = min(meta.df[ , tree_meta_path_names], na.rm = TRUE), 
                    to = max(meta.df[ , tree_meta_path_names], na.rm = TRUE), 
                    length.out = 5
                )
                assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::scale_size(
                    name = tree_meta_path_names,
                    breaks = tempo.scale,
                    labels = formatC(tempo.scale),
                    limits = range(c(meta.df[ , tree_meta_path_names]), na.rm = TRUE),
                    range = c(0, 5), 
                    guide = ggplot2::guide_legend(
                        override.aes = list(
                            fill = ifelse(is.null(tree_leaf_color), "tomato", tree_leaf_color)
                        )
                    )
                ))
                assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::scale_discrete_manual(
                    aesthetics = "fill", 
                    name = NULL, 
                    values = c("black", ifelse(is.null(tree_leaf_color), "tomato", tree_leaf_color)),
                    labels = c("Germline", "Seq")
                ))
                assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::guides(fill = "none")) # never legend for fill in this context
            }else{
                assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::geom_tippoint(
                    ggplot2::aes_string(
                        fill = if( ! is.null(get(tempo.col.values))){tempo.col.values}else{NA}
                    ),
                    pch = tree_leaf_shape, 
                    size = tree_leaf_size
                ))
            }
        }
        if(any(tree_kind %in% c("rectangular", "roundrect", "slanted", "ellipse"))){
            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::geom_tiplab(
                hjust = tree_label_hjust,
                size = tree_label_size,
                as_ylab = ifelse(tree_label_rigth == "TRUE" & tree_kind == "rectangular", TRUE, FALSE)
            ))
        }else if(any(tree_kind %in% c("circular", "fan", "equal_angle", "daylight"))){
            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::geom_tiplab(
                ggplot2::aes(angle = angle), 
                hjust = tree_label_hjust,
                size = tree_label_size
            ))
        }else{
            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::geom_tiplab(
                size = tree_label_size
            ))
        }
        if( ! any(tree_kind %in% c("circular", "fan"))){
            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::coord_cartesian( # to extend the text display outside of the plot region
                clip = ifelse(tree_label_outside == "TRUE", "off", "on")
            ))
        }
        if( ! any(tree_kind %in% c("circular", "fan", "equal_angle", "daylight"))){
            assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggtree::geom_treescale(width = 0.01, offset = 0.05))
        }
        assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::theme(plot.margin = ggplot2::margin(t = 0.25, l = 0.1, b = 0.1, r = tree_right_margin, unit = "in")))

        # end ggplot building
        # legend
        bef.final.plot <- ggplot2::ggplot_build(eval(parse(text = paste(paste0(tempo.gg.name, 1:tempo.gg.count), collapse = " + "))))
        legend.final <- fun_gg_get_legend(ggplot_built = bef.final.plot) # get legend
        assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::guides(fill = "none", color = "none", size = "none")) # inactivate the initial legend
        if(tree_legend == "FALSE"){ # even if any(unlist(legend.disp)) is TRUE
            legend.final <- ggplot2::ggplot()+ggplot2::theme_void() # empty graph instead of legend
        }
        if(is.null(meta_file)){
            legend.width = 0
        }else{
            legend.width = 0.15 # single proportion (between 0 and 1) indicating the relative width of the legend sector (on the right of the plot) relative to the width of the plot. Value 1 means that the window device width is split in 2, half for the plot and half for the legend. Value 0 means no room for the legend, which will overlay the plot region. Write NULL to inactivate the legend sector. In such case, ggplot2 will manage the room required for the legend display, meaning that the width of the plotting region can vary between graphs, depending on the text in the legend
        }
        # end legend
        # title
        tempo.v <- trees$data[[i3]]@v_gene
        tempo.j <- trees$data[[i3]]@j_gene
        chain <- substr(tempo.j, 1, 3) # extract the IGH or IGK name
        if(chain != substr(tempo.v, 1, 3)){
            stop(paste0("\n\n============\n\nERROR IN tree_vizu.R\nTHE CHAIN OF THE clone_id ", trees$data[[i3]]@clone, " IS NOT THE SAME BETWEEN V (", tempo.v, ") AND J (", tempo.j, ")\n\n============\n\n"), call. = FALSE)
        }
        tempo.v <- substring(tempo.v, 4)
        tempo.j <- substring(tempo.j, 4)
        clone.name <- paste0(tempo.v, "_", tempo.j)
        clone.id <-  trees$data[[i3]]@clone
        removed.seq <- NULL
        add.text <- NULL
        if(tree_duplicate_seq == "TRUE" & nrow(trees$data[[i3]]@data) != nrow(db.list[[i3]])){
            stop(paste0("\n\n============\n\nINTERNAL CODE ERROR 5 IN tree_vizu.R for clone ID ", clone.id, "\nTHE tree_duplicate_seq PARAMETER IS SET TO \"TRUE\"\nBUT THE NUMBER OF ROWS IN trees$data[[i3]]@data (n=", nrow(trees$data[[i3]]@data), ")\nIS DIFFERENT FROM THE NUMBER OF ROWS IN db (n=", nrow(db.list[[i3]]), ")\nAS IF SOME SEQUENCES WHERE REMOVED\n\n============\n\n"), call. = FALSE)
        }else if(tree_duplicate_seq == "FALSE" & nrow(trees$data[[i3]]@data) == nrow(db.list[[i3]])){
            add.text <- "All sequences of the tree displayed"
        }else if(tree_duplicate_seq == "FALSE" & nrow(trees$data[[i3]]@data) != nrow(db.list[[i3]])){
            # get removed sequences info
            duplic.seq.log <- ! db.list[[i3]][[1]] %in% trees$data[[i3]]@data[[1]]
            if( ! any(duplic.seq.log)){
                stop(paste0("\n\n============\n\nINTERNAL CODE ERROR 6 IN tree_vizu.R for clone ID ", clone.id, "\nTHE tree_duplicate_seq PARAMETER IS SET TO \"FALSE\"\nBUT NO SEQ NAMES REMOVED FROM THE TREE IN trees$data[[i3]]@data[[1]] IS DIFFERENT FROM THE NUMBER OF ROWS IN db (n=", nrow(db.list[[i3]]), ")\ntrees$data[[i3]]@data[[1]]: ", paste(trees$data[[i3]]@data[[1]], collapse = " "), "\ndb.list[[i3]][[1]]: ", paste(db.list[[i3]][[1]], collapse = " "), "\n\n============\n\n"), call. = FALSE)
            }else{
                removed.seq <- db.list[[i3]][[1]][duplic.seq.log]
                add.text <- "Warning: sequences removed from the display (Parameter tree_duplicate_seq == \"FALSE\". See seq_not_displayed.tsv)"
                identical.seq <- vector("character", length(removed.seq))
                tempo.pos <- which(duplic.seq.log)
                for(i4 in 1:length(tempo.pos)){
                    for(i5 in trees$data[[i3]]@data$sequence_id){
                        tempo.log <- db.list[[i3]]$d_identity[tempo.pos[i4]] != db.list[[i3]]$d_identity[db.list[[i3]]$sequence_id == i5] | db.list[[i3]]$j_identity[tempo.pos[i4]] != db.list[[i3]]$j_identity[db.list[[i3]]$sequence_id == i5]
                        if(tempo.log){
                            identical.seq[i4] <- i5
                        }
                    }
                }
                if(any(identical.seq == "")){
                    stop(paste0("\n\n============\n\nINTERNAL CODE ERROR 7 IN tree_vizu.R for clone ID ", clone.id, "\nidentical.seq SHOULD HAVE ", length(tempo.pos), " SEQUENCES NAMES (NO REMAINING EMPTY SLOT): ", paste(identical.seq, collapse = " "), "\n\n============\n\n"), call. = FALSE)
                }
                tempo.df <- data.frame(sequence_id = removed.seq, clone_id = clone.id, clone_name = clone.name, chain = chain, identical_to = identical.seq)
                write.table(tempo.df, file = paste0("./", clone.id, "_seq_not_displayed.tsv"), row.names = FALSE, col.name = TRUE, sep = "\t")
                # end get removed sequences info
            }
        }

        tempo.title <- paste0(
            "Clonal Group: ", clone.name, "\n",
            "Chain: ", chain, "\n", 
            "Clonal Group full name: ", trees$data[[i3]]@v_gene, "_", trees$data[[i3]]@j_gene, "\n",
            "Clone ID: ", trees$data[[i3]]@clone, "\n",
            "CDR3 junction length: ", trees$data[[i3]]@junc_len, "\n",
            "Number of leafs: ", nrow(trees$data[[i3]]@data), "\n",
            "Number of sequences in the clonal group: ", nrow(db.list[[i3]]), 
            ifelse(is.null(add.text), "", paste0("\n", add.text))
        )
        title.grob <- grid::textGrob(
            label = tempo.title,
            x = grid::unit(0, "lines"), 
            y = grid::unit(0, "lines"),
            hjust = 0,
            vjust = 0,
            gp = grid::gpar(fontsize = 6)
        )
        # end title
        pdf(NULL)
        final.plot <- suppressMessages(suppressWarnings(gridExtra::arrangeGrob(grobs = list(eval(parse(text = paste(paste0(tempo.gg.name, 1:tempo.gg.count), collapse = " + "))), legend.final), ncol=2, widths=c(1, legend.width), top = title.grob)))

################ end Plotting tree

################ saving plots


        plots[[i3]] <- final.plot
        ggplot2::ggsave(filename = paste0(trees$clone_id[i3], ".png"), plot = plots[[i3]], device = "png", path = ".", width = 5, height = 5, units = "in", dpi = 300)
        ggplot2::ggsave(filename = paste0(trees$clone_id[i3], ".svg"), plot = plots[[i3]], device = "svg", path = ".", width = 5, height = 5, units = "in", dpi = 300)
        ggplot2::ggsave(filename = paste0(trees$clone_id[i3], ".pdf"), plot = plots[[i3]], device = "pdf", path = ".", width = 5, height = 5, units = "in", dpi = 300)
    }
    # dowser::treesToPDF(plots, file = "trees.pdf", nrow=2, ncol=2)
    qpdf::pdf_combine(input = list.files(path = ".", pattern = ".pdf$"), output = "./trees.pdf")
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
    fun_report(data = paste0("\n\n", warn), output = log, path = "./", overwrite = FALSE)
}else{
    fun_report(data = paste0("\n\nNO WARNING MESSAGE TO REPORT"), output = log, path = "./", overwrite = FALSE)
}


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
