#!/usr/bin/Rscript
#########################################################################
##                                                                     ##
##     get_tree.R                                                      ##
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
if(version$version.string != "R version 4.0.5 (2021-03-31)"){
    stop(paste0("\n\n================\n\nERROR IN get_tree.R\n", version$version.string, " IS NOT THE 4.0.5 RECOMMANDED\n\n================\n\n"))
}
# other initializations
erase.objects = TRUE # write TRUE to erase all the existing objects in R before starting the algorithm and FALSE otherwise. Beginners should use TRUE
if(erase.objects == TRUE){
    rm(list = ls(all.names = TRUE))
    erase.objects = TRUE
}
erase.graphs = TRUE # write TRUE to erase all the graphic windows in R before starting the algorithm and FALSE otherwise
script <- "get_tree"


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
        stop(paste0("\n\n================\n\nERROR IN get_tree.R\nTHE args OBJECT HAS NA\n\n================\n\n"), call. = FALSE)
    }
    tempo.arg.names <- c(
        "seq_name_remplacement_ch", 
        "meta_file", 
        "clone_nb_seq", 
        "tree_duplicate_seq", 
        "igphylm_exe_path", 
        "cute", 
        "log"
    ) # objects names exactly in the same order as in the bash code and recovered in args. Here only one, because only the path of the config file to indicate after the get_tree.R script execution
    if(length(args) != length(tempo.arg.names)){
        stop(paste0("\n\n================\n\nERROR IN get_tree.R\nTHE NUMBER OF ELEMENTS IN args (", length(args),") IS DIFFERENT FROM THE NUMBER OF ELEMENTS IN tempo.arg.names (", length(tempo.arg.names),")\nargs:", paste0(args, collapse = ","), "\ntempo.arg.names:", paste0(tempo.arg.names, collapse = ","), "\n\n================\n\n"), call. = FALSE)
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


# Can only be run on Lunix because igphylm cannot be installed on Windows
# To run locally on WSL2
# copy-paste the output in /home/gael/20211126_dejoux/dataset/
# Then:
# sudo docker run -ti -v /home/gael/20211126_dejoux/:/caca gmillot/immcantation_v1.2:gitlab_v10.1
# Inside the image:
# cd /caca/dataset/fd7f938fe6ebcb2bdbefb2334d45c6/
# R
# Inside R:
# getwd() # to check that we are indeed in /caca/dataset/fd7f938fe6ebcb2bdbefb2334d45c6/
# copy-paste this code:

# seq_name_remplacement_ch = "17_renamed_seq.tsv"
# tree_duplicate_seq = "FALSE" 
# meta_file = "metadata_sort1.tsv"
# clone_nb_seq = "3"
# igphylm_exe_path = "/usr/local/share/igphyml/src/igphyml"
# cute = "cute_little_R_functions.R"
# log = "get_tree.log"




################################ end Test

################################ Recording of the initial parameters


param.list <- c(
    "erase.objects", 
    "erase.graphs", 
    "script", 
    "run.way",
    "tempo.arg.names", 
    if(run.way == "SCRIPT"){"command"}, 
    "seq_name_remplacement_ch", 
    "meta_file", 
    "clone_nb_seq", 
    "tree_duplicate_seq", 
    "igphylm_exe_path", 
    "cute", 
    "log"
)
if(any(duplicated(param.list))){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 1 IN get_tree.R\nTHE param.list OBJECT CONTAINS DUPLICATED ELEMENTS:\n", paste(param.list[duplicated(param.list)], collapse = " "), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE) # message for developers
}
if(erase.objects == TRUE){
    created.object.control <- ls()[ ! ls() %in% "param.list"]
    if( ! (all(created.object.control %in% param.list) & all(param.list %in% created.object.control))){
        stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 2 IN get_tree.R\nINCONSISTENCIES BETWEEN THE ARGUMENTS USED AND THE PARAMETERS REQUIRED IN THE EXECUTABLE CODE FILE\nTHE ARGUMENTS NOT PRESENT IN THE EXECUTABLE FILE (get_tree.R) ARE:\n", paste(created.object.control[ ! created.object.control %in% param.list], collapse = " "), "\nTHE PARAMETERS OF THE EXECUTABLE FILE (get_tree.R) NOT PRESENT IN THE ARGUMENTS ARE:\n", paste(param.list[ ! param.list %in% created.object.control], collapse = " "), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE) # message for developers
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
    stop(paste0("\n\n============\n\nERROR IN get_tree.R\ncute PARAMETER MUST BE LENGTH 1: ", paste(cute, collapse = " "), "\n\n============\n\n"), call. = FALSE)
}else if(grepl(x = cute, pattern = "^http")){
    tempo.try <- try(suppressWarnings(suppressMessages(source(cute, local = .GlobalEnv))), silent = TRUE)
    if(any(grepl(x = tempo.try, pattern = "^[Ee]rror"))){
        stop(paste0("\n\n============\n\nERROR IN get_tree.R\nHTTP INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else if( ! grepl(x = cute, pattern = "^http")){
    if( ! file.exists(cute)){
        stop(paste0("\n\n============\n\nERROR IN get_tree.R\nFILE INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else{
    tempo.cat <- paste0("\n\n================\n\nINTERNAL CODE ERROR 3 IN get_tree.R: CODE HAS TO BE MODIFIED\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n============\n\n")
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
    tempo.cat <- paste0("ERROR IN get_tree.R\nREQUIRED cute FUNCTION", ifelse(length(tempo) > 1, "S ARE", " IS"), " MISSING IN THE R ENVIRONMENT:\n", paste0(tempo, collapse = "()\n"))
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end required function checking


################ local function: package import


# R Packages required
req.package.list <- c(
    "lubridate", 
    "dowser",
    "alakazam",
    "tibble"
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
# management of NULL arguments, WARNING: only for get_tree.R because NULL is "NULL" in the nextflow.config file
tempo.arg <-c(
    "seq_name_remplacement_ch", 
    "meta_file", 
    "clone_nb_seq", 
    "tree_duplicate_seq", 
    "igphylm_exe_path", 
    "log"
)
tempo.log <- sapply(lapply(tempo.arg, FUN = get, env = sys.nframe(), inherit = FALSE), FUN = is.null)
if(any(tempo.log) == TRUE){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN get_tree.R:\n", ifelse(sum(tempo.log, na.rm = TRUE) > 1, "THESE ARGUMENTS\n", "THIS ARGUMENT\n"), paste0(tempo.arg[tempo.log], collapse = "\n"),"\nCANNOT BE NULL")
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end management of NULL arguments, WARNING: only for get_tree.R because NULL is "NULL" in the nextflow.config file
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


if( ! file.exists(seq_name_remplacement_ch)){
    tempo.cat <- paste0("ERROR IN get_tree.R:\nTHE seq_name_remplacement_ch PARAMETER MUST BE A VALID PATH OF A FILE IF NOT \"NULL\"\nHERE IT IS: \n", paste0(seq_name_remplacement_ch, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}

if(meta_file == "NULL"){
    meta_file <- NULL
}else if( ! file.exists(meta_file)){
    tempo.cat <- paste0("ERROR IN get_tree.R:\nTHE meta_file PARAMETER MUST BE A VALID PATH OF A FILE IF NOT \"NULL\"\nHERE IT IS: \n", paste0(meta_file, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}

if(length(clone_nb_seq) != 1 & any(grepl(clone_nb_seq, pattern = "\\D"))){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN get_tree.R:\nTHE clone_nb_seq PARAMETER MUST BE A SINGLE INTEGER\nHERE IT IS: \n", paste0(clone_nb_seq, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}else{
    clone_nb_seq <- as.integer(clone_nb_seq)
    tempo <- fun_check(data = clone_nb_seq, class = "vector", typeof = "integer", length = 1) ; eval(ee)
}

if( ! (length(tree_duplicate_seq) == 1 & any(tree_duplicate_seq %in% c("TRUE", "FALSE")))){ # positive numeric
    tempo.cat <- paste0("ERROR IN get_tree.R:\nTHE tree_label_size PARAMETER MUST BE \"TRUE\" OR \"FALSE\"\nHERE IT IS: \n", paste0(tree_duplicate_seq, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}else if(tree_duplicate_seq == "TRUE"){
    tree_duplicate_seq <- TRUE
}else{
    tree_duplicate_seq <- FALSE
}


if( ! file.exists(igphylm_exe_path)){
    tempo.cat <- paste0("ERROR IN get_tree.R:\nTHE igphylm FILE IS NOT CORRECTLY IMPORTED INTO THE NEXTFLOW work DIRECTORY BECAUSE igphylm_exe_path PARAMETER MUST BE A VALID PATH\nHERE IT IS: \n", paste0(igphylm_exe_path, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
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


fun_report(data = paste0("\n\n################################################################ get_tree PROCESS\n\n"), output = log, path = "./", overwrite = FALSE)
ini.date <- Sys.time()
ini.time <- as.numeric(ini.date) # time of process begin, converted into seconds
fun_report(data = paste0("\n\n################################ RUNNING DATE AND STARTING TIME\n\n"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0(ini.date, "\n\n"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0("\n\n################################ RUNNING\n\n"), output = log, path = "./", overwrite = FALSE)


################ End ignition




################ Data import



clone.id <- strsplit(seq_name_remplacement_ch, split = "_")[[1]][1] # name of the file split into array with the first being the clone ID
db <- alakazam::readChangeoDb(seq_name_remplacement_ch)
db2 <- tibble::add_column(db, collapsed = db[[1]]) # add a column to detected collapsed identical sequences names
clones <- dowser::formatClones(
    data = db2, 
    seq = "sequence_alignment", 
    clone = "clone_id", #  All entries in this column should be identical
    minseq = 1, # return an error if the number of seq in db (i.e., number of rows) is lower than minseq. REcontroled here but already controled by the i [[]] above
    heavy = NULL, # name of heavy chain locus (default = "IGH")
    text_fields = "collapsed", # column of db named "collapsed" (columns to retain during duplicate removal, that will be present in the clones tibble but collapsed depending on the line removal) 
    dup_singles = FALSE, # Duplicate sequences in singleton clones to include them as trees? Always use FALSE. Otherwise, it will create an artificial duplicated sequences in thr tree building so that we will have two instead of one sequence in the tree, after identical seq removal. See https://rdrr.io/cran/dowser/src/R/Clones.R
    trait = if(tree_duplicate_seq){names(db2)[1]}else{NULL} # control the removal of identical sequences but different sequence name
    # max_mask: maximum number of characters to mask at the leading and trailing sequence ends. If NULL then the upper masking bound will be automatically determined from the maximum number of observed leading or trailing Ns amongst all sequences. If set to 0 (default) then masking will not be performed
)
if(nrow(clones$data[[1]]@data) < clone_nb_seq){
    cat(paste0("\nLESS THAN ", clone_nb_seq, " DIFFERENT SEQUENCES SET BY THE clone_nb_seq PARAMETER OF THE repertoire_profiler.config FILE: NO TREE COMPUTED\n"))
    file.copy(from = paste0("./", seq_name_remplacement_ch), to = "./tree_dismissed_seq.tsv")
    write.table(matrix(names(db), nrow = 1), file = "./seq_for_tree.tsv", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t") # empty file
    file.create("./tree_clone_id.tsv") # empty file
    write.table(clone.id, file = "./tree_dismissed_clone_id.tsv", row.names = FALSE, col.names = FALSE, quote = FALSE)
    save(list = c("db", "clones"), file = paste0("./", clones$clone_id, "_get_tree_cloneID.RData"))
}else{
    file.copy(from = paste0("./", seq_name_remplacement_ch), to = "./seq_for_tree.tsv")
    file.create("./tree_dismissed_clone_id.tsv") # empty file
    write.table(matrix(names(db), nrow = 1), file = "./tree_dismissed_seq.tsv", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t") # empty file
    # add the metadata names in the first column of clones$data[[1]]@data to have them in final trees
    if(( ! is.null(meta_file)) & ( ! tree_duplicate_seq)){
        meta <- read.table(meta_file, header = TRUE, sep = "\t")
        tempo.df <- clones$data[[1]]@data # it is "class data.frame"
        for(i4 in 1:nrow(tempo.df)){
            tempo <- strsplit(tempo.df$collapsed[i4], split = ",")[[1]] # recover collapsed names
            tempo.log <- tempo %in% meta$Name
            if(any(tempo %in% meta$Name)){
                tempo.df$sequence_id[i4] <- paste(tempo[tempo.log], collapse = "_") # all the annotated names are in the first column now # info not lost in tempo.df$sequence_id[i4] because all the names are in tempo.df$collapsed
            }
        }
        clones$data[[1]]@data <- tempo.df
    }
    # end add the metadata names in the first column of clones$data[[1]]@data to have them in final trees
    trees <- dowser::getTrees(clones, build = "igphyml", exec = igphylm_exe_path, nproc = 10)
    # assign(paste0("c", trees$clone_id, "_trees"), trees)
    save(list = c("trees", "db", "clones"), file = paste0("./", trees$clone_id, "_get_tree_cloneID.RData"))
    write.table(trees$clone_id, file = "./tree_clone_id.tsv", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

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
    tempo.cat <- paste0("IN get_tree.R OF THE NEXFLOW EXECUTION:\n\n", warn)
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
