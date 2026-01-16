#!/usr/bin/Rscript
#########################################################################
##                                                                     ##
##     Igblast_chain_check.R                                           ##
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
    stop(paste0("\n\n================\n\nERROR IN repertoire.R\n", version$version.string, " IS NOT THE 4.1.2 RECOMMANDED\n\n================\n\n"))
}
# other initializations
erase.objects = TRUE # write TRUE to erase all the existing objects in R before starting the algorithm and FALSE otherwise. Beginners should use TRUE
if(erase.objects == TRUE){
    rm(list = ls(all.names = TRUE))
    erase.objects = TRUE
}
erase.graphs = TRUE # write TRUE to erase all the graphic windows in R before starting the algorithm and FALSE otherwise
script <- "igblast_chain_check"


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
        stop(paste0("\n\n================\n\nERROR IN repertoire.R\nTHE args OBJECT HAS NA\n\n================\n\n"), call. = FALSE)
    }
    tempo.arg.names <- c(
        "productive_seq", 
        "v_ref_files", 
        "d_ref_files", 
        "j_ref_files", 
        "c_ref_files", 
        "cute", 
        "log"
    ) # objects names exactly in the same order as in the bash code and recovered in args. Here only one, because only the path of the config file to indicate after the repertoire.R script execution
    if(length(args) != length(tempo.arg.names)){
        stop(paste0("\n\n================\n\nERROR IN repertoire.R\nTHE NUMBER OF ELEMENTS IN args (", length(args),") IS DIFFERENT FROM THE NUMBER OF ELEMENTS IN tempo.arg.names (", length(tempo.arg.names),")\nargs:", paste0(args, collapse = ","), "\ntempo.arg.names:", paste0(tempo.arg.names, collapse = ","), "\n\n================\n\n"), call. = FALSE)
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

# setwd("C:/Users/gmillot/Documents/Git_projects/repertoire_profiler/work/95/654c68e7d6a0d7fb5e5e56841741d7")
# productive_seq <- "productive_seq_init.tss"
# v_ref_files <- "imgt_human_IGHV.tss"
# d_ref_files <- "imgt_human_IGHD.tss"
# j_ref_files <- "imgt_human_IGHJ.tss"
# c_ref_files <- "imgt_human_IGHC.tss"
# cute = "cute_little_R_functions_v12.8.R"
# log = "repertoire.log"


################################ end Test

################################ Recording of the initial parameters


param.list <- c(
    "erase.objects", 
    "erase.graphs", 
    "script", 
    "run.way",
    "tempo.arg.names", 
    if(run.way == "SCRIPT"){"command"},   
    "productive_seq", 
    "v_ref_files", 
    "d_ref_files", 
    "j_ref_files", 
    "c_ref_files", 
    "cute", 
    "log"
)
if(any(duplicated(param.list))){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 1 IN repertoire.R\nTHE param.list OBJECT CONTAINS DUPLICATED ELEMENTS:\n", paste(param.list[duplicated(param.list)], collapse = " "), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE) # message for developers
}
if(erase.objects == TRUE){
    created.object.control <- ls()[ ! ls() %in% "param.list"]
    if( ! (all(created.object.control %in% param.list) & all(param.list %in% created.object.control))){
        stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 2 IN repertoire.R\nINCONSISTENCIES BETWEEN THE ARGUMENTS USED AND THE PARAMETERS REQUIRED IN THE EXECUTABLE CODE FILE\nTHE ARGUMENTS NOT PRESENT IN THE EXECUTABLE FILE (repertoire.R) ARE:\n", paste(created.object.control[ ! created.object.control %in% param.list], collapse = " "), "\nTHE PARAMETERS OF THE EXECUTABLE FILE (repertoire.R) NOT PRESENT IN THE ARGUMENTS ARE:\n", paste(param.list[ ! param.list %in% created.object.control], collapse = " "), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE) # message for developers
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
    stop(paste0("\n\n============\n\nERROR IN repertoire.R\ncute PARAMETER MUST BE LENGTH 1: ", paste(cute, collapse = " "), "\n\n============\n\n"), call. = FALSE)
}else if(grepl(x = cute, pattern = "^http")){
    tempo.try <- try(suppressWarnings(suppressMessages(source(cute, local = .GlobalEnv))), silent = TRUE)
    if(any(grepl(x = tempo.try, pattern = "^[Ee]rror"))){
        stop(paste0("\n\n============\n\nERROR IN repertoire.R\nHTTP INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else if( ! grepl(x = cute, pattern = "^http")){
    if( ! file.exists(cute)){
        stop(paste0("\n\n============\n\nERROR IN repertoire.R\nFILE INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else{
    tempo.cat <- paste0("\n\n================\n\nINTERNAL CODE ERROR 3 IN repertoire.R:\nCODE HAS TO BE MODIFIED\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n============\n\n")
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
    tempo.cat <- paste0("ERROR IN repertoire.R\nREQUIRED cute FUNCTION", ifelse(length(tempo) > 1, "S ARE", " IS"), " MISSING IN THE R ENVIRONMENT:\n", paste0(tempo, collapse = "()\n"))
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end required function checking


################ local function: package import


# R Packages required
req.package.list <- c(
    "lubridate"
)
# for(i in 1:length(req.package.list)){suppressMessages(library(req.package.list[i], character.only = TRUE))}
fun_pack(req.package = req.package.list, load = TRUE, lib_path = NULL) # packages are imported even if inside functions are written as package.name::function() in the present code


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
# management of NULL arguments, WARNING: only for repertoire.R because NULL is "NULL" in the nextflow.config file
tempo.arg <-c( 
    "productive_seq", 
    "v_ref_files", 
    "d_ref_files", 
    "j_ref_files", 
    "c_ref_files", 
    "log"
)
tempo.log <- sapply(lapply(tempo.arg, FUN = get, env = sys.nframe(), inherit = FALSE), FUN = is.null)
if(any(tempo.log) == TRUE){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN repertoire.R:\n", ifelse(sum(tempo.log, na.rm = TRUE) > 1, "THESE ARGUMENTS\n", "THIS ARGUMENT\n"), paste0(tempo.arg[tempo.log], collapse = "\n"),"\nCANNOT BE NULL")
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end management of NULL arguments, WARNING: only for repertoire.R because NULL is "NULL" in the nextflow.config file
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

v_file_names <- strsplit(v_ref_files, split = " ")[[1]]
if(length(v_file_names) == 0){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 4 IN repertoire.R:\nPROBLEM WITH v_ref_files: ", paste(v_ref_files, collapse = " "), "PLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
}
d_file_names <- strsplit(d_ref_files, split = " ")[[1]]
if(length(v_file_names) == 0){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 4 IN repertoire.R:\nPROBLEM WITH d_ref_files: ", paste(d_ref_files, collapse = " "), "PLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
}
j_file_names <- strsplit(j_ref_files, split = " ")[[1]]
if(length(j_file_names) == 0){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 4 IN repertoire.R:\nPROBLEM WITH j_ref_files: ", paste(j_ref_files, collapse = " "), "PLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
}
c_file_names <- strsplit(c_ref_files, split = " ")[[1]]
if(length(c_file_names) == 0){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 4 IN repertoire.R:\nPROBLEM WITH c_ref_files: ", paste(c_ref_files, collapse = " "), "PLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
}
# other checkings (not full checked because already checked in the .nf file)
# reserved word checking
# end reserved word checking
# end second round of checking and data preparation
# package checking
# end package checking


################################ End pre-ignition checking


################################ Main code

################ internal variables


################ end internal variables

################ Ignition


fun_report(data = paste0("\n\n################################################################ repertoire PROCESS\n\n"), output = log, path = "./", overwrite = FALSE)
ini.date <- Sys.time()
ini.time <- as.numeric(ini.date) # time of process begin, converted into seconds
fun_report(data = paste0("\n\n################################ RUNNING DATE AND STARTING TIME\n\n"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0(ini.date, "\n\n"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0("\n\n################################ RUNNING\n\n"), output = log, path = "./", overwrite = FALSE)


################ End ignition


################ Graphical parameter initialization


################ End graphical parameter initialization


################ Data import


df <- read.table(productive_seq, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
v_vec <- NULL
for(i0 in 1:length(v_file_names)){
    tempo <- sub(x = v_file_names[i0], pattern = ".fasta$", replacement = ".tsv")
    v_vec <- c(v_vec, scan(file = tempo, what = character()))
}
if(d_file_names != "NULL"){
    d_vec <- NULL
    for(i0 in 1:length(d_file_names)){
        tempo <- sub(x = d_file_names[i0], pattern = ".fasta$", replacement = ".tsv")
        d_vec <- c(d_vec, scan(file = tempo, what = character()))
    }
}
j_vec <- NULL
for(i0 in 1:length(j_file_names)){
    tempo <- sub(x = j_file_names[i0], pattern = ".fasta$", replacement = ".tsv")
    j_vec <- c(j_vec, scan(file = tempo, what = character()))
}
c_vec <- NULL
for(i0 in 1:length(c_file_names)){
    tempo <- sub(x = c_file_names[i0], pattern = ".fasta$", replacement = ".tsv")
    c_vec <- c(c_vec, scan(file = tempo, what = character()))
}
fun_report(data = paste0("\n\nSEQ NAME:\n", paste0(df[ , 1], collapse = "\n"), "\n\n"), output = log, path = "./", overwrite = FALSE)

################ End Data import




################ data modification, plotting and saving


# names of the columns of the .tsv to deal with
if(d_file_names != "NULL"){
    col_names <- c("v_call", "d_call", "j_call", "c_call")
    file_names <- c("v_vec", "d_vec", "j_vec", "c_vec") # warning : same number and same order !
}else{
    col_names <- c("v_call", "j_call", "c_call")
    file_names <- c("v_vec", "j_vec", "c_vec") # warning : same number and same order !
}

if( ! all(col_names %in% names(df))){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 7 IN repertoire.R:\nPROBLEM WITH THE NAMES OF productive_seq THAT MUST CONTAIN \"v_call\", \"d_call\", \"j_call\", \"c_call\".\nHERE IT IS:\n", paste(names(df), collapse = "\n"), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
}

not_good_line_nb <- NULL
for(i0 in 1:length(col_names)){
    tempo_list <- lapply(X = df[ , col_names[i0]], FUN = function(x){
        if(is.na(x)){
            NA
        }else{
            strsplit(x, split = ",")[[1]]
        }
    }) # if several alleles
    tempo_log <- sapply(X = tempo_list, FUN = function(x){
        if(any(is.na(x))){
            FALSE
        }else{
            ! all(x %in% get(file_names[i0]), na.rm = TRUE)
        }
    })
    if(sum(tempo_log, na.rm = TRUE) > 0){
        not_good_line_nb <- c(not_good_line_nb, which(tempo_log))
    }
}


if(is.null(not_good_line_nb)){
    selected_df <- df
    not_selected_df <- df[ -(1:nrow(df)), ]
}else{
    not_good_line_nb <- unique(not_good_line_nb)
    selected_df <- df[ - not_good_line_nb, ]
    not_selected_df <- df[not_good_line_nb, ]
}

if(nrow(selected_df) > 0){ # required. Otherwive empty line data frame sent into the next parallelized process -> error
    write.table(selected_df, file = paste0("./selected.tsv"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}
write.table(not_selected_df, file = paste0("./not_selected.tsv"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


################ end data modification, plotting and saving


################ Pdf window closing


################ end Pdf window closing


################ Seeding inactivation


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
    tempo.cat <- paste0("IN repertoire.R OF THE NEXFLOW EXECUTION:\n\n", warn)
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

