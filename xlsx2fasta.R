#!/usr/bin/env Rscript

#########################################################################
##                                                                     ##
##     xlsx2tsv.R                                                      ##
##                                                                     ##
##     Gael A. Millot                                                  ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Computational Biology Department                                ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################




################################ Aim


# Creates all the fasta files from a .xlsx file.


################################ End Aim


################################ Introduction


# The xlsx imported file must contain these column names: "ClusterId", "VH_Full_length_DNA", "VL_Full_length_DNA", "Roc_spe", "Vec_spe" and "Panc_spe"
# And must not have an empty sequence in "VH_Full_length_DNA" or "VL_Full_length_DNA" columns
# To use xlsx2tsv.R, 1) open it, 2) complete the "Parameters that need to be set by the user" section below, 3) save the modifications and 4) run.


################################ End Introduction


################################ Acknowlegments


################################ End Acknowlegments


################################ Initialization


# R version checking
if(version$version.string != "R version 4.1.2 (2021-11-01)"){
    stop(paste0("\n\n================\n\nERROR IN xlsx2tsv.R\n", version$version.string, " IS NOT THE 4.1.2 RECOMMANDED\n\n================\n\n"))
}
# other initializations
erase.objects = TRUE # write TRUE to erase all the existing objects in R before starting the algorithm and FALSE otherwise. Beginners should use TRUE
if(erase.objects == TRUE){
    rm(list = ls(all.names = TRUE))
    erase.objects = TRUE
}
erase.graphs = TRUE # write TRUE to erase all the graphic windows in R before starting the algorithm and FALSE otherwise
script <- "xlsx2tsv"


################################ End Initialization


################################ Parameters that need to be set by the user


path <- "Z://ROCURONIUM PROJECT//01 Primary data//04.Repertoire analysis//SORT2//SORT2 Seq-refiltered//ROC-Sort2(421)-20reads_specificity.xlsx" # true values of x_lim if x_lim = "whole", as imported through nextflow
out.path <- "./"
cute <- "https://gitlab.pasteur.fr/gmillot/cute_little_R_functions/-/raw/v11.4.0/cute_little_R_functions.R" # 
log <- "xlsx2tsv.log"


################################ End Parameters that need to be set by the user


################################ Config import


tempo.cat <- "KIND OF RUN (SCRIPT, COPY-PASTE OR SOURCE): "
if(interactive() == FALSE){ # if(grepl(x = commandArgs(trailingOnly = FALSE), pattern = "R\\.exe$|\\/R$|Rcmd\\.exe$|Rcmd$|Rgui\\.exe$|Rgui$|Rscript\\.exe$|Rscript$|Rterm\\.exe$|Rterm$")){ # detection of script usage
    run.way <- "SCRIPT"
    cat(paste0("\n\n", tempo.cat, run.way, "\n"))
    command <- paste0(commandArgs(trailingOnly = FALSE), collapse = ",") # recover the full command
    args <- commandArgs(trailingOnly = TRUE) # recover arguments written after the call of the R script
    if(any(is.na(args))){
        stop(paste0("\n\n================\n\nERROR IN xlsx2tsv.R\nTHE args OBJECT HAS NA\n\n================\n\n"), call. = FALSE)
    }
    tempo.arg.names <- c(
        "path", 
        "out.path", 
        "cute", 
        "log"
    ) # objects names exactly in the same order as in the bash code and recovered in args. Here only one, because only the path of the config file to indicate after the xlsx2tsv.R script execution
    if(length(args) != length(tempo.arg.names)){
        stop(paste0("\n\n================\n\nERROR IN xlsx2tsv.R\nTHE NUMBER OF ELEMENTS IN args (", length(args),") IS DIFFERENT FROM THE NUMBER OF ELEMENTS IN tempo.arg.names (", length(tempo.arg.names),")\nargs:", paste0(args, collapse = ","), "\ntempo.arg.names:", paste0(tempo.arg.names, collapse = ","), "\n\n================\n\n"), call. = FALSE)
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

# path <- "C:/Users/gael/Documents/Git_projects/fisher_for_vcf/dataset/fisher.tsv"
# cute <- "https://gitlab.pasteur.fr/gmillot/cute_little_R_functions/-/raw/v11.4.0/cute_little_R_functions.R" 
# log <- "xlsx2tsv_report.txt"


################################ end Test

################################ Recording of the initial parameters


param.list <- c(
    "erase.objects", 
    "erase.graphs", 
    "script", 
    "run.way",
    if(run.way == "SCRIPT"){"command"}, 
    "path", 
    "out.path",
    "cute", 
    "log"
)
if(any(duplicated(param.list))){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 1 IN xlsx2tsv.R\nTHE param.list OBJECT CONTAINS DUPLICATED ELEMENTS:\n", paste(param.list[duplicated(param.list)], collapse = " "), "\n\n================\n\n"), call. = FALSE) # message for developers
}
if(erase.objects == TRUE){
    created.object.control <- ls()[ ! ls() %in% "param.list"]
    if( ! (all(created.object.control %in% param.list) & all(param.list %in% created.object.control))){
        stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 2 IN xlsx2tsv.R\nINCONSISTENCIES BETWEEN THE ARGUMENTS USED AND THE PARAMETERS REQUIRED IN THE EXECUTABLE CODE FILE\nTHE ARGUMENTS NOT PRESENT IN THE EXECUTABLE FILE (xlsx2tsv.R) ARE:\n", paste(created.object.control[ ! created.object.control %in% param.list], collapse = " "), "\nTHE PARAMETERS OF THE EXECUTABLE FILE (xlsx2tsv.R) NOT PRESENT IN THE ARGUMENTS ARE:\n", paste(param.list[ ! param.list %in% created.object.control], collapse = " "), "\n\n================\n\n"), call. = FALSE) # message for developers
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
    stop(paste0("\n\n============\n\nERROR IN xlsx2tsv.R\ncute PARAMETER MUST BE LENGTH 1: ", paste(cute, collapse = " "), "\n\n============\n\n"), call. = FALSE)
}else if(grepl(x = cute, pattern = "^http")){
    tempo.try <- try(suppressWarnings(suppressMessages(source(cute, local = .GlobalEnv))), silent = TRUE)
    if(any(grepl(x = tempo.try, pattern = "^[Ee]rror"))){
        stop(paste0("\n\n============\n\nERROR IN xlsx2tsv.R\nHTTP INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else if( ! grepl(x = cute, pattern = "^http")){
    if( ! file.exists(cute)){
        stop(paste0("\n\n============\n\nERROR IN xlsx2tsv.R\nFILE INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else{
    tempo.cat <- paste0("\n\n================\n\nINTERNAL CODE ERROR 3 IN xlsx2tsv.R: CODE HAS TO BE MODIFIED\n\n============\n\n")
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
    tempo.cat <- paste0("ERROR IN xlsx2tsv.R\nREQUIRED cute FUNCTION", ifelse(length(tempo) > 1, "S ARE", " IS"), " MISSING IN THE R ENVIRONMENT:\n", paste0(tempo, collapse = "()\n"))
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end required function checking

################ end import functions from cute little functions toolbox

################ local function: package import

# R Packages required
req.package.list <- c(
    "lubridate", 
    "rJava", # dependency of xlsx. Warning: requires the install of java 64.bit, not 32, using https://www.java.com/fr/download/manual.jsp. Otherwise, error message: le chargement du package ou de l'espace de noms a échoué pour ‘xlsx’ :.onLoad a échoué dans loadNamespace() pour 'rJava', détails
    "xlsx"
)
# for(i in 1:length(req.package.list)){suppressMessages(library(req.package.list[i], character.only = TRUE))}
fun_pack(req.package = req.package.list, load = TRUE, lib.path = NULL) # packages are imported even if inside functions are written as package.name::function() in the present code

################ end local function: package import

################ other functions

################ end other functions

################################ End Functions


################################ Pre-ignition checking


# reserved words
# end reserved words
# argument primary checking
arg.check <- NULL #
text.check <- NULL #
checked.arg.names <- NULL # for function debbuging: used by r_debugging_tools
ee <- expression(arg.check <- c(arg.check, tempo$problem) , text.check <- c(text.check, tempo$text) , checked.arg.names <- c(checked.arg.names, tempo$object.name))
tempo <- fun_check(data = path, class = "vector", typeof = "character", length = 1) ; eval(ee)
# tempo <- fun_check(data = cute, class = "vector", typeof = "character", length = 1) ; eval(ee)
# cute already tested above
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
    "path", 
    "out.path",
    "cute", 
    "log"
)
tempo.log <- sapply(lapply(tempo.arg, FUN = get, env = sys.nframe(), inherit = FALSE), FUN = is.null)
if(any(tempo.log) == TRUE){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN xlsx2tsv.R:\n", ifelse(sum(tempo.log, na.rm = TRUE) > 1, "THESE ARGUMENTS\n", "THIS ARGUMENT\n"), paste0(tempo.arg[tempo.log], collapse = "\n"),"\nCANNOT BE NULL")
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end management of NULL arguments
# management of ""
tempo.arg <-c(
    "path", 
    "out.path",
    "cute", 
    "log"
)
tempo.log <- sapply(lapply(tempo.arg, FUN = get, env = sys.nframe(), inherit = FALSE), FUN = function(x){any(x == "")})
if(any(tempo.log) == TRUE){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN xlsx2tsv.R:\n", ifelse(sum(tempo.log, na.rm = TRUE) > 1, "THESE ARGUMENTS\n", "THIS ARGUMENT\n"), paste0(tempo.arg[tempo.log], collapse = "\n"),"\nCANNOT BE \"\"")
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end management of ""
# code that protects set.seed() in the global environment
# end code that protects set.seed() in the global environment
# warning initiation
ini.warning.length <- options()$warning.length
options(warning.length = 8170)
warn <- NULL
# warn.count <- 0 # not required
# end warning initiation
# other checkings

if( ! file.exists(path)){
    stop(paste0("\n\n============\n\nERROR IN xlsx2tsv.R\nFILE INDICATED IN THE path PARAMETER DOES NOT EXISTS:\n", path, "\n\n============\n\n"), call. = FALSE)
}

if( ! dir.exists(out.path)){
    stop(paste0("\n\n============\n\nERROR IN xlsx2tsv.R\nTHE DIRECTORY INDICATED IN THE out.path PARAMETER DOES NOT EXISTS:\n", out.path, "\n\n============\n\n"), call. = FALSE)
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

ini.date <- Sys.time()
ini.time <- as.numeric(ini.date) # time of process begin, converted into seconds
out.path <- paste0(out.path, "/xlsx_to_fasta_", as.integer(ini.time))
dir.create(out.path)
fun_report(data = paste0("\n\n################################################################ xlsx2tsv PROCESS\n\n"), output = log, path = out.path, overwrite = TRUE)
fun_report(data = paste0("\n\n################################ RUNNING DATE AND STARTING TIME\n\n"), output = log, path = out.path, overwrite = FALSE)
fun_report(data = paste0(ini.date, "\n\n"), output = log, path = out.path, overwrite = FALSE)
fun_report(data = paste0("\n\n################################ RUNNING\n\n"), output = log, path = out.path, overwrite = FALSE)


################ End ignition


################ Graphical parameter initialization


################ End graphical parameter initialization


################ Data import


obs <- xlsx::read.xlsx(path, sheetIndex = 1, header = TRUE)


################ end Data import


############ check

if(length(obs) == 0 || nrow(obs) == 0){
    stop(paste0("\n\n============\n\nERROR IN xlsx2tsv.R\nFILE INDICATED IN THE path PARAMETER IS EMPTY:\n", path, "\n\n============\n\n"), call. = FALSE)
}

if( ! all(c("ClusterId", "VH_Full_length_DNA", "VL_Full_length_DNA") %in% names(obs))){
    stop(paste0("\n\n============\n\nERROR IN xlsx2tsv.R\nIMPORTED FILE:\n", path, "\nMUST CONTAIN THESE COLUMN NAMES:\n\"ClusterId\", \"VH_Full_length_DNA\", \"VL_Full_length_DNA\", \"Roc_spe\", \"Vec_spe\" AND \"Panc_spe\"\n\n============\n\n"), call. = FALSE)
}

for(i0 in c("VH_Full_length_DNA", "VL_Full_length_DNA")){
    tempo.log <- is.na(obs[ , i0]) | obs[ , i0] == ""
    if(any(tempo.log)){
        stop(paste0("\n\n============\n\nERROR IN xlsx2tsv.R\nIMPORTED FILE:\n", path, "\nHAS AN EMPTY SEQUENCE IN THE \"VH_Full_length_DNA\" OR \"VL_Full_length_DNA\" COLUMNS IN LINES:\n", paste(which(tempo.log), collapse = ", "), "\n\n============\n\n"), call. = FALSE)
    }
}


############ end check


############ main


path.all <- paste0(out.path, "/All")
dir.create(path.all)
for(i0 in 1:nrow(obs)){
    tempo.cat <- paste0(">", obs$ClusterId[i0], "_VH\n", obs$VH_Full_length_DNA[i0])
    cat(tempo.cat, file = paste0(path.all, "/", obs$ClusterId[i0], "_VH_all.fasta"))
    tempo.cat <- paste0(">", obs$ClusterId[i0], "_VL\n", obs$VL_Full_length_DNA[i0])
    cat(tempo.cat, file = paste0(path.all, "/", obs$ClusterId[i0], "_VL_all.fasta"))
}

for(i0 in c("Roc_spe", "Vec_spe", "Panc_spe")){
    tempo.path <- paste0(out.path, "/", i0)
    dir.create(tempo.path)
    tempo.log <- obs[ ! is.na(obs[ , i0]), i0] == "T"
    if(all(is.na(obs[ , i0])) | ! (length(tempo.log) != 0 & any(tempo.log))){
        tempo.warn <- paste0(i0, " COLUMN HAS NO \"T\"")
        cat(paste0("\n\nWARNING: ", tempo.warn, "\n\n"))
        fun_report(data = paste0("WARNING\n", tempo.warn), output = log, path = out.path, overwrite = FALSE)
        warn <- paste0(ifelse(is.null(warn), tempo.warn, paste0(warn, "\n\n", tempo.warn)))
    }else{
        for(i2 in which(tempo.log)){
            tempo.cat <- paste0(">", obs$ClusterId[i2], "_VH\n", obs$VH_Full_length_DNA[i2])
            cat(tempo.cat, file = paste0(tempo.path, "/", obs$ClusterId[i2], "_VH_", i0, ".fasta"))
            tempo.cat <- paste0(">", obs$ClusterId[i2], "_VL\n", obs$VL_Full_length_DNA[i2])
            cat(tempo.cat, file = paste0(tempo.path, "/", obs$ClusterId[i2], "_VL_", i0, ".fasta"))
        }
    }
}


############ end main


################ Seeding inactivation


################ end Seeding inactivation


################ Environment saving


fun_report(data = paste0("\n\n################################ RUNNING END"), output = log, path = out.path, overwrite = FALSE)
end.date <- Sys.time()
end.time <- as.numeric(end.date)
total.lapse <- round(lubridate::seconds_to_period(end.time - ini.time))
fun_report(data = paste0("\n\nEND TIME: ", end.date), output = log, path = out.path, overwrite = FALSE)
fun_report(data = paste0("\n\nTOTAL TIME LAPSE: ", total.lapse), output = log, path = out.path, overwrite = FALSE)
fun_report(data = paste0("\n\nALL DATA SAVED IN all_objects.RData"), output = log, path = out.path, overwrite = FALSE)


################ end Environment saving


################ Warning messages


fun_report(data = paste0("\n\n################################ RECAPITULATION OF WARNING MESSAGES"), output = log, path = out.path, overwrite = FALSE)
if( ! is.null(warn)){
    fun_report(data = paste0("\n\n", warn), output = log, path = out.path, overwrite = FALSE)
}else{
    fun_report(data = paste0("\n\nNO WARNING MESSAGE TO REPORT"), output = log, path = out.path, overwrite = FALSE)
}


################ end Warning messages


################ Parameter printing


fun_report(data = paste0("\n\n################################ INITIAL SETTINGS OF PARAMETERS"), output = log, path = out.path, overwrite = FALSE)
fun_report(data = param.ini.settings, output = log, path = out.path, overwrite = FALSE, , vector.cat = TRUE)
fun_report(data = paste0("\n\n################################ R SYSTEM AND PACKAGES"), output = log, path = out.path, overwrite = FALSE)
tempo <- sessionInfo()
tempo$otherPkgs <- tempo$otherPkgs[order(names(tempo$otherPkgs))] # sort the packages
tempo$loadedOnly <- tempo$loadedOnly[order(names(tempo$loadedOnly))] # sort the packages
fun_report(data = tempo, output = log, path = out.path, overwrite = FALSE, , vector.cat = TRUE)
fun_report(data = paste0("\n\n################################ JOB END\n\nTIME: ", end.date, "\n\nTOTAL TIME LAPSE: ", total.lapse, "\n"), output = log, path = out.path, overwrite = FALSE)


################ end Parameter printing


################################ End Main code







