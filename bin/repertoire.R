#!/usr/bin/Rscript
#########################################################################
##                                                                     ##
##     repertoire.R                                                     ##
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
script <- "repertoire"


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
        "igblast_database_path", 
        "file_assembly_ch", 
        "repertoire_names_ch", 
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

# setwd("C:/Users/gael/Documents/Git_projects/ig_clustering/work/c0/15b238a2ee5928549fc971d797bf57")
# igblast_database_path = "germlines/imgt/mouse/vdj"
# file_assembly_ch <- "productive_seq.tsv"
# repertoire_names_ch <- "imgt_mouse_IGHD.tsv imgt_mouse_IGHJ.tsv imgt_mouse_IGHV.tsv"
# cute = "https://gitlab.pasteur.fr/gmillot/cute_little_R_functions/-/raw/v12.4.0/cute_little_R_functions.R"
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
    "igblast_database_path", 
    "file_assembly_ch", 
    "repertoire_names_ch", 
    "cute", 
    "log"
)
if(any(duplicated(param.list))){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 1 IN repertoire.R\nTHE param.list OBJECT CONTAINS DUPLICATED ELEMENTS:\n", paste(param.list[duplicated(param.list)], collapse = " "), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/ig_clustering OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE) # message for developers
}
if(erase.objects == TRUE){
    created.object.control <- ls()[ ! ls() %in% "param.list"]
    if( ! (all(created.object.control %in% param.list) & all(param.list %in% created.object.control))){
        stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 2 IN repertoire.R\nINCONSISTENCIES BETWEEN THE ARGUMENTS USED AND THE PARAMETERS REQUIRED IN THE EXECUTABLE CODE FILE\nTHE ARGUMENTS NOT PRESENT IN THE EXECUTABLE FILE (repertoire.R) ARE:\n", paste(created.object.control[ ! created.object.control %in% param.list], collapse = " "), "\nTHE PARAMETERS OF THE EXECUTABLE FILE (repertoire.R) NOT PRESENT IN THE ARGUMENTS ARE:\n", paste(param.list[ ! param.list %in% created.object.control], collapse = " "), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/ig_clustering OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE) # message for developers
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
    tempo.cat <- paste0("\n\n================\n\nINTERNAL CODE ERROR 3 IN repertoire.R:\nCODE HAS TO BE MODIFIED\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/ig_clustering OR REPORT AT gael.millot@pasteur.fr\n\n============\n\n")
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
    "lubridate", 
    "ggplot2"
)
# for(i in 1:length(req.package.list)){suppressMessages(library(req.package.list[i], character.only = TRUE))}
fun_pack(req.package = req.package.list, load = TRUE, lib.path = NULL) # packages are imported even if inside functions are written as package.name::function() in the present code


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
    "igblast_database_path", 
    "file_assembly_ch", 
    "repertoire_names_ch", 
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

arg.check2 <- NULL #
text.check2 <- NULL #
checked.arg.names2 <- NULL # for function debbuging: used by r_debugging_tools
ee <- expression(arg.check2 <- c(arg.check2, tempo$problem) , text.check2 <- c(text.check2, tempo$text) , checked.arg.names2 <- c(checked.arg.names2, tempo$object.name))

if( ! grepl(igblast_database_path, pattern = "^.+vdj$")){ # positive prop
    tempo.cat <- paste0("ERROR IN repertoire.R:\nRIGHT NOW, THE PIPELINE CAN ONLY WORK WITH THE VDJ DATASET OF THE imgt DATABASE (igblast_database_path PARAMETER FINISHING BY \"vdj\")\nHERE IT IS: \n", igblast_database_path)
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}

rep_file_names <- strsplit(repertoire_names_ch, split = " ")[[1]]
if(length(rep_file_names) == 0){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 4 IN repertoire.R:\nPROBLEM WITH repertoire_names_ch: ", paste(repertoire_names_ch, collapse = " "), "PLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/ig_clustering OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
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

################ internal variables

var1 <- c("v_call", "j_call") # names of the columns to deal with

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


df <- read.table(file_assembly_ch, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

alleles <- vector(mode = "list", length = length(rep_file_names))
for(i1 in 1:length(alleles)){
    tempo <- strsplit(rep_file_names[i1], split = "[_.]")[[1]]
    names(alleles)[i1] <- tempo[length(tempo) - 1]
    alleles[[i1]] <- scan(rep_file_names[i1], what = "character", quiet = TRUE)
}



################ End Data import


################ data verification

if( ! all(var1 %in% c("v_call", "j_call"))){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 5 IN repertoire.R:\nPROBLEM WITH THE var1 INTERNAL VARIABLE THAT MUST BE \"v_call\" AND \"j_call\"\nHERE IT IS:\n", paste(var1, collapse = "\n"), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/ig_clustering OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
}


if( ! all(var1 %in% names(df))){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 6 IN repertoire.R:\nPROBLEM WITH THE NAMES OF file_assembly_ch THAT MUST CONTAIN \"v_call\" AND \"j_call\"\nHERE IT IS:\n", paste(names(df), collapse = "\n"), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/ig_clustering OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
}

for(i1 in 1:length(alleles)){
    if(is.null(names(alleles)[i1])){
        stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 7 IN repertoire.R:\nPROBLEM WITH rep_file_names:\n", paste(rep_file_names, collapse = "\n"), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/ig_clustering OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
    }
}

################ end data verification


################ data modification and saving

# first loop with v_call and second with j_call
allele.kind <- tolower(substring(names(alleles), nchar(names(alleles)))) # get V from IGHV
for(i0 in 1:length(var1)){
    tempo2 <- tolower(substr(var1[i0], 1, 1))
    tempo.log <- allele.kind == tempo2
    if(any(is.na(tempo.log)) | sum(tempo.log, na.rm = TRUE) != 1){
        stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 8 IN repertoire.R:\nPROBLEM WITH tempo.log:\n", paste(tempo.log, collapse = "\n"), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/ig_clustering OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
    }else{
        tempo.pos <- which(tempo.log)
        df[ , names(df) == var1[i0]] <- factor(df[ , names(df) == var1[i0]], levels = alleles[[tempo.pos]])
        tempo.table <- table(df[ , names(df) == var1[i0]])
        write.table(tempo.table, file = paste0("./rep_", names(alleles)[tempo.pos], ".tsv"), row.names = FALSE, col.names = FALSE, sep = "\t") # separate repertoires
        # plot
        tempo.table.gg <- data.frame(as.data.frame(tempo.table), X = 1)
        tempo.gg.name <- "gg.indiv.plot."
        tempo.gg.count <- 0
        assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::ggplot(
            data = tempo.table.gg,
             aes(x = X, y = var1, fill= Freq)
        ))
        assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::geom_tile())
        assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::scale_fill_gradient(low="white", high="blue"))
        # assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), hrbrthemes::theme_ipsum())
        final.plot <- eval(parse(text = paste(paste0(tempo.gg.name, 1:tempo.gg.count), collapse = " + ")))
        ggplot2::ggsave(filename = paste0(names(alleles)[tempo.pos], ".png"), plot = final.plot, device = "png", path = ".", width = 5, height = 5, units = "in", dpi = 300)
        ggplot2::ggsave(filename = paste0(names(alleles)[tempo.pos], ".svg"), plot = final.plot, device = "svg", path = ".", width = 5, height = 5, units = "in", dpi = 300)
        ggplot2::ggsave(filename = paste0(names(alleles)[tempo.pos], ".pdf"), plot = final.plot, device = "pdf", path = ".", width = 5, height = 5, units = "in", dpi = 300)
    }
}


################ end data modification and saving

################ saving


# combined repertoires
for(i0 in 1:(length(var1) - 1)){
    for(i1 in 2:length(var1)){
        tempo1 <- tolower(substr(var1[i0], 1, 1))
        tempo2 <- tolower(substr(var1[i1], 1, 1))
        tempo.log1 <- allele.kind == tempo1
        tempo.log2 <- allele.kind == tempo2
        if(any(is.na(tempo.log1)) | sum(tempo.log1, na.rm = TRUE) != 1 | any(is.na(tempo.log2)) | sum(tempo.log2, na.rm = TRUE) != 1){
            stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 9 IN repertoire.R:\nPROBLEM WITH tempo.log1:\n", paste(tempo.log1, collapse = "\n"), " OR tempo.log2:\n", paste(tempo.log2, collapse = "\n"), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/ig_clustering OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
        }else{
            tempo.pos1 <- which(tempo.log1)
            tempo.pos2 <- which(tempo.log2)
            write.table(table(df[names(df) %in% c(var1[i0], var1[i1])]), file = paste0("./rep_", names(alleles)[tempo.pos1], "x", names(alleles)[tempo.pos2], ".tsv"), row.names = TRUE, col.names = NA, sep = "\t") # separate repertoires
        }
    }
}


################ end saving



################ Plotting







################ end Plotting tree


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
