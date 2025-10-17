#!/usr/bin/env Rscript

#########################################################################
##                                                                     ##
##     Tsv2fastaGff.R                                                     ##
##                                                                     ##
##     Gael A. Millot                                                  ##
##     Chloe Taurel                                                    ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Computational Biology Department                                ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################




################################ Aim


# The original aime of this program was to create all the fasta files from a .xlsx file.
# It has been slightly modified to create fasta files from a .tsv file, instead of a .xlsx file.
# It also now creates a gff file from the tsv one containing region coordinates. Goalign tool takes feature information, but only from a gff format.
# Note that for the creation of the .gff the tsv file must contain following columns => sequence_id, fr1_start, fr1_end, cdr1_start, cdr1_end, fr2_start, fr2_end, cdr2_start, cdr2_end, fr3_start, fr3_end, cdr3_start, cdr3_end

################################ End Aim


################################ Introduction


################################ End Introduction


################################ Acknowlegments


################################ End Acknowlegments


################################ Initialization


# R version checking
if(version$version.string != "R version 4.3.1 (2023-06-16 ucrt)"){
    cat(paste0("\n\nWARNING: THE ", version$version.string, " IS NOT THE 4.3.1 (2023-06-16 ucrt) RECOMMANDED\n\n"))
}
# other initializations
erase.objects = TRUE # write TRUE to erase all the existing objects in R before starting the algorithm and FALSE otherwise. Beginners should use TRUE
if(erase.objects == TRUE){
    rm(list = ls(all.names = TRUE))
    erase.objects = TRUE
}
erase.graphs = TRUE # write TRUE to erase all the graphic windows in R before starting the algorithm and FALSE otherwise
script <- "Tsv2fastaGff"
#cute <- "https://gitlab.pasteur.fr/gmillot/cute_little_R_functions/-/raw/v11.4.0/cute_little_R_functions.R" # single character string indicating the path of the cute_little_R_functions.R file required for this script. Example: cute <- "https://gitlab.pasteur.fr/gmillot/cute_little_R_functions/-/raw/v11.4.0/cute_little_R_functions.R"
#log <- "xlsx2fasta.log" # single character string indicating the name of the log file. Example: log <- "xlsx2fasta.log"


################################ End Initialization



################################ EXAMPLES FOR TEST AND EXPLANATION OF ARGUMENTS

# script <- "Tsv2fastaGff"

### Arguments : 

# path <- "C:/Users/gmillot/Documents/Git_projects/repertoire_profiler/work/b6/fe498b4e1e2a9cbdaa89cb9685d6cc/10_productive_seq_clone-pass_germ-pass_germ-seq-trans_germ-pass_shm-pass.tsv"      # tsv file containing data. needs to have all columns in Name, Seq and Germline
# Name <- "sequence_id"                # name of the column containing the sequence ids
# Seq <- "trimmed"        # name of the columns containing the sequences to put in the fasta file (can be a single string or several strings seperated by "," if several columns are needed. the fastas will then be created in different folders)
# clone_nb_seq <- 3                    # Minimum number of rows in the tsv file. The program expects this to be respected, otherwise raises an error.
# cute <- "C:/Users/gmillot/Documents/Git_projects/repertoire_profiler/bin/cute_little_R_functions_v12.8.R"
# seq_kind <- "CLONE"
# log <- "Tsv2fastaGff.log"


################################# End test




################################ Config import


tempo.cat <- "KIND OF RUN (SCRIPT, COPY-PASTE OR SOURCE): "
if(interactive() == FALSE){ # if(grepl(x = commandArgs(trailingOnly = FALSE), pattern = "R\\.exe$|\\/R$|Rcmd\\.exe$|Rcmd$|Rgui\\.exe$|Rgui$|Rscript\\.exe$|Rscript$|Rterm\\.exe$|Rterm$")){ # detection of script usage
    run.way <- "SCRIPT"
    cat(paste0("\n\n", tempo.cat, run.way, "\n\n"))
    command <- paste0(commandArgs(trailingOnly = FALSE), collapse = ",") # recover the full command
    args <- commandArgs(trailingOnly = TRUE) # recover arguments written after the call of the R script
    if(any(is.na(args))){
        stop(paste0("\n\n================\n\nERROR IN ", script, ".R\nTHE args OBJECT HAS NA\n\n================\n\n"), call. = FALSE)
    }
    tempo.arg.names <- c(
        "path", 
        "Name", 
        "Seq",
        "clone_nb_seq",
        "cute", 
        "seq_kind", 
        "log"
    ) # objects names exactly in the same order as in the bash code and recovered in args. Here only one, because only the path of the config file to indicate after the xlsx2fasta.R script execution
    if(length(args) != length(tempo.arg.names)){
        stop(paste0("\n\n================\n\nERROR IN ", script, ".R\nTHE NUMBER OF ELEMENTS IN args (", length(args),") IS DIFFERENT FROM THE NUMBER OF ELEMENTS IN tempo.arg.names (", length(tempo.arg.names),")\nargs:", paste0(args, collapse = ","), "\ntempo.arg.names:", paste0(tempo.arg.names, collapse = ","), "\n\n================\n\n"), call. = FALSE)
    }
    for(i1 in 1:length(tempo.arg.names)){
        assign(tempo.arg.names[i1], args[i1])
    }
    rm(tempo.arg.names, args, i1)
}else if(sys.nframe() == 0L){ # detection of copy-paste/direct execution (for debugging). With script it is also 0, with source, it is 4
    run.way <- "COPY-PASTE"
    cat(paste0("\n\n", tempo.cat, run.way, "\n\n"))
}else{
    run.way <- "SOURCE" # using source(), sys.nframe() is 4
    cat(paste0("\n\n", tempo.cat, run.way, "\n\n"))
}
rm(tempo.cat)


################################ End Config import

################################ Test

# No need because the parameter settings are not in another file but above

################################ end Test

################################ Recording of the initial parameters


param.list <- c(
    "erase.objects", 
    "erase.graphs", 
    "script", 
    "run.way",
    if(run.way == "SCRIPT"){"command"}, 
    "path", 
    "Name", 
    "Seq",
    "clone_nb_seq",
    "cute", 
    "seq_kind", 
    "log"
)
if(any(duplicated(param.list))){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 1 IN ", script, ".R\nTHE param.list OBJECT CONTAINS DUPLICATED ELEMENTS:\n", paste(param.list[duplicated(param.list)], collapse = " "), "\n\n================\n\n"), call. = FALSE) # message for developers
}
if(erase.objects == TRUE){
    created.object.control <- ls()[ ! ls() %in% "param.list"]
    if( ! (all(created.object.control %in% param.list) & all(param.list %in% created.object.control))){
        stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 2 IN ", script, ".R\nINCONSISTENCIES BETWEEN THE ARGUMENTS USED AND THE PARAMETERS REQUIRED IN THE EXECUTABLE CODE FILE\nTHE ARGUMENTS NOT PRESENT IN THE EXECUTABLE FILE (", script, ".R) ARE:\n", paste(created.object.control[ ! created.object.control %in% param.list], collapse = " "), "\nTHE PARAMETERS OF THE EXECUTABLE FILE (", script, ".R) NOT PRESENT IN THE ARGUMENTS ARE:\n", paste(param.list[ ! param.list %in% created.object.control], collapse = " "), "\n\n================\n\n"), call. = FALSE) # message for developers
    }
}
char.length <- nchar(param.list)
space.add <- max(char.length) - char.length + 5
param.ini.settings <- character(length = length(param.list))
for(i in 1:length(param.list)){
    param.ini.settings[i] <- paste0("\n", param.list[i], paste0(rep(" ", space.add[i]), collapse = ""), paste0(get(param.list[i]), collapse = ",")) # no env = sys.nframe(), inherit = FALSE in get() because look for function in the classical scope
}
clone_nb_seq <- as.numeric(clone_nb_seq) # numeric string already checked by nextflow


################################ End Recording of the initial parameters


################################ Functions


# Functions are built such that they should have no direct use of Global objects (going through the R scope), and only use function arguments
# 1) Cute little function is sourced for the moment into the .GlobalEnv environment, but may be interesting to put it into a new environement just above .GlobalEnv environment. See https://stackoverflow.com/questions/9002544/how-to-add-functions-in-an-existing-environment
# 2) Argument names of each function must not be a name of Global objects (error message otherwise)
# 3) Argument name of each function ends with "_fun" in the first function, "_2fun" in the second, etc. This prevent conflicts with the argument partial names when using these functions, notably when they are imbricated

fun_source_test <- function(path, script){ # do not write script = script: can produce recursive error if script argument is not specified thenafter
# AIM
# Test if one path exists (url or local)
# ARGUMENTS
# path: single character string of the path to test
# script: single character string of the current script file
# RETURN
# An error message if the path does not exists, nothing otherwise
# EXAMPLES
# fun_source_test(path = "caca", script = "test")
# DEBUGGING
# path = "https://zenodo.org/records/10814482/files/ig_sequences.xlsx" ; script = "test"
    name <- deparse(substitute(path))
    if(length(path) != 1){
        stop(paste0("\n\n============\n\nERROR IN ", script, ".R\n", name, " PARAMETER MUST BE LENGTH 1: ", paste(path, collapse = " "), "\n\n============\n\n"), call. = FALSE)
    }else if(grepl(x = path, pattern = "^http")){
        tempo.name <- paste0("tmp_xlsx2fasta.R_", as.integer(Sys.time()))
        if(file.exists(tempo.name)){
            stop(paste0("\n\n============\n\nERROR IN ", script, ".R\nTHE TEMPORARY FILE USED BY THE ", script, " SCRIPT ALREADY EXISTS:\n", file.path(tempo.name), "\n\n. PLEASE, RERUN.\n\n============\n\n"), call. = FALSE)
        }else{
            tempo.try <- try(suppressWarnings(suppressMessages(download.file(path, tempo.name, method="auto", quiet=TRUE))), silent = TRUE)
            if(file.exists(tempo.name)){
                file.remove(tempo.name)
            }
            if(any(grepl(x = tempo.try, pattern = "^[Ee]rror"))){
                stop(paste0("\n\n============\n\nERROR IN ", script, ".R\nHTTP INDICATED IN THE ", name, " PARAMETER DOES NOT EXISTS: ", path, "\n\n============\n\n"), call. = FALSE)
            }
        }
    }else if( ! grepl(x = path, pattern = "^http")){
        if( ! file.exists(path)){
            stop(paste0("\n\n============\n\nERROR IN ", script, ".R\nFILE INDICATED IN THE ", name, " PARAMETER DOES NOT EXISTS: ", path, "\n\n============\n\n"), call. = FALSE)
        }
    }else{
        tempo.cat <- paste0("\n\n================\n\nINTERNAL CODE ERROR 3 IN ", script, ".R: CODE HAS TO BE MODIFIED\n\n============\n\n")
        stop(tempo.cat, call. = FALSE)
    }
}



################ import functions from cute little functions toolbox

fun_source_test(path = cute, script = script)
source(cute, local = .GlobalEnv)

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
    tempo.cat <- paste0("ERROR IN ", script, ".R\nREQUIRED cute FUNCTION", ifelse(length(tempo) > 1, "S ARE", " IS"), " MISSING IN THE R ENVIRONMENT:\n", paste0(tempo, collapse = "()\n"))
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end required function checking

################ end import functions from cute little functions toolbox

################ local function: package import

################ end local function: package import

################ other functions

################ end other functions

################################ End Functions

################################ Ignition

ini.date <- Sys.time()
ini.time <- as.numeric(ini.date) # time of process begin, converted into seconds

################################ End ignition

################################ Checking


# reserved words
warn_file <- "warning.txt" # warning that will be displayed on the terminal

# end reserved words
# argument primary checking
arg.check <- NULL #
text.check <- NULL #
checked.arg.names <- NULL # for function debbuging: used by r_debugging_tools
ee <- expression(arg.check <- c(arg.check, tempo$problem) , text.check <- c(text.check, tempo$text) , checked.arg.names <- c(checked.arg.names, tempo$object.name))
tempo <- fun_check(data = path, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = Name, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = Seq, options = c("query", "igblast_full", "trimmed", "fwr1", "fwr2", "fwr3", "fwr4", "cdr1", "cdr2", "cdr3", "junction", "sequence_alignment", "v_sequence_alignment", "d_sequence_alignment", "j_sequence_alignment", "c_sequence_alignment", "germline_alignment", "v_germline_alignment", "d_germline_alignment", "j_germline_alignment", "c_germline_alignment")) ; eval(ee)
tempo <- fun_check(data = clone_nb_seq, class = "vector", typeof = "integer", length = 1, double.as.integer.allowed = TRUE) ; eval(ee)
tempo <- fun_check(data = seq_kind, options = c("ALL", "CLONE")) ; eval(ee)
# cute already tested above
tempo <- fun_check(data = log, class = "vector", typeof = "character", length = 1) ; eval(ee)
if(any(arg.check) == TRUE){ # normally no NA
    stop(paste0("\n\n================\n\n", paste(text.check[arg.check], collapse = "\n"), "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between == #
}

if(Seq == "query"){Seq2 <- c("sequence", "query_sequence_aa")}
if(Seq == "igblast_full"){Seq2 <- c("sequence", "sequence_aa")}
if(Seq == "trimmed"){Seq2 <- c("trimmed_sequence", "trimmed_sequence_aa")}
if(Seq == "fwr1"){Seq2 <- c("fwr1", "fwr1_aa")}
if(Seq == "fwr2"){Seq2 <- c("fwr2", "fwr2_aa")}
if(Seq == "fwr3"){Seq2 <- c("fwr3", "fwr3_aa")}
if(Seq == "fwr4"){Seq2 <- c("fwr4", "fwr4_aa")}
if(Seq == "cdr1"){Seq2 <- c("cdr1", "cdr1_aa")}
if(Seq == "cdr2"){Seq2 <- c("cdr2", "cdr2_aa")}
if(Seq == "cdr3"){Seq2 <- c("cdr3", "cdr3_aa")}
if(Seq == "junction"){Seq2 <- c("junction", "junction_aa")}
if(Seq == "sequence_alignment"){Seq2 <- c("sequence_alignment", "sequence_alignment_aa")}
if(Seq == "v_sequence_alignment"){Seq2 <- c("v_sequence_alignment", "v_sequence_alignment_aa")}
if(Seq == "d_sequence_alignment"){Seq2 <- c("d_sequence_alignment", "d_sequence_alignment_aa")}
if(Seq == "j_sequence_alignment"){Seq2 <- c("j_sequence_alignment", "j_sequence_alignment_aa")}
if(Seq == "c_sequence_alignment"){Seq2 <- c("c_sequence_alignment", "c_sequence_alignment_aa")}
if(Seq == "germline_alignment"){Seq2 <- c("germline_alignment", "germline_alignment_aa")}
if(Seq == "v_germline_alignment"){Seq2 <- c("v_germline_alignment", "v_germline_alignment_aa")}
if(Seq == "d_germline_alignment"){Seq2 <- c("d_germline_alignment", "d_germline_alignment_aa")}
if(Seq == "j_germline_alignment"){Seq2 <- c("j_germline_alignment", "j_germline_alignment_aa")}
if(Seq == "c_germline_alignment"){Seq2 <- c("c_germline_alignment", "c_germline_alignment_aa")}

Germline <- NULL
if(base::grepl(x = Seq, pattern = "sequence_alignment") & seq_kind == "CLONE" ){
    Germline <- base::sub(x = Seq2, pattern = "sequence", replacement = "germline")
}else if(Seq %in% c("query", "igblast_full", "trimmed") & seq_kind == "CLONE"){
    Germline <- c("germline_alignment", "germline_alignment_aa")
}
# end argument primary checking
# second round of checking and data preparation
# management of NA arguments
# end management of NA arguments
# management of NULL arguments
tempo.arg <-c(
    "path", 
    "cute", 
    "log"
)
tempo.log <- sapply(lapply(tempo.arg, FUN = get, env = sys.nframe(), inherit = FALSE), FUN = is.null)
if(any(tempo.log) == TRUE){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN ", script, ".R:\n", ifelse(sum(tempo.log, na.rm = TRUE) > 1, "THESE ARGUMENTS\n", "THIS ARGUMENT\n"), paste0(tempo.arg[tempo.log], collapse = "\n"),"\nCANNOT BE NULL")
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end management of NULL arguments
# management of ""
tempo.arg <-c(
    "path", 
    "cute", 
    "log"
)
tempo.log <- sapply(lapply(tempo.arg, FUN = get, env = sys.nframe(), inherit = FALSE), FUN = function(x){any(x == "")})
if(any(tempo.log) == TRUE){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN ", script, ".R:\n", ifelse(sum(tempo.log, na.rm = TRUE) > 1, "THESE ARGUMENTS\n", "THIS ARGUMENT\n"), paste0(tempo.arg[tempo.log], collapse = "\n"),"\nCANNOT BE \"\"")
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

fun_source_test(path = path, script = script)



# end other checkings
# reserved word checking
# end reserved word checking
# end second round of checking and data preparation
# package checking
# end package checking


################################ End Checking


################################ Main code


################ Ignition


fun_report(data = paste0("\n\n################################################################ ", script, ".R SCRIPT\n\n"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0("\n\n################################ RUNNING DATE AND STARTING TIME"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0(ini.date, "\n\n"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0("\n\n################################ RUNNING"), output = log, path = "./", overwrite = FALSE)


################ End ignition


################ Graphical parameter initialization


################ End graphical parameter initialization


################ Data import


obs <- read.table(path, header = TRUE, sep = "\t")


################ end Data import


############ check

# Inactivated because the tsv input file can be empty, there just won't be an output if so
# if(length(obs) == 0 || nrow(obs) == 0){
#    stop(paste0("\n\n============\n\nERROR IN ", script, ".R\nFILE INDICATED IN THE path PARAMETER IS EMPTY:\n", path, "\n\n============\n\n"), call. = FALSE)
#}

if( ! Name %in% names(obs)){
    stop(paste0("\n\n============\n\nERROR IN ", script, ".R\n\nTHE Name PARAMETER MUST BE A COLUMN NAME OF THE IMPORTED FILE:\n", path, "\n\nHERE IT IS Name:\n", Name, "\n\nCOLUMN NAMES:\n", paste(names(obs), collapse = "\n"), "\n\n============\n\n"), call. = FALSE)
}

if(any(duplicated(obs[, Name]))){
    stop(paste0("\n\n============\n\nERROR IN ", script, ".R\n\nDUPLICATED VALUE NOT AUTHORIZED IN THE COLUMN OF THE Name PARAMETER\n\nDUPLICATED VALUES ARE:\n", obs[duplicated(obs[, Name]), Name], "\n\nIN POSITIONS:\n", paste(which(obs[ , Name] %in% obs[duplicated(obs[, Name]), Name]), collapse = "\n"), "\n\n============\n\n"), call. = FALSE)
}

if( ! all(Seq2 %in% names(obs))){
    stop(paste0("\n\n============\n\nERROR IN ", script, ".R\n\nTHE Seq2 PARAMETER MUST BE COLUMN NAMES OF THE IMPORTED FILE:\n", path, "\n\nHERE IT IS Seq:\n", paste(Seq2, collapse = "\n"), "\n\nCOLUMN NAMES:\n", paste(names(obs), collapse = "\n"), "\n\n============\n\n"), call. = FALSE)
}

tempo.log <- is.na(obs[ , Name]) | obs[ , Name] == ""
if(any(tempo.log)){
    stop(paste0("\n\n============\n\nERROR IN ", script, ".R\nIMPORTED FILE:\n", path, "\nHAS AN EMPTY CELL IN THE ", Name, " COLUMN IN LINES:\n", paste(which(tempo.log), collapse = "\n"), "\n\n============\n\n"), call. = FALSE)
}

for(i0 in names(obs)){ # NA in xlsx file become "NA". Thus, has to be replaced by NA
    tempo.log <- obs[ , i0] == "NA" & ! is.na(obs[ , i0])
    if(any(tempo.log, na.rm = TRUE)){
        obs[tempo.log, i0] <- NA
    }
}


############ end check


############ main

## Create the fasta files :

count = 0
multiple_v_genes <- FALSE
multiple_j_genes <- FALSE

for(i0 in Seq2){
    count = count + 1
    tempo.log <- is.na(obs[ , i0]) | obs[ , i0] == ""
    if(sum(!tempo.log, na.rm = TRUE) >= clone_nb_seq){
        # Only create fasta files with at least <clone_nb_seq> sequences (Minimun number of non-identical sequences per clonal group for tree plotting)
        # NB : clone_nb_seq is defined in nextflow.config
        if(any(tempo.log)){
            tempo.warn <- paste0("IMPORTED FILE:\n", path, "\nHAS ", sum(tempo.log, na.rm = TRUE), " AMONG ", nrow(obs), " EMPTY SEQUENCES (NA OR \"\") IN THE ", i0, " COLUMN IN LINES:\n", paste(which(tempo.log), collapse = "\n"))
            cat(paste0("\nWARNING IN ", script, ".R\n", tempo.warn, "\n\n"))
            fun_report(data = paste0("WARNING\n", tempo.warn), output = log, path = "./", overwrite = FALSE)
            warn <- paste0(ifelse(is.null(warn), tempo.warn, paste0(warn, "\n\n", tempo.warn)))
        }else{
            tempo.cat <- paste0("\nIMPORTED FILE:\n", path, "\nHAS NO EMPTY SEQUENCES (NA OR \"\") AMONG ", nrow(obs), " IN THE ", i0, " COLUMN\n")
            fun_report(data = tempo.cat, output = log, path = "./", overwrite = FALSE)
        }
        tempo.df <- obs[ ! tempo.log, ]
        # Sequences in the same clone_assigned_seq.tsv file belong to the same clonal group and should have the same values in columns relative to clonal groups
        if(any(tempo.df$clone_id != tempo.df$clone_id[1])){
            stop(paste0("\n\n================\n\nERROR IN ", script, ".R\nALL clone_id VALUES SHOULD BE THE SAME IN A clone_assigned_seq.tsv FILE, BUT THEY ARE NOT.\nHERE THEY ARE : ", paste0(tempo.df$clone_id, collapse = "\n"),"\n\n================\n\n"), call. = FALSE)
        }
        if(seq_kind == "ALL"){
            multiple_v_genes <- TRUE
        }else if(any(tempo.df$v_gene != tempo.df$v_gene[1])) {
            tempo.warn <- paste0(
                "Different genes for the V cassette were found in this clonal group\n",
                "Here they are:\n", paste0(tempo.df$v_gene, collapse = "\n"),
                "\nSet clone_strategy = \"first\" in nextflow.config to avoid this."
            )
            cat(paste0("\nWARNING IN ", script, ".R\n", tempo.warn, "\n\n"))
            fun_report(data = paste0("WARNING\n", tempo.warn), output = log, path = "./", overwrite = FALSE)
            warn <- paste0(ifelse(is.null(warn), tempo.warn, paste0(warn, "\n\n", tempo.warn)))
            multiple_v_genes <- TRUE
        }
        if(seq_kind == "ALL"){
            multiple_j_genes <- TRUE
        }else if(any(tempo.df$j_gene != tempo.df$j_gene[1])) {
             tempo.warn <- paste0(
                "Different genes for the V cassette were found this clonal group\n",
                "Here they are:\n", paste0(tempo.df$j_gene, collapse = "\n"),
                "\nSet clone_strategy = \"first\" in nextflow.config to avoid this."
            )
            cat(paste0("\nWARNING IN ", script, ".R\n", tempo.warn, "\n\n"))
            fun_report(data = paste0("WARNING\n", tempo.warn), output = log, path = "./", overwrite = FALSE)
            warn <- paste0(ifelse(is.null(warn), tempo.warn, paste0(warn, "\n\n", tempo.warn)))
            multiple_j_genes <- TRUE
        }
        # if(any(tempo.df$v_gene != tempo.df$v_gene[1])){
        #     stop(paste0("\n\n================\n\nERROR IN ", script, ".R\nALL v_gene VALUES SHOULD BE THE SAME IN A clone_assigned_seq.tsv FILE, BUT THEY ARE NOT.\nHERE THEY ARE : ", paste0(tempo.df$v_gene, collapse = "\n"),"\n\n================\n\n"), call. = FALSE)
        # }
        # if(any(tempo.df$j_gene != tempo.df$j_gene[1])){
        #     stop(paste0("\n\n================\n\nERROR IN ", script, ".R\nALL j_gene VALUES SHOULD BE THE SAME IN A clone_assigned_seq.tsv FILE, BUT THEY ARE NOT.\nHERE THEY ARE : ", paste0(tempo.df$j_gene, collapse = "\n"),"\n\n================\n\n"), call. = FALSE)
        # }
        if(seq_kind == "CLONE"){
            if(any(tempo.df$junction_length != tempo.df$junction_length[1])){
                stop(paste0("\n\n================\n\nERROR IN ", script, ".R\nALL junction_length VALUES SHOULD BE THE SAME IN A clone_assigned_seq.tsv FILE, BUT THEY ARE NOT.\nHERE THEY ARE : ", paste0(tempo.df$junction_length, collapse = "\n"),"\n\n================\n\n"), call. = FALSE)
            }
        }
        # End check of columns relative to clonal groups

        # Creation of the fasta file
        if(base::grepl(x = i0, pattern = "_aa$")){
            dir_name <- paste0("./for_alignment_aa")
        }else if(count > 2){
            stop(paste0("\n\n================\n\nERROR IN ", script, ".R\ncount CANNOT BE MORE THAN 2 : ", paste0(count, collapse = "\n"),"\n\n================\n\n"), call. = FALSE)
        }else{
            dir_name <- paste0("./for_alignment_nuc")
        }
        dir.create(dir_name, showWarnings = FALSE, recursive = TRUE)

        # Extract all different values of v_genes and j_genes in the clonal group (different genes for a same cassette can be present in a clonal group if the clone_strategy = "set" in nextflow.config)
        # If several genes in a column, they are separated by a comma (",")
        v_genes_all <- unique(unlist(strsplit(tempo.df$v_gene, ",")))
        j_genes_all <- unique(unlist(strsplit(tempo.df$j_gene, ",")))
        # Replace the "," characters by "-" to have appropriate file names
        v_gene_clean <- paste(sort(v_genes_all), collapse = "-") 
        j_gene_clean <- paste(sort(j_genes_all), collapse = "-")

        if(seq_kind == "CLONE"){
            tempo.name <- paste0(i0, "_clone_id_", tempo.df[1, "clone_id"], "_", v_gene_clean, "_", j_gene_clean, ".fasta") # Create the name of the file
        }else{
            tempo.name <- paste0(i0, ".fasta")
        }
        for(i1 in 1:nrow(tempo.df)) {
            tempo.cat <- paste0(">", tempo.df[i1, Name], "\n", tempo.df[i1, i0], "\n")
            cat(tempo.cat, file = file.path(dir_name, tempo.name), append = TRUE)
        }

        # Germline addition to the fasta
        # Sort the dataframe by the column specified in Name
        if( ! base::is.null(Germline)){
            germ = Germline[count]
            tempo.df <- tempo.df[order(tempo.df[[Name]]), ]

            # Get the germline values and compute frequencies
            germ_values <- tempo.df[[germ]]
            germ_table <- table(germ_values)
            max_freq <- max(germ_table)
            most_frequent_germs <- names(germ_table[germ_table == max_freq])

            # Choose the first matching germline by order in tempo.df (sorted by Name)
            selected_index <- which(tempo.df[[germ]] %in% most_frequent_germs)[1]
            selected_germline <- tempo.df[selected_index, germ]

            # Get cleaned v_gene and j_gene values for the chosen germline row
            germ_v_gene <- gsub(",", "-", tempo.df[selected_index, "germline_v_gene"])
            germ_j_gene <- gsub(",", "-", tempo.df[selected_index, "germline_j_gene"])

            # Write germline to the same fasta file as above
            germ_seq_name <- paste0(Germline[1], "_clone_id_", tempo.df[1, "clone_id"], "_", germ_v_gene, "_", germ_j_gene) # no choice to use Germline[1], otherwise goalign codonalign cannot align (because nuc and aa must have the same name)
            tempo.cat <- paste0(">", germ_seq_name, "\n", selected_germline, "\n")
            cat(tempo.cat, file = file.path(dir_name, tempo.name), append = TRUE)

            # If germlines differ in the group, issue a warning
            if (length(unique(germ_values)) > 1) {
                tempo.warn <- paste0(
                    "Multiple germline sequences found for clone_id ", tempo.df[1, "clone_id"], ".\n",
                    "Selected germline from sequence: ", tempo.df[selected_index, Name], "\n",
                    "Associated germline_v_gene: ", tempo.df[selected_index, "germline_v_gene"], "\n",
                    "Associated germline_j_gene: ", tempo.df[selected_index, "germline_j_gene"], "\n",
                    "Associated v_gene: ", tempo.df[selected_index, "v_gene"], "\n",
                    "Associated j_gene: ", tempo.df[selected_index, "j_gene"], "\n",
                    "Used germline sequence: ", selected_germline
                )
                cat(paste0("\nWARNING IN ", script, ".R\n", tempo.warn, "\n\n"))
                fun_report(data = paste0("WARNING\n", tempo.warn), output = log, path = "./", overwrite = FALSE)
                warn <- paste0(ifelse(is.null(warn), tempo.warn, paste0(warn, "\n\n", tempo.warn)))
            }
        }
    } else {
        stop(paste0("\n\n================\n\nERROR IN ", script, ".R\nNO FASTA FILE CREATED BECAUSE THE IMPORTED FILE:\n", path, "\nHAS MORE THAN ", clone_nb_seq, " EMPTY SEQUENCES (NA OR \"\") IN THE ", i0, " COLUMN\n\n================\n\n"))
    }

}

if(seq_kind == "ALL"){
    if (multiple_j_genes || multiple_v_genes){
        tempo.warn <- paste0(
            "\n[WARNING] Multiple V or J genes detected in a clonal group. ",
            "\nSet clone_strategy = \"first\" in nextflow.config to avoid this."
        )
        cat(paste0("\nWARNING IN ", script, ".R\n", tempo.warn, "\n\n"))
        fun_report(data = paste0("WARNING\n", tempo.warn), output = log, path = "./", overwrite = FALSE)
        warn <- paste0(ifelse(is.null(warn), tempo.warn, paste0(warn, "\n\n", tempo.warn)))
    }
}


## Enf of creating the fasta files

## Create the gff file
# https://shazam.readthedocs.io/en/stable/topics/setRegionBoundaries/
if( ! base::is.null(Germline)){
    seq_kind <- c("sequence", "germline")
    region_kind <- c("v", "d", "j", "c")
    seq_vdjc_features <- paste(region_kind, seq_kind[1], sep = "_")
    seq_vdjc_features_colors <- c("red", "green", "blue", "yellow")
    germline_vdjc_features <- paste(region_kind, seq_kind[2], sep = "_")
    germline_vdjc_features_colors <- c("red", "green", "blue", "yellow")
    fwr_cdr_features <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")
    fwr_cdr_features_colors <- c("yellow", "pink", "yellow", "pink", "yellow", "pink", "yellow")

    ## Convert all column names in obs to lowercase for comparison
    # colnames_lc <- tolower(colnames(obs))

    ## Check that all expected columns (case-sensitive) exist
    tempo_names <- c("seq_vdjc_features", "germline_vdjc_features", "fwr_cdr_features")
    missing_col <- character()
    non_unique_cols <- character()
    for (suf in c("_start", "_end")) {
        for(i1 in tempo_names){
            for(i2 in get(i1)) {
                col <- paste0(i2, suf)
                if (!(col %in% colnames(tempo.df))) {
                    missing_col <- c(missing_col, col)
                }else{
                    unique_vals <- unique(tempo.df[[col]])
                    if (length(unique_vals) != 1) {
                        non_unique_cols <- c(non_unique_cols, col)
                    }
                }
            }
        }
    }
    if (length(missing_col) > 0) {
        stop(paste0("\n\n================\n\nERROR IN ", script, ".R\nONE OR MORE COORDINATE COLUMNS MISSING FROM THE IMPORTED FILE:\n", path, "\nHERE IS THE MISSING COLUMN : ", paste0(missing_col, collapse = "\n"), "\n\n================\n\n"), call. = FALSE)
    }
    if (length(non_unique_cols) > 0) {
        tempo.warn <- paste0(
            "DATA ROWS IN IMPORTED FILE : ", path, "\n",
            "HAVE DIFFERENT VALUES FOR THE FOLLOWING REGION COORDINATES COLUMN.\n",
            "THE MOST FREQUENT VALUE WAS TAKEN FOR THESE COLUMNS:\n",
            paste(non_unique_cols, collapse = "\n")
        )
        cat(paste0("\nWARNING IN ", script, ".R\n", tempo.warn, "\n\n"))
        fun_report(data = paste0("WARNING\n", tempo.warn), output = log, path = "./", overwrite = FALSE)
        warn <- paste0(ifelse(is.null(warn), tempo.warn, paste0(warn, "\n\n", tempo.warn)))
    }
    ## end Check that all expected columns (case-sensitive) exist


    for(i0 in tempo_names){
        gff_rows <- list()
        gff_rows_convert <- list()
        for(i1 in 1:length(get(i0))){
            # Find actual column names (preserve original case)
            start_col <- colnames(tempo.df)[colnames(tempo.df) == paste0(get(i0)[i1], "_start")]
            end_col   <- colnames(tempo.df)[colnames(tempo.df) == paste0(get(i0)[i1], "_end")]

            # Take the most frequent value for coordinate if not unique, and handle the case when all values are NA in a column
            if (all(is.na(tempo.df[[start_col]]))) {
                start_val <- NA
            } else {
                start_val <- names(which.max(table(tempo.df[[start_col]])))[1] # [1] in case of equality in max()
            }
            if (all(is.na(tempo.df[[end_col]]))) {
                end_val <- NA
            } else {
                end_val <- names(which.max(table(tempo.df[[end_col]])))[1] # [1] in case of equality in max()
            }
            # Skip this feature if either coordinate is NA
            if (is.na(start_val) || is.na(end_val)) {
                tempo.warn <- paste0("Skipping feature ", get(i0)[i1], " in nuc GFF of ", i0, ", due to missing coordinates: start = ", start_val, ", end = ", end_val, "\n")
                cat(paste0("\nWARNING IN ", script, ".R\n", tempo.warn, "\n\n"))
                fun_report(data = paste0("WARNING\n", tempo.warn), output = log, path = "./", overwrite = FALSE)
                warn <- paste0(ifelse(is.null(warn), tempo.warn, paste0(warn, "\n\n", tempo.warn)))
            }else{
                seq_name <- paste0(Germline[1], "_clone_id_", tempo.df[1, "clone_id"], "_", germ_v_gene, "_", germ_j_gene)
                seq_name_nuc <- paste0(Seq2[1], "_clone_id_", tempo.df[1, "clone_id"], "_", germ_v_gene, "_", germ_j_gene)
                seq_name_aa <- paste0(Seq2[2], "_clone_id_", tempo.df[1, "clone_id"], "_", germ_v_gene, "_", germ_j_gene)
                tempo_start <- as.integer(start_val) # unique value
                tempo_end <- as.integer(end_val)   # unique value
                row <- c(
                    seq_name,
                    ".",
                    "gene",
                    tempo_start, 
                    tempo_end, 
                    ".",
                    ".",
                    ".",
                    paste0("Name=", get(i0)[i1], ";Color=", get(paste0(i0, "_colors"))[i1])
                )
                gff_rows[[length(gff_rows) + 1]] <- row
                # tol <- .Machine$double.eps^0.5 ; abs(tempo_start %% 3) < tol to use if tempo_start is a double
                if((tempo_start - 1) %% 3 == 0 && tempo_end %% 3 == 0){
                    row_aa <- c(
                        seq_name,
                        ".",
                        "gene",
                        as.integer((tempo_start - 1) / 3 + 1), # unique value
                        as.integer(tempo_end / 3),   # unique value
                        ".",
                        ".",
                        ".",
                        paste0("Name=", get(i0)[i1], ";Color=", get(paste0(i0, "_colors"))[i1])
                    )
                    gff_rows_convert[[length(gff_rows) + 1]] <- row_aa 
                }else{
                    tempo.warn <- paste0("Skipping feature ", get(i0)[i1], " in aa GFF of ", i0, ", due to due to coordinates not multiple of three : start = ", tempo_start, ", end = ", tempo_end, "\n")
                    cat(paste0("\nWARNING IN ", script, ".R\n", tempo.warn, "\n\n"))
                    fun_report(data = paste0("WARNING\n", tempo.warn), output = log, path = "./", overwrite = FALSE)
                    warn <- paste0(ifelse(is.null(warn), tempo.warn, paste0(warn, "\n\n", tempo.warn)))
                }
            }
        }
        # gff_rows is a list of GFF row vectors, as before.
        # Write to GFF file
        gff_table <- do.call(rbind, gff_rows)
        if(length(gff_table) > 0){
            gff_lines <- apply(gff_table, 1, function(x) paste(x, collapse="\t"))
        }else{
            gff_lines <- character()
        }
        gff_lines <- c("##gff-version 3", gff_lines)
        output_gff <- paste0(seq_name_nuc, "_", i0, "_nuc.gff")
        writeLines(gff_lines, con = output_gff)

        gff_table_convert <- do.call(rbind, gff_rows_convert)
        if(length(gff_table_convert) > 0){
            gff_lines_convert <- apply(gff_table_convert, 1, function(x) paste(x, collapse="\t"))
        }else{
            gff_lines_convert <- character()
        }
        gff_lines_convert <- c("##gff-version 3", gff_lines_convert)
        output_gff_convert <- paste0(seq_name_aa, "_", i0, "_aa.gff")
        writeLines(gff_lines_convert, con = output_gff_convert)
    }
}

## End of creating the gff file

fun_report(data = paste0("\n\n################################ RUNNING END"), output = log, path = "./", overwrite = FALSE)

############ end main


################ Seeding inactivation


################ end Seeding inactivation


################ Environment saving


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


fun_report(data = paste0("\n\n################################ INITIAL SETTINGS OF PARAMETERS"), output = log, path = "./", overwrite = FALSE, sep = 1)
fun_report(data = param.ini.settings, output = log, path = "./", overwrite = FALSE, , vector.cat = TRUE)
fun_report(data = paste0("\n\n################################ R SYSTEM AND PACKAGES"), output = log, path = "./", overwrite = FALSE)
tempo <- sessionInfo()
# tempo$otherPkgs <- tempo$otherPkgs[order(names(tempo$otherPkgs))] # sort the packages
tempo$loadedOnly <- tempo$loadedOnly[order(names(tempo$loadedOnly))] # sort the packages
fun_report(data = tempo, output = log, path = "./", overwrite = FALSE, , vector.cat = TRUE)
end.date <- Sys.time()
end.time <- as.numeric(end.date)
total.lapse <- round(lubridate::seconds_to_period(end.time - ini.time))
fun_report(data = paste0("\n\n################################ JOB END\n\n\nTIME: ", end.date, "\n\nTOTAL TIME LAPSE: ", total.lapse, "\n"), output = log, path = "./", overwrite = FALSE)


################ end Parameter printing


################################ End Main code




