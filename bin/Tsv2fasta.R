#!/usr/bin/env Rscript

#########################################################################
##                                                                     ##
##     Tsv2fasta.R                                                     ##
##                                                                     ##
##     Gael A. Millot                                                  ##
##     Chloe Taurel                                                    ##
##     Bioinformatics and Biostatistics Hub                            ##
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
script <- "Tsv2fasta"
#cute <- "https://gitlab.pasteur.fr/gmillot/cute_little_R_functions/-/raw/v11.4.0/cute_little_R_functions.R" # single character string indicating the path of the cute_little_R_functions.R file required for this script. Example: cute <- "https://gitlab.pasteur.fr/gmillot/cute_little_R_functions/-/raw/v11.4.0/cute_little_R_functions.R"
#log <- "xlsx2fasta.log" # single character string indicating the name of the log file. Example: log <- "xlsx2fasta.log"


################################ End Initialization



################################ EXAMPLES FOR TEST AND EXPLANATION OF ARGUMENTS

# script <- "Tsv2fastaGff"

### Arguments : 

# path <- "C:/Users/gmillot/Documents/Git_projects/repertoire_profiler/work/b6/fe498b4e1e2a9cbdaa89cb9685d6cc/10_productive_seq_clone-pass_germ-pass_germ-seq-trans_germ-pass_shm-pass.tsv"      # tsv file containing data. needs to have all columns in Name, align_seq and Germline
# Name <- "sequence_id"                # name of the column containing the sequence ids
# align_seq <- "sequence_alignment"        # name of the columns containing the sequences to put in the fasta file (can be a single string or several strings seperated by "," if several columns are needed. the fastas will then be created in different folders)
# clone_germline_kind <- "dmask"
# align_clone_nb <- 3                    # Minimum number of rows in the tsv file. The program expects this to be respected, otherwise raises an error.
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
        "Name", # sequence_id column
        "align_seq",
        "clone_germline_kind", 
        "align_clone_nb",
        "cute", 
        "seq_kind", # ALL or CLONE
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
    "align_seq", 
    "clone_germline_kind", 
    "align_clone_nb", 
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
align_clone_nb <- as.numeric(align_clone_nb) # numeric string already checked by nextflow


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
tempo <- fun_check(data = align_seq, options = c("query", "igblast_full", "trimmed", "fwr1", "fwr2", "fwr3", "fwr4", "cdr1", "cdr2", "cdr3", "junction", "sequence_alignment", "v_sequence_alignment", "d_sequence_alignment", "j_sequence_alignment", "c_sequence_alignment", "germline_alignment", "v_germline_alignment", "d_germline_alignment", "j_germline_alignment", "c_germline_alignment")) ; eval(ee)
tempo <- fun_check(data = clone_germline_kind, options = c("full","dmask","vonly")) ; eval(ee)
tempo <- fun_check(data = align_clone_nb, class = "vector", typeof = "integer", length = 1, double.as.integer.allowed = TRUE) ; eval(ee)
tempo <- fun_check(data = seq_kind, options = c("ALL", "CLONE", "IMGT")) ; eval(ee)
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
    "Name", 
    "align_seq",
    "clone_germline_kind", 
    "align_clone_nb",
    "cute", 
    "seq_kind", 
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
    "Name", 
    "align_seq",
    "clone_germline_kind", 
    "cute", 
    "seq_kind", 
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



# print the sequence_alignment_with_gaps sequences as aligned fasta
if(seq_kind == "IMGT"){
    tempo_name <- "sequence_alignment_with_gaps"
    tempo_name_aa <- "sequence_alignment_with_gaps_aa"
    if( ! tempo_name %in% names(obs)){
        stop(paste0("\n\n============\n\nERROR IN ", script, ".R\n\nsequence_alignment_with_gaps MUST BE COLUMN A COLUMN NAME OF THE IMPORTED FILE:\n", path, "\n\nCOLUMN NAMES:\n", paste(names(obs), collapse = "\n"), "\n\n============\n\n"), call. = FALSE)
    }
    if( ! tempo_name_aa %in% names(obs)){
        stop(paste0("\n\n============\n\nERROR IN ", script, ".R\n\nsequence_alignment_with_gaps_aa MUST BE COLUMN A COLUMN NAME OF THE IMPORTED FILE:\n", path, "\n\nCOLUMN NAMES:\n", paste(names(obs), collapse = "\n"), "\n\n============\n\n"), call. = FALSE)
    }
    # 1. Find the maximum sequence length
    # obs[[tempo_name]] <- sapply(X = obs[[tempo_name]], FUN = function(x){gsub(x = x, pattern = "-", replacement = "")}) # do not remove the hyphens already here because alignments are not good anymore
    max_len <- max(nchar(obs[[tempo_name]]))
    max_len_aa <- max(nchar(obs[[tempo_name_aa]]))
    # 2. Pad each sequence with hyphens at the end
    obs[[tempo_name]] <- sapply(X = obs[[tempo_name]], FUN = function(x){paste0(x, paste(rep("-", max_len - nchar(x)), collapse = ""))})
    obs[[tempo_name_aa]] <- sapply(X = obs[[tempo_name_aa]], FUN = function(x){paste0(x, paste(rep("-", max_len_aa - nchar(x)), collapse = ""))})
    for(i1 in 1:nrow(obs)){
        tempo_cat <- paste0(">", obs[i1, Name], "\n", obs[i1, tempo_name], "\n")
        cat(tempo_cat, file = file.path(paste0(tempo_name, "_imgt_nuc.fasta")), append = TRUE)
        tempo_cat <- paste0(">", obs[i1, Name], "\n", obs[i1, tempo_name_aa], "\n")
        cat(tempo_cat, file = file.path(paste0(tempo_name, "_imgt_aa.fasta")), append = TRUE) # tempo_name because name of nuc kept
    }
}else{
# end print the sequence_alignment_with_gaps sequences as aligned fasta
    if(align_seq == "query"){align_seq2 <- c("sequence", "query_sequence_aa")}
    if(align_seq == "igblast_full"){align_seq2 <- c("sequence", "sequence_aa")}
    if(align_seq == "trimmed"){align_seq2 <- c("trimmed_sequence", "trimmed_sequence_aa")}
    if(align_seq == "fwr1"){align_seq2 <- c("fwr1", "fwr1_aa")}
    if(align_seq == "fwr2"){align_seq2 <- c("fwr2", "fwr2_aa")}
    if(align_seq == "fwr3"){align_seq2 <- c("fwr3", "fwr3_aa")}
    if(align_seq == "fwr4"){align_seq2 <- c("fwr4", "fwr4_aa")}
    if(align_seq == "cdr1"){align_seq2 <- c("cdr1", "cdr1_aa")}
    if(align_seq == "cdr2"){align_seq2 <- c("cdr2", "cdr2_aa")}
    if(align_seq == "cdr3"){align_seq2 <- c("cdr3", "cdr3_aa")}
    if(align_seq == "junction"){align_seq2 <- c("junction", "junction_aa")}
    if(align_seq == "sequence_alignment"){align_seq2 <- c("sequence_alignment", "sequence_alignment_aa")}
    if(align_seq == "v_sequence_alignment"){align_seq2 <- c("v_sequence_alignment", "v_sequence_alignment_aa")}
    if(align_seq == "d_sequence_alignment"){align_seq2 <- c("d_sequence_alignment", "d_sequence_alignment_aa")}
    if(align_seq == "j_sequence_alignment"){align_seq2 <- c("j_sequence_alignment", "j_sequence_alignment_aa")}
    if(align_seq == "c_sequence_alignment"){align_seq2 <- c("c_sequence_alignment", "c_sequence_alignment_aa")}
    if(align_seq == "germline_alignment"){align_seq2 <- c("germline_alignment", "germline_alignment_aa")}
    if(align_seq == "v_germline_alignment"){align_seq2 <- c("v_germline_alignment", "v_germline_alignment_aa")}
    if(align_seq == "d_germline_alignment"){align_seq2 <- c("d_germline_alignment", "d_germline_alignment_aa")}
    if(align_seq == "j_germline_alignment"){align_seq2 <- c("j_germline_alignment", "j_germline_alignment_aa")}
    if(align_seq == "c_germline_alignment"){align_seq2 <- c("c_germline_alignment", "c_germline_alignment_aa")}

    if( ! all(align_seq2 %in% names(obs))){
        stop(paste0("\n\n============\n\nERROR IN ", script, ".R\n\nTHE align_seq2 PARAMETER MUST BE COLUMN NAMES OF THE IMPORTED FILE:\n", path, "\n\nHERE IT IS align_seq2:\n", paste(align_seq2, collapse = "\n"), "\n\nCOLUMN NAMES:\n", paste(names(obs), collapse = "\n"), "\n\n============\n\n"), call. = FALSE)
    }

    if(seq_kind == "CLONE"){
        Germline <- c("clonal_germline_sequence_no_gaps", "clonal_germline_sequence_aa")
        if( ! all(Germline %in% names(obs))){
            stop(paste0("\n\n============\n\nERROR IN ", script, ".R\n\nTHE align_seq2 PARAMETER MUST BE COLUMN NAMES OF THE IMPORTED FILE:\n", path, "\n\nHERE IT IS align_seq2:\n", paste(Germline, collapse = "\n"), "\n\nCOLUMN NAMES:\n", paste(names(obs), collapse = "\n"), "\n\n============\n\n"), call. = FALSE)
        }
    }

    count = 0
    multiple_v_genes <- FALSE
    multiple_j_genes <- FALSE
    for(i0 in align_seq2){
        count = count + 1
        tempo.log <- is.na(obs[ , i0]) | obs[ , i0] == ""
        if(sum(!tempo.log, na.rm = TRUE) >= align_clone_nb){
            # Only create fasta  files with at least <align_clone_nb> sequences (Minimun number of non-identical sequences per clonal group for tree plotting)
            # NB : align_clone_nb is defined in nextflow.config
            # Create the fasta files :
            if(any(tempo.log)){
                tempo.warn <- paste0("IMPORTED FILE:\n", path, "\nHAS ", sum(tempo.log, na.rm = TRUE), " AMONG ", nrow(obs), " EMPTY SEQUENCES (NA OR \"\") IN THE ", i0, " COLUMN IN LINES:\n", paste(which(tempo.log), collapse = "\n"))
                cat(paste0("\nWARNING IN ", script, ".R\n", tempo.warn, "\n\n"))
                fun_report(data = paste0("WARNING\n", tempo.warn), output = log, path = "./", overwrite = FALSE)
                warn <- paste0(ifelse(is.null(warn), tempo.warn, paste0(warn, "\n\n", tempo.warn)))
            }else{
                tempo.cat <- paste0("\nIMPORTED FILE:\n", path, "\nHAS NO EMPTY SEQUENCES (NA OR \"\") AMONG ", nrow(obs), " IN THE ", i0, " COLUMN\n")
                fun_report(data = tempo.cat, output = log, path = "./", overwrite = FALSE)
            }
            df <- obs[ ! tempo.log, ] # .tsv file with lines removed whem empty align_seq2
            if(seq_kind == "CLONE"){
                # Sequences in the same clone_assigned_seq.tsv file belong to the same clonal group and should have the same values in columns relative to clonal groups
                if(any(df$clone_id != df$clone_id[1])){
                    stop(paste0("\n\n================\n\nERROR IN ", script, ".R\nALL clone_id VALUES SHOULD BE THE SAME IN A clone_assigned_seq.tsv FILE, BUT THEY ARE NOT.\nHERE THEY ARE : ", paste0(df$clone_id, collapse = "\n"),"\n\n================\n\n"), call. = FALSE)
                }
            }
            if(seq_kind == "ALL"){
                multiple_v_genes <- TRUE
            }else if(any(df$v_gene != df$v_gene[1])) {
                tempo.warn <- paste0(
                    "Different genes for the V cassette were found in this clonal group\n",
                    "Here they are:\n", paste0(df$v_gene, collapse = "\n"),
                    "\nSet clone_strategy = \"first\" in nextflow.config to avoid this."
                )
                cat(paste0("\nWARNING IN ", script, ".R\n", tempo.warn, "\n\n"))
                fun_report(data = paste0("WARNING\n", tempo.warn), output = log, path = "./", overwrite = FALSE)
                warn <- paste0(ifelse(is.null(warn), tempo.warn, paste0(warn, "\n\n", tempo.warn)))
                multiple_v_genes <- TRUE
            }
            if(seq_kind == "ALL"){
                multiple_j_genes <- TRUE
            }else if(any(df$j_gene != df$j_gene[1])) {
                tempo.warn <- paste0(
                    "Different genes for the V cassette were found this clonal group\n",
                    "Here they are:\n", paste0(df$j_gene, collapse = "\n"),
                    "\nSet clone_strategy = \"first\" in nextflow.config to avoid this."
                )
                cat(paste0("\nWARNING IN ", script, ".R\n", tempo.warn, "\n\n"))
                fun_report(data = paste0("WARNING\n", tempo.warn), output = log, path = "./", overwrite = FALSE)
                warn <- paste0(ifelse(is.null(warn), tempo.warn, paste0(warn, "\n\n", tempo.warn)))
                multiple_j_genes <- TRUE
            }
            # if(any(df$v_gene != df$v_gene[1])){
            #     stop(paste0("\n\n================\n\nERROR IN ", script, ".R\nALL v_gene VALUES SHOULD BE THE SAME IN A clone_assigned_seq.tsv FILE, BUT THEY ARE NOT.\nHERE THEY ARE : ", paste0(df$v_gene, collapse = "\n"),"\n\n================\n\n"), call. = FALSE)
            # }
            # if(any(df$j_gene != df$j_gene[1])){
            #     stop(paste0("\n\n================\n\nERROR IN ", script, ".R\nALL j_gene VALUES SHOULD BE THE SAME IN A clone_assigned_seq.tsv FILE, BUT THEY ARE NOT.\nHERE THEY ARE : ", paste0(df$j_gene, collapse = "\n"),"\n\n================\n\n"), call. = FALSE)
            # }
            if(seq_kind == "CLONE"){
                if(any(df$junction_length != df$junction_length[1])){
                    stop(paste0("\n\n================\n\nERROR IN ", script, ".R\nALL junction_length VALUES SHOULD BE THE SAME IN A clone_assigned_seq.tsv FILE, BUT THEY ARE NOT.\nHERE THEY ARE : ", paste0(df$junction_length, collapse = "\n"),"\n\n================\n\n"), call. = FALSE)
                }
            }
            # End check of columns relative to clonal groups

            # Creation of the fasta file
            if(base::grepl(x = i0, pattern = "_aa$")){
                nuc_or_aa <- "aa"
            }else if(count > 2){
                stop(paste0("\n\n================\n\nERROR IN ", script, ".R\ncount CANNOT BE MORE THAN 2 : ", paste0(count, collapse = "\n"),"\n\n================\n\n"), call. = FALSE)
            }else{
                nuc_or_aa <- "nuc"
            }
            dir_name <- paste0("./for_alignment_", nuc_or_aa)
            dir.create(dir_name, showWarnings = FALSE, recursive = TRUE)

            # Extract all different values of v_genes and j_genes in the clonal group (different genes for a same cassette can be present in a clonal group if the clone_strategy = "set" in nextflow.config)
            # If several genes in a column, they are separated by a comma (",")
            v_genes_all <- unique(unlist(strsplit(df$v_gene, ",")))
            j_genes_all <- unique(unlist(strsplit(df$j_gene, ",")))
            # Replace the "," characters by "-" to have appropriate file names
            v_gene_clean <- paste(sort(v_genes_all), collapse = "-") 
            j_gene_clean <- paste(sort(j_genes_all), collapse = "-")

            if(seq_kind == "CLONE"){
                seq_name <- paste0(i0, "_clone_id_", df[1, "clone_id"], "_", v_gene_clean, "_", j_gene_clean) # Create the name of the file
            }else{
                seq_name <- i0
            }
            tempo.name <- paste0(seq_name, ".fasta")
            for(i1 in 1:nrow(df)) {
                tempo_seq <- gsub(x = df[i1, i0], pattern = "-", replacement = "")
                if(nuc_or_aa == "aa"){
                    tempo_seq <- gsub(x = tempo_seq, pattern = "\\*", replacement = "X")
                }
                tempo.cat <- paste0(">", df[i1, Name], "\n",tempo_seq, "\n")
                cat(tempo.cat, file = file.path(dir_name, tempo.name), append = TRUE)
            }

            # Germline addition to the fasta
            # this step is complicate in case of several germline sequences in the same clonal group but with clone_strategy = "first" in nextflow.config, it should be alright
            # Sort the dataframe by the column specified in Name
            if(seq_kind == "CLONE"){
                germ = Germline[count]
                df <- df[order(df[[Name]]), ]

                # Get the germline values and compute frequencies
                germ_seq <- df[[germ]] # sequences of germline_alignment_d_mask_no_gaps for instance
                germ_table <- table(germ_seq) # several or 1 sequence?
                max_freq <- max(germ_table)
                most_frequent_germs <- names(germ_table[germ_table == max_freq])

                # Choose the first matching germline by order in df (sorted by Name)
                selected_index <- which(df[[germ]] %in% most_frequent_germs)[1]
                selected_germline_seq <- gsub(x = df[selected_index, germ], pattern = "-", replacement = "")

                # Get cleaned v_gene and j_gene values for the chosen germline row
                germ_v_gene <- gsub(pattern = ",", replacement = "-", x = df[selected_index, "clonal_germline_v_gene"])
                germ_j_gene <- gsub(pattern = ",", replacement = "-", x = df[selected_index, "clonal_germline_j_gene"])

                # Write germline to the same fasta file as above
                germ_seq_name <- paste0(Germline[1], "_clone_id_", df[1, "clone_id"], "_", germ_v_gene, "_", germ_j_gene) # no choice to use Germline[1], otherwise goalign codonalign cannot align (because nuc and aa must have the same name)
                tempo.cat <- paste0(">", germ_seq_name, "\n", selected_germline_seq, "\n")
                cat(tempo.cat, file = file.path(dir_name, tempo.name), append = TRUE)

                # If germlines differ in the group, issue a warning
                if (length(unique(germ_seq)) > 1) {
                    tempo.warn <- paste0(
                        "Multiple germline sequences found for clone_id ", df[1, "clone_id"], ".\n",
                        "Selected germline from sequence: ", df[selected_index, Name], "\n",
                        "Associated clonal_germline_v_gene: ", df[selected_index, "clonal_germline_v_gene"], "\n",
                        "Associated clonal_germline_j_gene: ", df[selected_index, "clonal_germline_j_gene"], "\n",
                        "Associated v_gene: ", df[selected_index, "v_gene"], "\n",
                        "Associated j_gene: ", df[selected_index, "j_gene"], "\n",
                        "Used germline sequence: ", selected_germline_seq
                    )
                    cat(paste0("\nWARNING IN ", script, ".R\n", tempo.warn, "\n\n"))
                    fun_report(data = paste0("WARNING\n", tempo.warn), output = log, path = "./", overwrite = FALSE)
                    warn <- paste0(ifelse(is.null(warn), tempo.warn, paste0(warn, "\n\n", tempo.warn)))
                }
            }
        }else{
            stop(paste0("\n\n================\n\nERROR IN ", script, ".R\nNO FASTA AND GFF FILES CREATED BECAUSE THE IMPORTED FILE:\n", path, "\nHAS LESS THAN ", align_clone_nb, " SEQUENCES (TOO MUCH NA OR \"\") IN THE ", i0, " COLUMN\n\n================\n\n"))
        }
    }
}

if(seq_kind == "CLONE"){
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




