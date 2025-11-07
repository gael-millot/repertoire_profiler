#!/usr/bin/env Rscript

#########################################################################
##                                                                     ##
##     Gff.R                                                           ##
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
script <- "Gff"
#cute <- "https://gitlab.pasteur.fr/gmillot/cute_little_R_functions/-/raw/v11.4.0/cute_little_R_functions.R" # single character string indicating the path of the cute_little_R_functions.R file required for this script. Example: cute <- "https://gitlab.pasteur.fr/gmillot/cute_little_R_functions/-/raw/v11.4.0/cute_little_R_functions.R"
#log <- "xlsx2fasta.log" # single character string indicating the name of the log file. Example: log <- "xlsx2fasta.log"


################################ End Initialization



################################ EXAMPLES FOR TEST AND EXPLANATION OF ARGUMENTS

# script <- "Tsv2fastaGff"

### Arguments : 

# fasta_path <- "C:/Users/gmillot/Documents/Git_projects/repertoire_profiler/example_of_results/repertoire_profiler_1762523781/alignments/nuc/sequence_alignment_aligned_nuc.fasta"
# tsv_path <- "C:/Users/gmillot/Documents/Git_projects/repertoire_profiler/example_of_results/repertoire_profiler_1762523781/tsv/productive_seq.tsv"
# Name <- "sequence_id"                # name of the column containing the sequence ids
# align_seq <- "sequence_alignment"        # name of the columns containing the sequences to put in the fasta file (can be a single string or several strings seperated by "," if several columns are needed. the fastas will then be created in different folders)
# clone_germline_kind <- "dmask"
# align_clone_nb <- 3                    # Minimum number of rows in the tsv file. The program expects this to be respected, otherwise raises an error.
# cute <- "C:/Users/gmillot/Documents/Git_projects/repertoire_profiler/bin/cute_little_R_functions_v12.8.R"
# tag <- "CLONE"
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
        "fasta_path", 
        "tsv_path", 
        "Name", # sequence_id column
        "align_seq",
        "clone_germline_kind", 
        "align_clone_nb",
        "cute", 
        "tag", # ALL or CLONE
        "log"
    ) # objects names exactly in the same order as in the bash code and recovered in args. Here only one, because only the fasta_path of the config file to indicate after the xlsx2fasta.R script execution
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
    "fasta_path", 
    "tsv_path", 
    "Name", 
    "align_seq", 
    "clone_germline_kind", 
    "align_clone_nb", 
    "cute", 
    "tag", 
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
# Test if one fasta_path exists (url or local)
# ARGUMENTS
# fasta_path: single character string of the fasta_path to test
# script: single character string of the current script file
# RETURN
# An error message if the fasta_path does not exists, nothing otherwise
# EXAMPLES
# fun_source_test(path = "caca", script = "test")
# DEBUGGING
# fasta_path = "https://zenodo.org/records/10814482/files/ig_sequences.xlsx" ; script = "test"
    name <- deparse(substitute(fasta_path))
    if(length(fasta_path) != 1){
        stop(paste0("\n\n============\n\nERROR IN ", script, ".R\n", name, " PARAMETER MUST BE LENGTH 1: ", paste(fasta_path, collapse = " "), "\n\n============\n\n"), call. = FALSE)
    }else if(grepl(x = fasta_path, pattern = "^http")){
        tempo.name <- paste0("tmp_xlsx2fasta.R_", as.integer(Sys.time()))
        if(file.exists(tempo.name)){
            stop(paste0("\n\n============\n\nERROR IN ", script, ".R\nTHE TEMPORARY FILE USED BY THE ", script, " SCRIPT ALREADY EXISTS:\n", file.path(tempo.name), "\n\n. PLEASE, RERUN.\n\n============\n\n"), call. = FALSE)
        }else{
            tempo.try <- try(suppressWarnings(suppressMessages(download.file(fasta_path, tempo.name, method="auto", quiet=TRUE))), silent = TRUE)
            if(file.exists(tempo.name)){
                file.remove(tempo.name)
            }
            if(any(grepl(x = tempo.try, pattern = "^[Ee]rror"))){
                stop(paste0("\n\n============\n\nERROR IN ", script, ".R\nHTTP INDICATED IN THE ", name, " PARAMETER DOES NOT EXISTS: ", fasta_path, "\n\n============\n\n"), call. = FALSE)
            }
        }
    }else if( ! grepl(x = fasta_path, pattern = "^http")){
        if( ! file.exists(fasta_path)){
            stop(paste0("\n\n============\n\nERROR IN ", script, ".R\nFILE INDICATED IN THE ", name, " PARAMETER DOES NOT EXISTS: ", fasta_path, "\n\n============\n\n"), call. = FALSE)
        }
    }else{
        tempo.cat <- paste0("\n\n================\n\nINTERNAL CODE ERROR 3 IN ", script, ".R: CODE HAS TO BE MODIFIED\n\n============\n\n")
        stop(tempo.cat, call. = FALSE)
    }
}

map_ungapped_to_gapped <- function(seq_aligned, ungapped_pos) {
# AIM
# shift vdjc or fwr/cdr position provided in the .tsv file, according to the hyphens (--) inserted in the aligned sequences, to have the correct positions of the features in the aligned sequneces in the html file.
# ARGUMENTS
# seq_aligned: aligned sequence (with hyphens)
# ungapped_pos: initial positions in the unaligned sequence (without hyphens)
# RETURN
# shifted position
# EXAMPLES
# map_ungapped_to_gapped("atg--t---g", c(3,4,5)) # initial sequence is atgtg, input coords are third g, fourth t and fifth g
# DEBUGGING
# seq_aligned = "atg--t---g" ; ungapped_pos = c(3,4,5)
  # Remove gaps but keep index positions of non-gaps
  aligned_chars <- strsplit(seq_aligned, "")[[1]]
  non_gap_indices <- which(aligned_chars != "-")
  # Return the gapped index corresponding to the ungapped position
  return(non_gap_indices[ungapped_pos])
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

# end reserved words
# argument primary checking
arg.check <- NULL #
text.check <- NULL #
checked.arg.names <- NULL # for function debbuging: used by r_debugging_tools
ee <- expression(arg.check <- c(arg.check, tempo$problem) , text.check <- c(text.check, tempo$text) , checked.arg.names <- c(checked.arg.names, tempo$object.name))
tempo <- fun_check(data = fasta_path, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = Name, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = align_seq, options = c("query", "igblast_full", "trimmed", "fwr1", "fwr2", "fwr3", "fwr4", "cdr1", "cdr2", "cdr3", "junction", "sequence_alignment", "v_sequence_alignment", "d_sequence_alignment", "j_sequence_alignment", "c_sequence_alignment", "germline_alignment", "v_germline_alignment", "d_germline_alignment", "j_germline_alignment", "c_germline_alignment")) ; eval(ee)
tempo <- fun_check(data = clone_germline_kind, options = c("full","dmask","vonly")) ; eval(ee)
tempo <- fun_check(data = align_clone_nb, class = "vector", typeof = "integer", length = 1, double.as.integer.allowed = TRUE) ; eval(ee)
tempo <- fun_check(data = tag, options = c("ALL", "CLONE")) ; eval(ee)
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
    "fasta_path", 
    "tsv_path", 
    "Name", 
    "align_seq", 
    "clone_germline_kind", 
    "align_clone_nb",
    "cute", 
    "tag", 
    "log"
)
tempo_log <- sapply(lapply(tempo.arg, FUN = get, env = sys.nframe(), inherit = FALSE), FUN = is.null)
if(any(tempo_log) == TRUE){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN ", script, ".R:\n", ifelse(sum(tempo_log, na.rm = TRUE) > 1, "THESE ARGUMENTS\n", "THIS ARGUMENT\n"), paste0(tempo.arg[tempo_log], collapse = "\n"),"\nCANNOT BE NULL")
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end management of NULL arguments
# management of ""
tempo.arg <-c(
    "fasta_path", 
    "tsv_path", 
    "Name", 
    "align_seq", 
    "clone_germline_kind", 
    "cute", 
    "tag", 
    "log"
)
tempo_log <- sapply(lapply(tempo.arg, FUN = get, env = sys.nframe(), inherit = FALSE), FUN = function(x){any(x == "")})
if(any(tempo_log) == TRUE){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN ", script, ".R:\n", ifelse(sum(tempo_log, na.rm = TRUE) > 1, "THESE ARGUMENTS\n", "THIS ARGUMENT\n"), paste0(tempo.arg[tempo_log], collapse = "\n"),"\nCANNOT BE \"\"")
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
fun_source_test(path = fasta_path, script = script)
fun_source_test(path = tsv_path, script = script)
if(grepl(x = tsv_path, pattern = ".*productive_seq.tsv$") & tag == "CLONE"){
    stop(paste0("\n\n============\n\nERROR IN ", script, ".R\n\nTag CANNOT BE \"CLONE\" IF tsv_path POINTS TO productive_seq.tsv FILE\n\n============\n\n"), call. = FALSE)
}
if(grepl(x = tsv_path, pattern = ".*clone_assigned_seq.tsv$") & tag == "ALL"){
    stop(paste0("\n\n============\n\nERROR IN ", script, ".R\n\nTag CANNOT BE \"CLONE\" IF tsv_path POINTS TO clone_assigned_seq.tsv FILE\n\n============\n\n"), call. = FALSE)
}

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
fun_report(data = paste0("\n\n#########\n\n", fasta_path, "\n\n#########\n\n"), output = log, path = "./", overwrite = FALSE)

################ End ignition


################ Graphical parameter initialization


################ End graphical parameter initialization


################ Data import

df <- read.table(tsv_path, header = TRUE, sep = "\t")
fasta <- readLines(fasta_path)

################ end Data import


############ check


if( ! Name %in% names(df)){
    stop(paste0("\n\n============\n\nERROR IN ", script, ".R\n\nTHE Name PARAMETER MUST BE A COLUMN NAME OF THE IMPORTED FILE:\n", tsv_path, "\n\nHERE IT IS Name:\n", Name, "\n\nCOLUMN NAMES:\n", paste(names(df), collapse = "\n"), "\n\n============\n\n"), call. = FALSE)
}

if(any(duplicated(df[, Name]))){
    stop(paste0("\n\n============\n\nERROR IN ", script, ".R\n\nDUPLICATED VALUE NOT AUTHORIZED IN THE COLUMN OF THE Name PARAMETER\n\nDUPLICATED VALUES ARE:\n", df[duplicated(df[, Name]), Name], "\n\nIN POSITIONS:\n", paste(which(df[ , Name] %in% df[duplicated(df[, Name]), Name]), collapse = "\n"), "\n\n============\n\n"), call. = FALSE)
}


tempo_log <- is.na(df[ , Name]) | df[ , Name] == ""
if(any(tempo_log)){
    stop(paste0("\n\n============\n\nERROR IN ", script, ".R\nIMPORTED FILE:\n", tsv_path, "\nHAS AN EMPTY CELL IN THE ", Name, " COLUMN IN LINES:\n", paste(which(tempo_log), collapse = "\n"), "\n\n============\n\n"), call. = FALSE)
}

for(i0 in names(df)){ # NA in xlsx file become "NA". Thus, has to be replaced by NA
    tempo_log <- df[ , i0] == "NA" & ! is.na(df[ , i0])
    if(any(tempo_log, na.rm = TRUE)){
        df[tempo_log, i0] <- NA
    }
}


############ end check


############ main


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

Germline <- NULL
if(clone_germline_kind == "dmask" & tag == "CLONE" ){
    Germline <- c("germline_alignment_d_mask_no_gaps", "germline_alignment_d_mask_aa_no_gaps")
}else if(clone_germline_kind == "vonly" & tag == "CLONE" ){
    Germline <- c("germline_alignment_v_region_no_gaps", "germline_alignment_v_region_aa_no_gaps")
}else if(clone_germline_kind == "full" & tag == "CLONE" ){
    Germline <- c("germline_alignment_full_no_gaps", "germline_alignment_full_aa_no_gaps")
}


if( ! all(align_seq2 %in% names(df))){
    stop(paste0("\n\n============\n\nERROR IN ", script, ".R\n\nTHE align_seq2 PARAMETER MUST BE COLUMN NAMES OF THE IMPORTED FILE:\n", tsv_path, "\n\nHERE IT IS align_seq2:\n", paste(align_seq2, collapse = "\n"), "\n\nCOLUMN NAMES:\n", paste(names(df), collapse = "\n"), "\n\n============\n\n"), call. = FALSE)
}

if(base::grepl(x = fasta_path, pattern = "*aa.fasta$")){
    nuc_or_aa <- "aa"
    if(is.null(Germline)){
        align_seq3 <- align_seq2[2]
    }else{
        align_seq3 <- Germline[2]
    }
}else if(base::grepl(x = fasta_path, pattern = "*nuc.fasta$")){
    nuc_or_aa <- "nuc"
    if(is.null(Germline)){
        align_seq3 <- align_seq2[1]
    }else{
        align_seq3 <- Germline[1]
    }
}

tempo_log <- is.na(df[ , align_seq3]) | df[ , align_seq3] == ""
if(sum(!tempo_log, na.rm = TRUE) >= align_clone_nb){
    # Create the gff file
    # https://shazam.readthedocs.io/en/stable/topics/setRegionBoundaries/
    # define the column names of the coordinates in the tsv file
    vdjc_features <- c("v", "d", "j", "c")
    fwr_cdr_features <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4") 
    vdjc_features_colors <- c("red", "green", "blue", "yellow")
    fwr_cdr_features_colors <- c("yellow", "pink", "yellow", "pink", "yellow", "pink", "yellow")
    if(tag == "CLONE"){ # means that the column taken are coordinates of the germline sequence only. All the columns in the clone_assign_seq.tsv start by "clonal_germline_"
        vdjc_column_start <- paste0("clonal_germline_", vdjc_features, "_start")
        vdjc_column_end <- paste0("clonal_germline_", vdjc_features, "_end")
        fwr_cdr_column_start <- paste0("clonal_germline_", fwr_cdr_features, "_start")
        fwr_cdr_column_end <- paste0("clonal_germline_", fwr_cdr_features, "_end")
    }else{ # tag == "ALL"
        if(grepl(x = align_seq, pattern = ".*_alignment")){ # for align_seq %in% sequence_alignment|v_sequence_alignment|d_sequence_alignment|j_sequence_alignment|c_sequence_alignment|germline_alignment|v_germline_alignment|d_germline_alignment|j_germline_alignment|c_germline_alignment. These coordinates are for both aligned germline and sequence 
            vdjc_column_start <- paste0(vdjc_features, "_alignment_start")
            vdjc_column_end <- paste0(vdjc_features, "_alignment_end")
            fwr_cdr_column_start <- NULL
            fwr_cdr_column_end <- NULL
        }else if(align_seq %in% c("query", "igblast_full", "trimmed")){
            vdjc_column_start <- paste0(vdjc_features, "_sequence_start")
            vdjc_column_end <- paste0(vdjc_features, "_sequence_end")
            fwr_cdr_column_start <- paste0(fwr_cdr_features, "_start")
            fwr_cdr_column_end <- paste0(fwr_cdr_features, "_end")
        }else{ # for align_seq %in% fwr1|fwr2|fwr3|fwr4|cdr1|cdr2|cdr3|junction
            vdjc_column_start <- NULL
            vdjc_column_end <- NULL
            fwr_cdr_column_start <- NULL
            fwr_cdr_column_end <- NULL
        }
    }
    # end define the column names of the coordinates in the tsv file
    # Check that all expected columns exist
    not_exist <- NULL
    for(i2 in c("vdjc_column_start", "vdjc_column_end", "fwr_cdr_column_start", "fwr_cdr_column_end")){
        if( ! is.null(get(i2))){
            tempo_log <- ! get(i2) %in% names(df)
            if(any(tempo_log, na.rm = TRUE)){
                not_exist <- c(not_exist, get(i2)[tempo_log])
            }
        }
    }
    if (length(not_exist) > 0) {
        stop(paste0("\n\n================\n\nERROR IN ", script, ".R\nONE OR MORE COORDINATE COLUMNS MISSING FROM THE IMPORTED FILE:\n", tsv_path, "\nHERE IS THE MISSING COLUMN : ", paste0(not_exist, collapse = "\n"), "\n\n================\n\n"), call. = FALSE)
    }
    for(i2 in c("vdjc_column", "fwr_cdr_column")){
        if((is.null(paste0(i2, "_start")) & ! is.null(paste0(i2, "_end"))) | ( ! is.null(paste0(i2, "_start")) & is.null(paste0(i2, "_end")))){
            stop(paste0("\n\n================\n\nERROR IN ", script, ".R\nNEITHER OR BOTH ", i2, " _start AND _end OBJECTS MUST BE NULL", paste0(not_exist, collapse = "\n"), "\n\n================\n\n"), call. = FALSE)
        }
    }
    # end Check that all expected columns exist
    # get the names of sequences
    seq_names <- fasta[grepl(x = fasta, pattern = "^>")]
    seq_aligned <- fasta[ ! grepl(x = fasta, pattern = "^>")]
    seq_names_for_tsv <- sub(x = seq_names, pattern = "^>", replacement = "")
    selected_index <- match(seq_names_for_tsv, df$sequence_id) # line number in df
    if(tag == "CLONE"){ # coordinates only for germline seq
        tempo_log <- grepl(x = seq_names_for_tsv, pattern = paste0("^>", align_seq3))
        seq_names <- seq_names[tempo_log]
        seq_aligned <- seq_aligned[tempo_log]
        seq_names_for_tsv <- seq_names_for_tsv[tempo_log]
        clone_id <- as.integer(sub(pattern = ".*clone_id_([0-9]+)_.*", replacement = "\\1", x = seq_names_for_tsv))
        selected_index <- which(df$clone_id == clone_id)[1] # first line taken for coordinates
    }
    # end get the names of sequences
    # get coordinates
    for(i2 in c("vdjc", "fwr_cdr")){
        gff_rows <- list()
        for(i3 in 1:length(selected_index)){
            for(i4 in 1:length(get(paste0(i2, "_features")))){
                if( ! is.null(get(paste0(i2, "_column_start")))){
                    # nuc coordinates
                    start_coord <- df[selected_index[i3], get(paste0(i2, "_column_start"))[i4]]
                    end_coord <- df[selected_index[i3], get(paste0(i2, "_column_end"))[i4]]
                    if(nuc_or_aa == "nuc"){
                        if( ! is.na(start_coord)){
                            # coord subtraction for the trimmed sequence, because coordinates are those of the query sequence
                            if(tag == "ALL" & align_seq == "trimmed"){
                                rm_seq <- df[selected_index[i3], "removed_sequence"]
                                if(is.na(rm_seq)){
                                    length_to_rm <- 0
                                }else{
                                    length_to_rm <- length(rm_seq)
                                }
                                start_coord <- start_coord - length_to_rm
                                if (start_coord < 1) {
                                    stop(paste0("\n\n================\n\nERROR IN ", script, ".R\nAFTER CORRECTION OF THE START COORDINATE ACCORDING TO THE LENGTH OF THE REMOVED PART IN THE TRIMMED SEQ, IT IS LESS THAN 1:\n", start_coord, "\nSEQ NAME:\n", seq_names_for_tsv[selected_index[i3]], "\nALIGNED SEQ:\n", seq_aligned[selected_index[i3]], "REMOVED SEQ:\n", rm_seq, "\n\n================\n\n"), call. = FALSE)
                                }
                            }
                            # end coord subtraction for the trimmed sequence, because coordinates are those of the query sequence
                                # nuc column
                            # shifed coordinates due to hyphens in the aligned seq
                            start_coord <- map_ungapped_to_gapped(seq_aligned = seq_aligned[selected_index[i3]], ungapped_pos = start_coord)
                            # end shifed coordinates due to hyphens in the aligned seq
                        }
                        if( ! is.na(end_coord)){
                            # coord subtraction for the trimmed sequence, because coordinates are those of the query sequence
                            if(tag == "ALL" & align_seq == "trimmed"){
                                rm_seq <- df[selected_index[i3], "removed_sequence"]
                                if(is.na(rm_seq)){
                                    length_to_rm <- 0
                                }else{
                                    length_to_rm <- length(rm_seq)
                                }
                                end_coord <- end_coord - length_to_rm
                                if (end_coord < 1) {
                                    stop(paste0("\n\n================\n\nERROR IN ", script, ".R\nAFTER CORRECTION OF THE END COORDINATE ACCORDING TO THE LENGTH OF THE REMOVED PART IN THE TRIMMED SEQ, IT IS LESS THAN 1 (AND START SEQ IS NOT !?):\n", end_coord, "\nSEQ NAME:\n", seq_names_for_tsv[selected_index[i3]], "\nALIGNED SEQ:\n", seq_aligned[selected_index[i3]], "REMOVED SEQ:\n", rm_seq, "\n\n================\n\n"), call. = FALSE)
                                }
                            }
                            # end coord subtraction for the trimmed sequence, because coordinates are those of the query sequence
                                # nuc column
                            # shifed coordinates due to hyphens in the aligned seq
                            end_coord <- map_ungapped_to_gapped(seq_aligned = seq_aligned[selected_index[i3]], ungapped_pos = end_coord)
                            # end shifed coordinates due to hyphens in the aligned seq
                        }
                        # end nuc coordinates
                    # aa coordinates
                    }else if(nuc_or_aa == "aa"){
                        for(i7 in c("start", "end")){
                            tempo_coord_name <- paste0(i7, "_coord") # start_coord or stop_coord
                            if(is.na(get(tempo_coord_name))){
                                assign(tempo_coord_name, NA)
                            }else{
                                if(i7 == "start"){
                                    coord <- as.integer(get(tempo_coord_name)) - 1
                                }else{
                                    coord <- as.integer(get(tempo_coord_name))
                                }
                                if(as.integer(coord) %% 3 == 0){
                                    if(i7 == "start"){
                                        assign(tempo_coord_name, as.integer(coord / 3 + 1))
                                    }else{
                                        assign(tempo_coord_name, as.integer(coord / 3))
                                    }
                                }else{
                                    if(i7 == "start"){
                                        assign(tempo_coord_name, as.integer(trunc(coord / 3) + 1))
                                    }else{
                                        assign(tempo_coord_name, as.integer(trunc(coord / 3)))
                                    }
                                    tempo.warn <- paste0("APPROXIMATE AA COORDINATES FOR ", df[selected_index[i3], Name], ", FOR THE ", get(paste0(i2, "_column_", i7))[i4], " COLUMN, SINCE NUC COORDINATES ARE NOT MULTIPLE OF 3.")
                                    cat(paste0("\nWARNING IN ", script, ".R\n", tempo.warn, "\n\n"))
                                    fun_report(data = paste0("WARNING\n", tempo.warn), output = log, path = "./", overwrite = FALSE)
                                    warn <- paste0(ifelse(is.null(warn), tempo.warn, paste0(warn, "\n\n", tempo.warn)))
                                }
                            }
                        }
                    }
                    # end aa coordinates
                    gff_rows[[length(gff_rows) + 1]] <- c(
                        if(tag == "ALL"){
                            df[selected_index[i3], Name]
                        }else{
                            seq_names_for_tsv
                        },
                        ".",
                        "gene",
                        start_coord, 
                        end_coord, 
                        ".",
                        ".",
                        ".",
                        paste0("Name=", get(paste0(i2, "_features"))[i4], ";Color=", get(paste0(i2, "_features_colors"))[i4])
                    )
                }
            }
        }
        gff_table <- do.call(rbind, gff_rows)
        if(length(gff_table) > 0){
            gff_lines <- apply(gff_table, 1, function(x) paste(x, collapse="\t"))
        }else{
            gff_lines <- character()
        }
        gff_lines <- c("##gff-version 3", gff_lines)
        output_gff <- paste0(i2, "_", sub(x = basename(fasta_path), pattern = ".fasta$", replacement = ""), ".gff")
        writeLines(gff_lines, con = output_gff)
    }
    # end Create the gff file
} else {
    stop(paste0("\n\n================\n\nERROR IN ", script, ".R\nNO FASTA AND GFF FILES CREATED BECAUSE THE IMPORTED FILE:\n", tsv_path, "\nHAS LESS THAN ", align_clone_nb, " SEQUENCES (TOO MUCH NA OR \"\") IN THE ", align_seq3, " COLUMN\n\n================\n\n"))
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




