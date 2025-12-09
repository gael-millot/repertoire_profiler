#!/usr/bin/env Rscript

#########################################################################
##                                                                     ##
##     Gff_imgt.R                                                      ##
##                                                                     ##
##     Gael A. Millot                                                  ##
##     Chloe Taurel                                                    ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Computational Biology Department                                ##
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
script <- "Gff_imgt"

################################ End Initialization



################################ EXAMPLES FOR TEST AND EXPLANATION OF ARGUMENTS

# script <- "Tsv2fastaGff"

### Arguments : 

# fasta_path <- "C:/Users/gmillot/Documents/Git_projects/repertoire_profiler/results/repertoire_profiler_1765195115/alignments/nuc/imgt/sequence_alignment_with_gaps_imgt_nuc.fasta"
# tsv_path <- "C:/Users/gmillot/Documents/Git_projects/repertoire_profiler/results/repertoire_profiler_1765195115/tsv/productive_seq.tsv"
# Name <- "sequence_id"                # name of the column containing the sequence ids
# align_seq <- "query"                # kind of seq
# cute <- "C:/Users/gmillot/Documents/Git_projects/repertoire_profiler/bin/cute_little_R_functions_v12.8.R"
# log <- "Gff_imgt.log"


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
        "fasta_aa_path", 
        "tsv_path", 
        "Name", # sequence_id column
        "align_seq", 
        "cute", 
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
    "fasta_aa_path", 
    "tsv_path", 
    "Name", 
    "align_seq", 
    "cute", 
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


map_gapped_to_ungapped <- function(seq_aligned, gapped_pos, script) {
# AIM
# shift vdjc or fwr/cdr position provided in the .tsv file, according to the removal of hyphens (--) and dots (...) inserted in the aligned sequences, to have the correct positions of the features in the aligned sequences in jalview.
# ARGUMENTS
# seq_aligned: aligned sequence (with hyphens)
# ungapped_pos: initial positions in the aligned sequence (hyphens and dots)
# RETURN
# shifted position
# EXAMPLES
# map_gapped_to_ungapped("atg-.t---g", c(3, 6, 7), "test") # initial sequence is atgtg, input coords are third g, fourth t and a hyphen
# DEBUGGING
# seq_aligned = "atg--t---g" ; gapped_pos = c(3, 6, 7) ; script = "test"
    # Split sequence into characters
    seq_chars <- unlist(strsplit(seq_aligned, ""))
    # Check if any coordinate points to '.' or '-'
    if(any(seq_chars[gapped_pos] %in% c(".", "-"))) {
        bad <- gapped_pos[seq_chars[gapped_pos] %in% c(".", "-")]
        # tempo.cat <- paste0("ERROR IN ", script, ".R\nThe following coordinates point to '.' or '-' in the original sequence:\n", paste0(bad, collapse = "\n"))
        # stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
    }
    # Create a map from original positions to cleaned sequence positions
    valid_idx <- which(!seq_chars %in% c(".", "-"))
    mapping <- seq_along(valid_idx)
    names(mapping) <- valid_idx
    # Convert original coordinates to coordinates in cleaned sequence
    new_coords <- as.integer(mapping[as.character(gapped_pos)])
    # Return the gapped index corresponding to the ungapped position
    return(new_coords)
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
tempo <- fun_check(data = fasta_aa_path, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = Name, class = "vector", typeof = "character", length = 1) ; eval(ee)
# cute already tested above
tempo <- fun_check(data = align_seq, options = c("query", "igblast_full", "trimmed", "fwr1", "fwr2", "fwr3", "fwr4", "cdr1", "cdr2", "cdr3", "junction", "sequence_alignment", "v_sequence_alignment", "d_sequence_alignment", "j_sequence_alignment", "c_sequence_alignment", "germline_alignment", "v_germline_alignment", "d_germline_alignment", "j_germline_alignment", "c_germline_alignment")) ; eval(ee)
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
    "fasta_aa_path", 
    "tsv_path", 
    "Name", 
    "align_seq", 
    "cute", 
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
    "fasta_aa_path", 
    "tsv_path", 
    "Name", 
    "align_seq", 
    "cute", 
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
fun_source_test(path = fasta_aa_path, script = script)
fun_source_test(path = tsv_path, script = script)

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
fasta_aa <- readLines(fasta_aa_path)

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

# test that all non trimmed sequences start at 1 for fwr1_start, or NA
tempo_log <- df$fwr1_start != 1 & df$is_sequence_trimmed == FALSE
if(any(tempo_log, na.rm = TRUE)){ # NA not considered
    stop(paste0("\n\n============\n\nERROR IN ", script, ".R\nIMPORTED FILE:\n", tsv_path, "\nHAS fwr1_start COLUMN OF VALUE NOT 1, WHILE is_sequence_trimmed COLUMN IS TRUE, IN LINES:\n", paste(which(tempo_log), collapse = "\n"), "\n\n============\n\n"), call. = FALSE)
}
tempo_log <- df$fwr1_start == 1 & ! is.na(df$removed_sequence)
if(any(tempo_log, na.rm = TRUE)){ # NA not considered
    stop(paste0("\n\n============\n\nERROR IN ", script, ".R\nIMPORTED FILE:\n", tsv_path, "\nHAS fwr1_start COLUMN OF VALUE 1, WHILE removed_sequence COLUMN IS NOT NA, IN LINES:\n", paste(which(tempo_log), collapse = "\n"), "\n\n============\n\n"), call. = FALSE)
}
tempo_log <- df$fwr1_start == 1 & ! is.na(df$removed_sequence)
if(any(tempo_log, na.rm = TRUE)){ # NA not considered
    stop(paste0("\n\n============\n\nERROR IN ", script, ".R\nIMPORTED FILE:\n", tsv_path, "\nHAS fwr1_start COLUMN OF VALUE 1, WHILE removed_sequence COLUMN IS NOT NA, IN LINES:\n", paste(which(tempo_log), collapse = "\n"), "\n\n============\n\n"), call. = FALSE)
}
# end test that all non trimmed sequences start at 1 for fwr1_start, or NA


############ end check


############ main


if( ! "sequence_alignment_with_gaps" %in% names(df)){
    stop(paste0("\n\n============\n\nERROR IN ", script, ".R\n\nsequence_alignment_with_gaps MUST BE COLUMN A COLUMN NAME OF THE IMPORTED FILE:\n", path, "\n\nCOLUMN NAMES:\n", paste(names(df), collapse = "\n"), "\n\n============\n\n"), call. = FALSE)
}


# tempo_log <- is.na(df$sequence_alignment_with_gaps) | dfdf$sequence_alignment_with_gaps == ""
# Create the gff file
# https://shazam.readthedocs.io/en/stable/topics/setRegionBoundaries/
# define the column names of the coordinates in the tsv file
vdjc_features <- c("v", "d", "j", "c")
fwr_cdr_features <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4") 
    vdjc_features_colors <- c("#ffb6C1", "#90ee90", "#ffc878", "#0395c5")
    fwr_cdr_features_colors <- c("#ffb6C1", "#90ee90", "#ffc878", "#0395c5", "#f17cf1", "#02e6e6", "#d1d104")

vdjc_column_start <- paste0(vdjc_features, "_alignment_with_gaps_start") 
vdjc_column_end <- paste0(vdjc_features, "_alignment_with_gaps_end") 
fwr_cdr_column_start <- paste0(fwr_cdr_features, "_alignment_with_gaps_start") 
fwr_cdr_column_end <- paste0(fwr_cdr_features, "_alignment_with_gaps_end") 

vdjc_column_start_jalv <- paste0(vdjc_features, "_alignment_start") 
vdjc_column_end_jalv <- paste0(vdjc_features, "_alignment_end") 
fwr_cdr_column_start_jalv <- paste0(fwr_cdr_features, "_alignment_start") 
fwr_cdr_column_end_jalv <- paste0(fwr_cdr_features, "_alignment_end") 

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
selected_index <- match(seq_names_for_tsv, df$sequence_id) # line number in df. Warning: only use this in df
# end get the names of sequences
# get the aa sequences
seq_names_aa <- fasta_aa[grepl(x = fasta_aa, pattern = "^>")]
seq_aligned_aa <- fasta_aa[ ! grepl(x = fasta_aa, pattern = "^>")]
if( ! all(seq_names == seq_names_aa)){
    stop(paste0("\n\n================\n\nERROR IN ", script, ".R\nNAMES OF SEQUENCES IN THE NUC AND AA ALIGNMENT FASTA FILES SHOULD BE IDENTICAL AND IN THE SAME ORDER.\n\nHERE NUC ARE:", paste0(seq_names, collapse = "\n"), "\n\nAA ARE:", paste0(seq_names_aa, collapse = "\n"), "\n\n================\n\n"), call. = FALSE)
}
# end get the aa sequences
# get coordinates
for(i2 in c("vdjc", "fwr_cdr")){
    gff_rows <- list()
    gff_aa_rows <- list()
    gff_rows_jalv <- list()
    gff_aa_rows_jalv <- list()
    for(i3 in 1:length(selected_index)){
        for(i4 in 1:length(get(paste0(i2, "_features")))){
            if( ! is.null(get(paste0(i2, "_column_start")))){
                # nuc coordinates
                start_coord <- df[selected_index[i3], get(paste0(i2, "_column_start"))[i4]]
                end_coord <- df[selected_index[i3], get(paste0(i2, "_column_end"))[i4]]
                start_coord_jalv <- df[selected_index[i3], get(paste0(i2, "_column_start_jalv"))[i4]]
                end_coord_jalv <- df[selected_index[i3], get(paste0(i2, "_column_end_jalv"))[i4]]
                gff_rows[[length(gff_rows) + 1]] <- c(
                    df[selected_index[i3], Name],
                    ".",
                    "gene",
                    start_coord, 
                    end_coord, 
                    ".",
                    ".",
                    ".",
                    paste0("Name=", get(paste0(i2, "_features"))[i4], ";Color=", get(paste0(i2, "_features_colors"))[i4])
                )
                gff_rows_jalv[[length(gff_rows_jalv) + 1]] <- c(
                    df[selected_index[i3], Name],
                    ".",
                    get(paste0(i2, "_features"))[i4],
                    start_coord_jalv, 
                    end_coord_jalv, 
                    ".",
                    ".",
                    ".",
                    "."
                )
                # end nuc coordinates
                # aa coordinates
                if(is.na(start_coord)){
                    start_coord_aa <- NA
                }else{
                    start_coord_aa <- as.integer(trunc(start_coord / 3) + 1)
                    if((as.integer(start_coord) - 1) %% 3 != 0){ # == 0 means start_coord is at starting nuc codon 1, 4, 7, 10, etc., != 0 means start_coord is not at starting nuc codon -> shift to the next codon
                        start_coord_aa <- start_coord_aa + 1
                    }
                }
                if(is.na(end_coord)){
                    end_coord_aa <- NA
                }else{
                    end_coord_aa <- as.integer(trunc(end_coord / 3)) # if(as.integer(end_coord) %% 3 != 0) not required here because 1) == 0 means end_coord is at ending nuc codon 3, 6, 9, etc., 2) != 0 means end_coord is not at ending nuc codon -> translation put a hyphen for this codon. But 8 (middle of 3rd codon) -> 8/3 -> 2nd aa: if the 3rd codon is hyphen, then still ok to end at aa 2. So nothing to change
                }
                gff_aa_rows[[length(gff_aa_rows) + 1]] <- c(
                    df[selected_index[i3], Name],
                    ".",
                    "gene",
                    start_coord_aa, 
                    end_coord_aa, 
                    ".",
                    ".",
                    ".",
                    paste0("Name=", get(paste0(i2, "_features"))[i4], ";Color=", get(paste0(i2, "_features_colors"))[i4])
                )
                if(is.na(start_coord_aa)){
                    start_coord_jalv_aa <- NA
                }else{
                    start_coord_jalv_aa <- map_gapped_to_ungapped(seq_aligned = seq_aligned_aa[i3], gapped_pos = start_coord_aa, script = script)
                }
                if(is.na(end_coord_aa)){
                    end_coord_jalv_aa <- NA
                }else{
                    end_coord_jalv_aa <- map_gapped_to_ungapped(seq_aligned = seq_aligned_aa[i3], gapped_pos = end_coord_aa, script = script)
                }
                gff_aa_rows_jalv[[length(gff_aa_rows_jalv) + 1]] <- c(
                    df[selected_index[i3], Name],
                    ".",
                    get(paste0(i2, "_features"))[i4],
                    start_coord_jalv_aa, 
                    end_coord_jalv_aa, 
                    ".",
                    ".",
                    ".",
                    "."
                )
                # end aa coordinates
            }
        }
    }
    # nuc coordinates
    gff_table <- do.call(rbind, gff_rows)
    if(length(gff_table) > 0){
        gff_lines <- apply(gff_table, 1, function(x) paste(x, collapse="\t"))
    }else{
        gff_lines <- character()
    }
    gff_lines <- c("##gff-version 3", gff_lines)
    output_gff <- paste0(i2, "_", sub(x = basename(fasta_path), pattern = ".fasta$", replacement = ""), "_biojs.gff")
    writeLines(gff_lines, con = output_gff)
    # end nuc coordinates
    # aa coordinates
    gff_aa_table <- do.call(rbind, gff_aa_rows)
    if(length(gff_aa_table) > 0){
        gff_aa_lines <- apply(gff_aa_table, 1, function(x) paste(x, collapse="\t"))
    }else{
        gff_aa_lines <- character()
    }
    gff_aa_lines <- c("##gff_aa-version 3", gff_aa_lines)
    output_gff_aa <- paste0(i2, "_", sub(x = basename(fasta_path), pattern = "nuc.fasta$", replacement = ""), "aa_biojs.gff")
    writeLines(gff_aa_lines, con = output_gff_aa)


    gff_table_jalv <- do.call(rbind, gff_rows_jalv)
    if(length(gff_rows_jalv) > 0){
        gff_lines_jalv <- apply(gff_table_jalv, 1, function(x) paste(x, collapse="\t"))
    }else{
        gff_lines_jalv <- character()
    }
    tempo_text <- paste0(get(paste0(i2, "_features")), "\t", sub(x = get(paste0(i2, "_features_colors")), pattern = "#", replacement = ""))
    gff_lines_jalv <- c(tempo_text, "\n", "GFF", gff_lines_jalv)
    output_gff_jalv <- paste0(i2, "_sequence_alignment_with_gaps_imgt_nuc_jalview.gff")
    writeLines(gff_lines_jalv, con = output_gff_jalv)

    gff_aa_table_jalv <- do.call(rbind, gff_aa_rows_jalv)
    if(length(gff_aa_rows_jalv) > 0){
        gff_aa_lines_jalv <- apply(gff_aa_table_jalv, 1, function(x) paste(x, collapse="\t"))
    }else{
        gff_aa_lines_jalv <- character()
    }
    gff_aa_lines_jalv <- c(tempo_text, "\n", "GFF", gff_aa_lines_jalv)
    output_gff_aa_jalv <- paste0(i2, "_sequence_alignment_with_gaps_imgt_aa_jalview.gff")
    writeLines(gff_aa_lines_jalv, con = output_gff_aa_jalv)

    # end aa coordinates
}
# end Create the gff file





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




