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
        "igblast_data_check_ch", 
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


# setwd("Z:/thomas_derenne/repertoire_profiler-master/work/6b/2d1f7c09424fd2b0ff522192dd3c9d")
# all_passed_seq <- "productive_seq.tsv"
# igblast_data_check_ch <- "imgt_human_IGHD.tsv imgt_human_IGHJ.tsv imgt_human_IGHV.tsv imgt_human_IGHC.tsv"
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
    "productive_seq", 
    "igblast_data_check_ch", 
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
    "igblast_data_check_ch", 
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

rep_file_names <- strsplit(igblast_data_check_ch, split = " ")[[1]]
if(length(rep_file_names) == 0){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 4 IN repertoire.R:\nPROBLEM WITH igblast_data_check_ch: ", paste(igblast_data_check_ch, collapse = " "), "PLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
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

# names of the columns to deal with in single allele repertoire
var_allele_obs <- c("v_call", "j_call")
const_allele_obs <- c( "c_call") 
allele_obs <- c(var_allele_obs, const_allele_obs)

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

# processing of the imgt database files
tempo <- strsplit(rep_file_names, split = "[_.]")
allele <- vector(mode = "list", length = length(rep_file_names)) # list that will contain all the imgt alleles
for(i1 in 1:length(allele)){
    tempo <- strsplit(rep_file_names[i1], split = "[_.]")[[1]]
    names(allele)[i1] <- tempo[length(tempo) - 1]
    allele[[i1]] <- scan(rep_file_names[i1], what = "character", quiet = TRUE)
}



# end processing of the imgt database files

################ End Data import


################ data verification

if( ! all(allele_obs %in% c("v_call", "j_call", "c_call"))){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 5 IN repertoire.R:\nPROBLEM WITH THE allele_obs INTERNAL VARIABLE THAT MUST BE \"v_call\", \"j_call\", \"c_call\"\nHERE IT IS:\n", paste(allele_obs, collapse = "\n"), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
}

if( ! all(allele_obs %in% names(df))){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 7 IN repertoire.R:\nPROBLEM WITH THE NAMES OF productive_seq THAT MUST CONTAIN \"v_call\", \"j_call\", \"c_call\"\nHERE IT IS:\n", paste(names(df), collapse = "\n"), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
}

if(is.null(names(allele))){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 9 IN repertoire.R:\nNO NAMES OF allele CAN BE NULL:", paste(names(allele), collapse = "\n"), "\nSEE rep_file_names:\n", paste(rep_file_names, collapse = "\n"), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
}


################ end data verification


################ data modification, plotting and saving

######## observed data

if(grepl(x = names(df)[2], pattern = "^initial_")){ # means that fonctional annotations are present
    annotation.log <- df[ , 1] == df[ , 2] 
    all.annotation.log <- all(annotation.log) # if one difference between df[ , 1] == df[ , 2], then all.annotation.log is FALSE
}else{
    all.annotation.log <- TRUE
}
if(all.annotation.log == FALSE){
    kind <- c(kind, "annotated")
}

# take the first annotation if several in the tempo column
if(productive_seq == "productive_seq.tsv"){
    tempo.several.annot.log <- logical(length = nrow(df)) # only false
    tempo <- c(allele_obs, gene_obs)
    for(i1 in 1:length(tempo)){
        col_data <- df[ , names(df) == tempo[i1]]
        tempo.nb <- sapply(strsplit(as.character(col_data), split = ","), FUN = function(x){
            if (all(is.na(x))) return(0) else return(length(x))
        })
        tempo.several.annot.log[tempo.nb > 1] <- TRUE
    }
    if(any(tempo.several.annot.log)){
        write.table(df[tempo.several.annot.log, ], file = paste0("./productive_seq_several_annot_igmt.tsv"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE) # separate repertoires
    }else{
        write.table("", file = paste0("./productive_seq_several_annot_igmt.tsv"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    }
    for(i1 in 1:length(tempo)){
        col_data <- df[ , names(df) == tempo[i1]]
        df[ , names(df) == tempo[i1]] <- sapply(strsplit(as.character(col_data), split = ","), FUN = function(x){
            if (all(is.na(x))) return(NA) else return(x[1])
        })
    }
}
# end take the first annotation if several in the tempo column


######## observed data

######## removal of the common character on the left of strings

# in imgt files
tempo <- left_common_chars(s = allele)
allele_common <- tempo$common
allele_trunk <- tempo$rest
tempo <- left_common_chars(s = gene)
gene_common <- tempo$common
gene_trunk <- tempo$rest
if(allele_common != gene_common){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 11 IN repertoire.R:\nallele_common AND gene_common SHOULD BE IDENTICAL.\nHERE allele_common IS:\n", allele_common, "\nAND gene_common IS:\n", gene_common, "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
}
# concordance
pos <- NULL
names_allele <- names(allele_trunk)

for(i0 in 1:length(check_concordance_imgt)){
    tempo_log <- names_allele %in% check_concordance_imgt[[i0]]
    # Check that a concordance is found
    if(sum(tempo_log) != 1){ 
        stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 12 IN repertoire.R:\nnames_allele NOT PRESENT IN check_concordance_imgt NUMBER ", i0, ".\nHERE names_allele IS:\n", paste(names_allele, collapse = "\n"), "\nAND check_concordance_imgt NUMBER ", i0, " IS:\n", paste(check_concordance_imgt[[i0]], collapse = "\n"), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
    }else{
        pos <- c(pos, which(tempo_log))
    }
}
allele_trunk <- allele_trunk[pos]
gene_trunk <- gene_trunk[pos]
# end concordance
# end in imgt files

# in data
tempo_list <- as.list(df[, allele_obs])
tempo <- left_common_chars(s = tempo_list)
allele_obs_common <- tempo$common
allele_obs_trunk <- tempo$rest
tempo_list <- as.list(df[, gene_obs])
tempo <- left_common_chars(s = tempo_list)
gene_obs_common <- tempo$common
gene_obs_trunk <- tempo$rest
if( ! (allele_obs_common == allele_common & allele_obs_common == gene_obs_common)){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 15 IN repertoire.R:\nallele_obs_common, gene_obs_common AND allele_common SHOULD BE IDENTICAL.\nHERE allele_obs_common IS:\n", allele_obs_common, "\ngene_obs_common IS:\n", gene_obs_common, "\nAND allele_common IS:\n", allele_common, "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE)
}
# end data


######## end removal of the common character on the left of strings


# Chain is determined with igblast_ref_files because if we studied IGL and IGK but only one was found, both still need to be in the title
loci <- sapply(igblast_data_check_ch, extract_chain)
chain <- unique(loci[!is.na(loci)]) # extract the IGH or IGK name (ignoring NA values)

# simple repertoire
type <- c("allele", "gene")
for(i0 in type){
    # concordance between imgt and obs data already made above
    # simple repertoires
    for(i1 in 1:length(get(paste0(i0, "_obs_trunk")))){
        tempo.df <- get(paste0(i0, "_obs_trunk"))[[i1]]
        tempo.df <- tempo.df[ ! is.na(tempo.df)]
        tempo.df <- factor(tempo.df, levels = unique(get(paste0(i0, "_trunk"))[[i1]]))
        # plot
        for(i2 in kind){
            if(i2 == "non-zero"){
                tempo.table <- table(tempo.df)
                tempo.table.gg <- tempo.table[tempo.table > 0]
            }else if(i2 == "annotated"){
                tempo.table <- table(tempo.df[ ! annotation.log])
                tempo.table.gg <- tempo.table[tempo.table > 0]
            }else{
                tempo.table <- table(tempo.df)
                write.table(tempo.table, file = paste0("./rep_", i0, "_", names(get(paste0(i0, "_trunk")))[i1], ".tsv"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) # separate repertoires
                tempo.table.gg <- tempo.table
            }
            if(sum(tempo.table.gg, na.rm = TRUE) > 0){
                if(length(tempo.table.gg) == 1){
                    tempo.table.gg <- data.frame(Var1 = names(tempo.table.gg), Count = tempo.table.gg, row.names = NULL)
                    names(tempo.table.gg)[1] <- i0
                }else{
                    tempo.table.gg <- as.data.frame(tempo.table.gg)
                    names(tempo.table.gg)[1] <- i0
                    names(tempo.table.gg)[names(tempo.table.gg) == "Freq"] <- "Count"
                }
                tempo.table.gg$Count[tempo.table.gg$Count == 0] <- NA
                tempo.title <- paste0(
                    "Locus: ", names(check_concordance_imgt)[i1], "\n",
                    "Kind: ", i2, "\n",
                    "Type: ", i0, "\n",
                    "Chain: ", paste(chain, collapse = " "), "\n",
                    "Total count: ", sum(tempo.table.gg$Count, na.rm = TRUE)
                )
                # label.size <- 25 / log2(nrow(tempo.table.gg) + 1) # use https://www.wolframalpha.com/widgets/view.jsp?id=f995c9aeb1565edd78adb37d2993d66
                label.size <- 25 - nrow(tempo.table.gg)*0.3 
                final.plot <- fun_gg_heatmap2(
                    data1 = tempo.table.gg,
                    x = NULL,
                    y = i0,
                    z = "Count",
                    label.size = ifelse(label.size <= 1, 1, label.size),
                    color.low = "white",
                    color.high = "blue",
                    zero.color = grey(0.95),
                    cell.value = TRUE,
                    cell.value.size = ifelse(label.size <= 0, 1, label.size) / 3,
                    title = tempo.title,
                    title.size = 12
                )
            }else{
                # no need to use pdf(NULL) with fun_gg_empty_graph()
                final.plot <- fun_gg_empty_graph(text = paste0("NO GRAPH PLOTTED FOR ", names(get(paste0(i0, "_trunk")))[i1], "\nNO ALLELE/GENE DETECTED"), text.size = 3)
            }
            ggplot2::ggsave(filename = paste0(names(get(paste0(i0, "_trunk")))[i1], "_", i0, "_", i2, ".png"), plot = final.plot, device = "png", path = ".", width = 4, height = 10, units = "in", dpi = 300, bg = "white") # do not modify width and height. Otherwise impair axis.text.y, axis.ticks.y, panel.border sizes
            ggplot2::ggsave(filename = paste0(names(get(paste0(i0, "_trunk")))[i1], "_", i0, "_", i2, ".svg"), plot = final.plot, device = "svg", path = ".", width = 4, height = 10, units = "in", dpi = 300, bg = "white")
            ggplot2::ggsave(filename = paste0(names(get(paste0(i0, "_trunk")))[i1], "_", i0, "_", i2, ".pdf"), plot = final.plot, device = "pdf", path = ".", width = 4, height = 10, units = "in", dpi = 300, bg = "white")
        }
        # end plot
    }
    # end simple repertoires
    # combined repertoires
    # V x J
    tempo.df <- get(paste0(i0 , "_obs_trunk"))[get(paste0("var_", i0, "_obs"))]
    tempo.df[[1]] <- factor(tempo.df[[1]], levels = unique(get(paste0(i0, "_trunk"))[[1]]))
    tempo.df[[2]] <- factor(tempo.df[[2]], levels = unique(get(paste0(i0, "_trunk"))[[2]]))
    # plot
    for(i1 in kind){
        # here the work is using data frames because it keeps the structure even if one cell
        if(i1 == "non-zero"){
            tempo.table2 <- as.data.frame.matrix(table(tempo.df))
            tempo.log <- apply(tempo.table2, 1, sum, na.rm = TRUE) > 0
            tempo.table3 <- tempo.table2[tempo.log, ]
            tempo.log <- apply(tempo.table3, 2, sum, na.rm = TRUE) > 0
            tempo.table3 <- tempo.table3[tempo.log]
        }else if(i1 == "annotated"){
            tempo <- lapply(X = tempo.df, FUN = function(x){x[ ! annotation.log]})
            tempo.table2 <- as.data.frame.matrix(table(tempo))
            tempo.log <- apply(tempo.table2, 1, sum, na.rm = TRUE) > 0
            tempo.table3 <- tempo.table2[tempo.log, ]
            tempo.log <- apply(tempo.table3, 2, sum, na.rm = TRUE) > 0
            tempo.table3 <- tempo.table3[tempo.log]
        }else{
            tempo.table2 <- as.data.frame.matrix(table(tempo.df))
            write.table(tempo.table2, file = paste0("./rep_", i0, "_", names(get(paste0(i0, "_trunk")))[1], "_", names(get(paste0(i0, "_trunk")))[2], ".tsv"), row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE) # separate repertoires
            tempo.table3 <- tempo.table2
        }
        # end here the work is using data frames because it keeps the structure even if one cell
        if(sum(tempo.table3, na.rm = TRUE) > 0){
            if(nrow(tempo.table3) > 1 & ncol(tempo.table3) > 1){
                tempo.table.gg <- as.data.frame(as.table(as.matrix(tempo.table3)))
                label.size <- (25 - nrow(tempo.table3) * 0.3) / max(1, 0.5 * ncol(tempo.table3))
            }else{
                tempo.table.gg <- data.frame(Var1 = rownames(tempo.table3), Var2 = colnames(tempo.table3), Freq = tempo.table3, row.names = NULL)
                label.size <- 25 - nrow(tempo.table.gg) * 0.3
            }
            names(tempo.table.gg) <- c(names(check_concordance_imgt)[1], names(check_concordance_imgt)[2], "Count")
            tempo.table.gg$Count[tempo.table.gg$Count == 0] <- NA
            tempo.title <- paste0(
                "Locus: ", names(check_concordance_imgt)[1], "x", names(check_concordance_imgt)[2], "\n",
                "Kind: ", i1, "\n",
                "Type: ", i0, "\n",
                "Chain: ", paste(chain, collapse = " "), "\n",
                "Total count: ", sum(tempo.table.gg$Count, na.rm = TRUE)
            )
            #label.size <- -5/132 * nrow(tempo.table.gg) + 1081/66 # use https://www.wolframalpha.com/widgets/view.jsp?id=f995c9aeb1565edd78adb37d2993d66
            
            final.plot <- fun_gg_heatmap2(
                data1 = tempo.table.gg,
                x = names(check_concordance_imgt)[2],
                y = names(check_concordance_imgt)[1],
                z = "Count",
                label.size = ifelse(label.size <= 1, 1, label.size),
                color.low = "white",
                color.high = "blue",
                zero.color = grey(0.95),
                cell.value = TRUE,
                cell.value.size = ifelse(label.size <= 0, 1, label.size) / 3,
                title = tempo.title,
                title.size = 12
            )
        }else{
            # no need to use pdf(NULL) with fun_gg_empty_graph()
            final.plot <- fun_gg_empty_graph(text = paste0("NO GRAPH PLOTTED FOR ", names( get(paste0(i0, "_trunk")))[i1], "\nNO ALLELE/GENE DETECTED"), text.size = 3)
        }
        ggplot2::ggsave(filename = paste0("./rep_", i0, "_", names(get(paste0(i0, "_trunk")))[1], "_", names(get(paste0(i0, "_trunk")))[2], "_",  i1, ".png"), plot = final.plot, device = "png", path = ".", width = 4, height = 10, units = "in", dpi = 300, bg = "white") # do not modify width and height. Otherwise impair axis.text.y, axis.ticks.y, panel.border sizes
        ggplot2::ggsave(filename = paste0("./rep_", i0, "_", names(get(paste0(i0, "_trunk")))[1], "_", names(get(paste0(i0, "_trunk")))[2], "_",  i1, ".svg"), plot = final.plot, device = "svg", path = ".", width = 4, height = 10, units = "in", dpi = 300, bg = "white")
        ggplot2::ggsave(filename = paste0("./rep_", i0, "_", names(get(paste0(i0, "_trunk")))[1], "_", names(get(paste0(i0, "_trunk")))[2], "_",  i1, ".pdf"), plot = final.plot, device = "pdf", path = ".", width = 4, height = 10, units = "in", dpi = 300, bg = "white")
    }
    # end plot
    # end V x J
    # for each C
    isotype_subclass_obs <- sort(unique(gene_obs_trunk[[3]]))
    for(i1 in isotype_subclass_obs){
        tempo.df <- lapply(X = get(paste0(i0 , "_obs_trunk"))[get(paste0("var_", i0, "_obs"))], 
                           FUN = function(x) {
                               c_values <- get(paste0(i0, "_obs_trunk"))[[get(paste0("const_", i0, "_obs"))]]
                               c_genes <- lapply(strsplit(c_values, "\\*"), function(x) x[[1]])
                               x[c_genes == i1]
                           })
        
        tempo.df[[1]] <- factor(tempo.df[[1]], levels = unique(get(paste0(i0, "_trunk"))[[1]]))
        tempo.df[[2]] <- factor(tempo.df[[2]], levels = unique(get(paste0(i0, "_trunk"))[[2]]))
        # plot
        for(i2 in kind){
            # here the work is using data frames because it keeps the structure even if one cell
            if(i2 == "non-zero"){
                tempo.table2 <- as.data.frame.matrix(table(tempo.df))
                tempo.log <- apply(tempo.table2, 1, sum, na.rm = TRUE) > 0
                tempo.table3 <- tempo.table2[tempo.log, ]
                tempo.log <- apply(tempo.table3, 2, sum, na.rm = TRUE) > 0
                tempo.table3 <- tempo.table3[tempo.log]
            }else if(i2 == "annotated"){
                tempo <- lapply(X = tempo.df, FUN = function(x){x[ ! annotation.log]})
                tempo.table2 <- as.data.frame.matrix(table(tempo))
                tempo.log <- apply(tempo.table2, 1, sum, na.rm = TRUE) > 0
                tempo.table3 <- tempo.table2[tempo.log, ]
                tempo.log <- apply(tempo.table3, 2, sum, na.rm = TRUE) > 0
                tempo.table3 <- tempo.table3[tempo.log]
            }else{
                tempo.table2 <- as.data.frame.matrix(table(tempo.df))
                write.table(tempo.table2, file = paste0("./rep_", i0, "_", i1, "_", names(get(paste0(i0, "_trunk")))[1], "_", names(get(paste0(i0, "_trunk")))[2], ".tsv"), row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE) # separate repertoires
                tempo.table3 <- tempo.table2
            }
            # end here the work is using data frames because it keeps the structure even if one cell
            if(sum(tempo.table3, na.rm = TRUE) > 0){
                if(nrow(tempo.table3) > 1 & ncol(tempo.table3) > 1){
                    tempo.table.gg <- as.data.frame(as.table(as.matrix(tempo.table3)))
                    label.size <- (25 - nrow(tempo.table3) * 0.3) / max(1, 0.5 * ncol(tempo.table3))
                }else{
                    tempo.table.gg <- data.frame(Var1 = rownames(tempo.table3), Var2 = colnames(tempo.table3), Freq = tempo.table3, row.names = NULL)
                    label.size <- 25 - nrow(tempo.table.gg)*0.3 
                }
                names(tempo.table.gg) <- c(names(check_concordance_imgt)[1], names(check_concordance_imgt)[2], "Count")
                tempo.table.gg$Count[tempo.table.gg$Count == 0] <- NA
                tempo.title <- paste0(
                    "Locus: ", names(check_concordance_imgt)[1], "x", names(check_concordance_imgt)[2], " for ", i1, "\n",
                    "Kind: ", i2, "\n",
                    "Type: ", i0, "\n",
                    "Chain: ", paste(chain, collapse = " "), "\n",
                    "Total count: ", sum(tempo.table.gg$Count, na.rm = TRUE)
                )
                # label.size <- -5/132 * nrow(tempo.table.gg) + 1081/66 # use https://www.wolframalpha.com/widgets/view.jsp?id=f995c9aeb1565edd78adb37d2993d66
                final.plot <- fun_gg_heatmap2(
                    data1 = tempo.table.gg,
                    x = names(check_concordance_imgt)[2],
                    y = names(check_concordance_imgt)[1],
                    z = "Count",
                    label.size = ifelse(label.size <= 1, 1, label.size),
                    color.low = "white",
                    color.high = "blue",
                    zero.color = grey(0.95),
                    cell.value = TRUE,
                    cell.value.size = ifelse(label.size <= 0, 1, label.size) / 3,
                    title = tempo.title,
                    title.size = 12
                )
            }else{
                # no need to use pdf(NULL) with fun_gg_empty_graph()
                final.plot <- fun_gg_empty_graph(title = tempo.title, title.size = 12, text = paste0("NO GRAPH PLOTTED FOR VxJ for ", i1, "\nNO ALLELE/GENE DETECTED"), text.size = 3)
            }
            ggplot2::ggsave(filename = paste0("./rep_", i0, "_for_", i1, "_", names(get(paste0(i0, "_trunk")))[1], "_", names(get(paste0(i0, "_trunk")))[2], "_",  i2, ".png"), plot = final.plot, device = "png", path = ".", width = 4, height = 10, units = "in", dpi = 300, bg = "white") # do not modify width and height. Otherwise impair axis.text.y, axis.ticks.y, panel.border sizes
            ggplot2::ggsave(filename = paste0("./rep_", i0, "_for_", i1, "_", names(get(paste0(i0, "_trunk")))[1], "_", names(get(paste0(i0, "_trunk")))[2], "_",  i2, ".svg"), plot = final.plot, device = "svg", path = ".", width = 4, height = 10, units = "in", dpi = 300, bg = "white")
            ggplot2::ggsave(filename = paste0("./rep_", i0, "_for_", i1, "_", names(get(paste0(i0, "_trunk")))[1], "_", names(get(paste0(i0, "_trunk")))[2], "_",  i2, ".pdf"), plot = final.plot, device = "pdf", path = ".", width = 4, height = 10, units = "in", dpi = 300, bg = "white")
        }
    }
    # end plot
    # end for each C
    # end combined repertoires
    files <- list.files(path = ".", pattern = paste0("^.*", i0, ".*\\.pdf$"))
    sorted_files <- files[order(grepl("for", files))]
    tempo <- qpdf::pdf_combine(input = sorted_files, output = paste0("./", i0, "_repertoire.pdf")) # assignation to prevent a returned element
    
}

################ end data modification, plotting and saving


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

