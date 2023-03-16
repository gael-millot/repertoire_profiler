#!/usr/bin/Rscript
#########################################################################
##                                                                     ##
##     donut.R                                                     ##
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
    stop(paste0("\n\n================\n\nERROR IN donut.R\n", version$version.string, " IS NOT THE 4.1.2 RECOMMANDED\n\n================\n\n"))
}
# other initializations
erase.objects = TRUE # write TRUE to erase all the existing objects in R before starting the algorithm and FALSE otherwise. Beginners should use TRUE
if(erase.objects == TRUE){
    rm(list = ls(all.names = TRUE))
    erase.objects = TRUE
}
erase.graphs = TRUE # write TRUE to erase all the graphic windows in R before starting the algorithm and FALSE otherwise
script <- "donut"


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
        stop(paste0("\n\n================\n\nERROR IN donut.R\nTHE args OBJECT HAS NA\n\n================\n\n"), call. = FALSE)
    }
    tempo.arg.names <- c(
        "file_name", 
        "kind", 
        "donut.hole.size", 
        "donut.colors", 
        "cute", 
        "log"
    ) # objects names exactly in the same order as in the bash code and recovered in args. Here only one, because only the path of the config file to indicate after the donut.R script execution
    if(length(args) != length(tempo.arg.names)){
        stop(paste0("\n\n================\n\nERROR IN donut.R\nTHE NUMBER OF ELEMENTS IN args (", length(args),") IS DIFFERENT FROM THE NUMBER OF ELEMENTS IN tempo.arg.names (", length(tempo.arg.names),")\nargs:", paste0(args, collapse = ","), "\ntempo.arg.names:", paste0(tempo.arg.names, collapse = ","), "\n\n================\n\n"), call. = FALSE)
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
# file_name = "./all_productive_before_tree_seq.tsv"
# kind = "all"
# donut.hole.size = "2"
# donut.colors = "NULL"
# cute = "https://gitlab.pasteur.fr/gmillot/cute_little_R_functions/-/raw/v11.8.0/cute_little_R_functions.R"
# log = "all_donut.log"



################################ end Test

################################ Recording of the initial parameters


param.list <- c(
    "erase.objects", 
    "erase.graphs", 
    "script", 
    "run.way",
    if(run.way == "SCRIPT"){"command"}, 
    "file_name", 
    "kind", 
    "donut.hole.size", 
    "donut.colors", 
    "cute", 
    "log"
)
if(any(duplicated(param.list))){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 1 IN donut.R\nTHE param.list OBJECT CONTAINS DUPLICATED ELEMENTS:\n", paste(param.list[duplicated(param.list)], collapse = " "), "\n\n================\n\n"), call. = FALSE) # message for developers
}
if(erase.objects == TRUE){
    created.object.control <- ls()[ ! ls() %in% "param.list"]
    if( ! (all(created.object.control %in% param.list) & all(param.list %in% created.object.control))){
        stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 2 IN donut.R\nINCONSISTENCIES BETWEEN THE ARGUMENTS USED AND THE PARAMETERS REQUIRED IN THE EXECUTABLE CODE FILE\nTHE ARGUMENTS NOT PRESENT IN THE EXECUTABLE FILE (donut.R) ARE:\n", paste(created.object.control[ ! created.object.control %in% param.list], collapse = " "), "\nTHE PARAMETERS OF THE EXECUTABLE FILE (donut.R) NOT PRESENT IN THE ARGUMENTS ARE:\n", paste(param.list[ ! param.list %in% created.object.control], collapse = " "), "\n\n================\n\n"), call. = FALSE) # message for developers
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
    stop(paste0("\n\n============\n\nERROR IN donut.R\ncute PARAMETER MUST BE LENGTH 1: ", paste(cute, collapse = " "), "\n\n============\n\n"), call. = FALSE)
}else if(grepl(x = cute, pattern = "^http")){
    tempo.try <- try(suppressWarnings(suppressMessages(source(cute, local = .GlobalEnv))), silent = TRUE)
    if(any(grepl(x = tempo.try, pattern = "^[Ee]rror"))){
        stop(paste0("\n\n============\n\nERROR IN donut.R\nHTTP INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else if( ! grepl(x = cute, pattern = "^http")){
    if( ! file.exists(cute)){
        stop(paste0("\n\n============\n\nERROR IN donut.R\nFILE INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else{
    tempo.cat <- paste0("\n\n================\n\nINTERNAL CODE ERROR 3 IN donut.R: CODE HAS TO BE MODIFIED\n\n============\n\n")
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
    tempo.cat <- paste0("ERROR IN donut.R\nREQUIRED cute FUNCTION", ifelse(length(tempo) > 1, "S ARE", " IS"), " MISSING IN THE R ENVIRONMENT:\n", paste0(tempo, collapse = "()\n"))
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end required function checking


################ local function: package import


# R Packages required
req.package.list <- c(
    "lubridate", 
    "ggplot2"
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
tempo <- fun_check(data = file_name, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = kind, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = donut.hole.size, class = "vector", typeof = "character", length = 1) ; eval(ee)
tempo <- fun_check(data = donut.colors, class = "vector", typeof = "character", length = 1) ; eval(ee)
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
    "file_name", 
    "kind", 
    "donut.hole.size", 
    "donut.colors", 
    "cute", 
    "log"
)
tempo.log <- sapply(lapply(tempo.arg, FUN = get, env = sys.nframe(), inherit = FALSE), FUN = is.null)
if(any(tempo.log) == TRUE){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN donut.R:\n", ifelse(sum(tempo.log, na.rm = TRUE) > 1, "THESE ARGUMENTS\n", "THIS ARGUMENT\n"), paste0(tempo.arg[tempo.log], collapse = "\n"),"\nCANNOT BE NULL")
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

tempo <- fun_check(data = kind, options = c("all", "tree"), length = 1) ; eval(ee)

if(length(donut.hole.size) != 1 & any(grepl(donut.hole.size, pattern = "^\\-{0,1}[0-9]+\\.{0,1}[0-9]*$"))){ # positive numeric
    tempo.cat <- paste0("ERROR IN donut.R:\nTHE donut.hole.size PARAMETER MUST BE A SINGLE NUMBER\nHERE IT IS: \n", paste0(donut.hole.size, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}else{
    donut.hole.size <- as.numeric(donut.hole.size)
    tempo <- fun_check(data = donut.hole.size, class = "vector", mode = "numeric", length = 1) ; eval(ee)
}

if(donut.colors == "NULL"){
    donut.colors <- NULL
}else{
    tempo1 <- fun_check(data = donut.colors, class = "vector", mode = "character", length = 1)
    tempo2 <- fun_check(data = donut.colors, class = "vector", typeof = "integer", double.as.integer.allowed = TRUE, length = 1)
    checked.arg.names <- c(checked.arg.names, tempo2$object.name)
    if(tempo1$problem == TRUE & tempo2$problem == TRUE){
        tempo.cat <- paste0("ERROR IN donut.R\ndonut.colors ARGUMENT MUST BE (1) A HEXADECIMAL COLOR STRING STARTING BY #, OR (2) A COLOR NAME GIVEN BY colors(), OR (3) AN INTEGER VALUE")
        text.check2 <- c(text.check2, tempo.cat)
        arg.check2 <- c(arg.check2, TRUE)
    }else if(tempo1$problem == FALSE & tempo2$problem == TRUE){
        if( ! all(donut.colors %in% colors() | grepl(pattern = "^#", donut.colors), na.rm = TRUE)){
            tempo.cat <- paste0("ERROR IN donut.R\ndonut.colors ARGUMENT MUST BE (1) A HEXADECIMAL COLOR STRING STARTING BY #, OR (2) A COLOR NAME GIVEN BY colors(), OR (3) AN INTEGER VALUE")
            text.check2 <- c(text.check2, tempo.cat)
            arg.check2 <- c(arg.check2, TRUE)
        }
    }else{
        # no fun_check test here, it is just for checked.arg.names
        tempo <- fun_check(data = donut.colors, class = "vector")
        checked.arg.names <- c(checked.arg.names, tempo$object.name)
    }
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


fun_report(data = paste0("\n\n################################################################ donut PROCESS\n\n"), output = log, path = "./", overwrite = FALSE)
ini.date <- Sys.time()
ini.time <- as.numeric(ini.date) # time of process begin, converted into seconds
fun_report(data = paste0("\n\n################################ RUNNING DATE AND STARTING TIME\n\n"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0(ini.date, "\n\n"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0("\n\n################################ RUNNING\n\n"), output = log, path = "./", overwrite = FALSE)


################ End ignition


################ Graphical parameter initialization


################ End graphical parameter initialization





################ Data import

obs <- read.table(file_name, sep = "\t", header = TRUE)

################ End Data import

################ data modification

tempo.v <- obs$germline_v_call
tempo.j <- obs$germline_j_call
chain <- unique(substr(tempo.j, 1, 3)) # extract the IGH or IGK name
if(length(chain) != 1){
    stop(paste0("\n\n============\n\nINTERNAL CODE ERROR 1 IN donut.R for kind ", kind, ": chain MUST BE A SINGLE VALUE.\nHERE IT IS: ", paste(chain, collapse = " "), "\n\n============\n\n"), call. = FALSE)
}else if(chain != unique(substr(tempo.v, 1, 3))){
    stop(paste0("\n\n============\n\nINTERNAL CODE ERROR 2 IN donut.R for kind ", kind, ": CHAIN OF J DIFFERENT FROM CHAIN OF V.\nCHAIN OF V: ", paste(chain, collapse = " "), "\nCHAIN OF J: ", paste(unique(substr(tempo.v, 1, 3)), collapse = " "), "\n\n============\n\n"), call. = FALSE)
}
tempo.v <- substring(tempo.v, 4)
tempo.j <- substring(tempo.j, 4)
clone.name <- paste0(tempo.v, "_", tempo.j)
clone.id <-  obs$clone_id


obs2 <- data.frame(table(clone.name))
names(obs2)[1] <- "Germline_Clone_Name"
obs2 <- data.frame(obs2, Prop = obs2$Freq / sum(obs2$Freq))
obs3 <- data.frame(obs2, x = donut.hole.size)

################ End data modification

################ Plotting tree

tempo.gg.name <- "gg.indiv.plot."
tempo.gg.count <- 0
assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::ggplot(
    data = obs3,
    mapping = ggplot2::aes(x = x, y = Freq, fill = Germline_Clone_Name), 
    color = "white"
))
assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::geom_col(color = "white", size = 1.5))
# assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::geom_text(
#     ggplot2::aes(label = Freq), 
#     position = ggplot2::position_stack(vjust = 0.5)
# ))
assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::annotate(
    geom = "text", 
    x = 0.2, 
    y = 0, 
    label = sum(obs3$Freq), 
    size = 15
))
assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::coord_polar(theta = "y", direction = -1))
if( ! is.null(donut.colors)){
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::scale_fill_brewer(palette = donut.colors))
}
assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::xlim(c(0.2, donut.hole.size + 0.5)))
assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::theme_void())
tempo.plot <- eval(parse(text = paste(paste0(tempo.gg.name, 1:tempo.gg.count), collapse = " + ")))


tempo.title <- paste0(
    "Kind of sequences: ", kind, "\n",
    "Chain: ", chain, "\n",
    "The total number of sequences is indicated in the center of the donut.\nSee the donutchart.tsv file for the stats of the donut."
)

title.grob <- grid::textGrob(
    label = tempo.title,
    x = grid::unit(0, "lines"), 
    y = grid::unit(0, "lines"),
    hjust = 0,
    vjust = 0,
    gp = grid::gpar(fontsize = 8)
)
pdf(NULL)
tempo.plot <- gridExtra::arrangeGrob(tempo.plot, top = title.grob)


################ end Plotting tree

################ saving plots

ggplot2::ggsave(filename = paste0(kind, "_donutchart.png"), plot = tempo.plot, device = "png", path = ".", width = 5, height = 5, units = "in", dpi = 300)
ggplot2::ggsave(filename = paste0(kind, "_donutchart.svg"), plot = tempo.plot, device = "svg", path = ".", width = 5, height = 5, units = "in", dpi = 300)
ggplot2::ggsave(filename = paste0(kind, "_donutchart.pdf"), plot = tempo.plot, device = "pdf", path = ".", width = 5, height = 5, units = "in", dpi = 300)


################ end saving plots

################ data saving

write.table(obs2, file = paste0("./", kind, "_donutchart.tsv"), row.names = FALSE, sep = "\t")

################ End data saving

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
