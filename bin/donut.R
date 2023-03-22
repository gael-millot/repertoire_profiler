#!/usr/bin/Rscript
#########################################################################
##                                                                     ##
##     donut.R                                                         ##
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
    stop(paste0("\n\n================\n\nERROR IN donut.R\n\n\n", version$version.string, " IS NOT THE 4.1.2 RECOMMANDED\n\n================\n\n"))
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
        stop(paste0("\n\n================\n\nERROR IN donut.R\n\n\nTHE args OBJECT HAS NA\n\n================\n\n"), call. = FALSE)
    }
    tempo.arg.names <- c(
        "file_name", 
        "kind", 
        "donut.hole.size", 
        "donut.colors", 
        "donut_limit_legend",
        "cute", 
        "log"
    ) # objects names exactly in the same order as in the bash code and recovered in args. Here only one, because only the path of the config file to indicate after the donut.R script execution
    if(length(args) != length(tempo.arg.names)){
        stop(paste0("\n\n================\n\nERROR IN donut.R\n\n\nTHE NUMBER OF ELEMENTS IN args (", length(args),") IS DIFFERENT FROM THE NUMBER OF ELEMENTS IN tempo.arg.names (", length(tempo.arg.names),")\nargs:", paste0(args, collapse = ","), "\ntempo.arg.names:", paste0(tempo.arg.names, collapse = ","), "\n\n================\n\n"), call. = FALSE)
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
# file_name = "./productive_seq.tsv"
# kind = "all"
# donut.hole.size = "2"
# donut.colors = "NULL"
# donut_limit_legend = "0.1"
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
    "donut_limit_legend", 
    "cute", 
    "log"
)
if(any(duplicated(param.list))){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 1 IN donut.R\n\nTHE param.list OBJECT CONTAINS DUPLICATED ELEMENTS:\n", paste(param.list[duplicated(param.list)], collapse = " "), "\n\n================\n\n"), call. = FALSE) # message for developers
}
if(erase.objects == TRUE){
    created.object.control <- ls()[ ! ls() %in% "param.list"]
    if( ! (all(created.object.control %in% param.list) & all(param.list %in% created.object.control))){
        stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 2 IN donut.R\n\nINCONSISTENCIES BETWEEN THE ARGUMENTS USED AND THE PARAMETERS REQUIRED IN THE EXECUTABLE CODE FILE\nTHE ARGUMENTS NOT PRESENT IN THE EXECUTABLE FILE (donut.R) ARE:\n", paste(created.object.control[ ! created.object.control %in% param.list], collapse = " "), "\nTHE PARAMETERS OF THE EXECUTABLE FILE (donut.R) NOT PRESENT IN THE ARGUMENTS ARE:\n", paste(param.list[ ! param.list %in% created.object.control], collapse = " "), "\n\n================\n\n"), call. = FALSE) # message for developers
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
    stop(paste0("\n\n============\n\nERROR IN donut.R\n\ncute PARAMETER MUST BE LENGTH 1: ", paste(cute, collapse = " "), "\n\n============\n\n"), call. = FALSE)
}else if(grepl(x = cute, pattern = "^http")){
    tempo.try <- try(suppressWarnings(suppressMessages(source(cute, local = .GlobalEnv))), silent = TRUE)
    if(any(grepl(x = tempo.try, pattern = "^[Ee]rror"))){
        stop(paste0("\n\n============\n\nERROR IN donut.R\n\nHTTP INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else if( ! grepl(x = cute, pattern = "^http")){
    if( ! file.exists(cute)){
        stop(paste0("\n\n============\n\nERROR IN donut.R\n\nFILE INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
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
    tempo.cat <- paste0("ERROR IN donut.R\n\nREQUIRED cute FUNCTION", ifelse(length(tempo) > 1, "S ARE", " IS"), " MISSING IN THE R ENVIRONMENT:\n", paste0(tempo, collapse = "()\n"))
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end required function checking


################ local function: package import


# R Packages required
req.package.list <- c(
    "lubridate", 
    "ggplot2",
    "ggrepel"
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
tempo <- fun_check(data = donut_limit_legend, class = "vector", typeof = "character", length = 1) ; eval(ee)
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
    "donut_limit_legend", 
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
    tempo <- fun_check(data = donut.colors, options = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral", "Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3", "Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"), length = 1) ; eval(ee)
# for the moment it is inactivated below. But I can reactivate it to set specific colors
#    tempo1 <- fun_check(data = donut.colors, class = "vector", mode = "character", length = 1)
#    tempo2 <- fun_check(data = donut.colors, class = "vector", typeof = "integer", double.as.integer.allowed = TRUE, length = 1)
#    checked.arg.names <- c(checked.arg.names, tempo2$object.name)
#    if(tempo1$problem == TRUE & tempo2$problem == TRUE){
#        tempo.cat <- paste0("ERROR IN donut.R\n\ndonut.colors ARGUMENT MUST BE (1) A HEXADECIMAL COLOR STRING STARTING BY #, OR (2) A COLOR NAME GIVEN BY colors(), OR (3) AN INTEGER VALUE")
#        text.check2 <- c(text.check2, tempo.cat)
#        arg.check2 <- c(arg.check2, TRUE)
#    }else if(tempo1$problem == FALSE & tempo2$problem == TRUE){
#        if( ! all(donut.colors %in% colors() | grepl(pattern = "^#", donut.colors), na.rm = TRUE)){
#            tempo.cat <- paste0("ERROR IN donut.R\n\ndonut.colors ARGUMENT MUST BE (1) A HEXADECIMAL COLOR STRING STARTING BY #, OR (2) A COLOR NAME GIVEN BY colors(), OR (3) AN INTEGER VALUE")
#            text.check2 <- c(text.check2, tempo.cat)
#            arg.check2 <- c(arg.check2, TRUE)
#        }
#    }else{
#        # no fun_check test here, it is just for checked.arg.names
#        tempo <- fun_check(data = donut.colors, class = "vector")
#        checked.arg.names <- c(checked.arg.names, tempo$object.name)
#    }

}

if(donut_limit_legend == "NULL"){
    donut_limit_legend <- NULL
}else{
    if(length(donut_limit_legend) != 1 & any(grepl(donut_limit_legend, pattern = "(^{0,1}0+\\.*[0-9]*$)|(^1$)"))){ # positive numeric
        tempo.cat <- paste0("ERROR IN donut.R:\nTHE donut_limit_legend PARAMETER MUST BE A SINGLE NUMBER\nHERE IT IS: \n", paste0(donut_limit_legend, collapse = " "))
        text.check2 <- c(text.check2, tempo.cat)
        arg.check2 <- c(arg.check2, TRUE)
    }else{
        donut_limit_legend <- as.numeric(donut_limit_legend)
        tempo <- fun_check(data = donut_limit_legend, prop = TRUE, length = 1) ; eval(ee)
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
obs2 <- data.frame(obs2, Prop = obs2$Freq / sum(obs2$Freq), kind = kind)
obs2$Germline_Clone_Name <- factor(obs2$Germline_Clone_Name, levels = obs2$Germline_Clone_Name[order(obs2$Prop, decreasing = TRUE)]) # reorder so that the donut is according to decreasing proportion starting at the top in a clockwise direction
obs2 <- obs2[order(as.numeric(obs2$Germline_Clone_Name), decreasing = FALSE), ] # obs2 with rows in decreasing order, according to Prop
tempo.log <- obs[ , 1] == obs[ , 2]
if( ! all(tempo.log)){ # warning: I can do that because I have duplicated the first column of the dataset in the second column, in order to change the name in the first column with metadata
    tempo.labels <- obs[ , 1]
    tempo.labels[tempo.log] <- NA
    tempo.data1 <- aggregate(tempo.labels ~ clone.name, FUN = function(x){paste(x, collapse = ",")})
    obs2 <- data.frame(obs2, labels = tempo.data1$tempo.labels[match(obs2$Germline_Clone_Name, tempo.data1$clone.name)])
}
obs3 <- data.frame(obs2, x = 0) # staked bar at the origin of the donut set to x = donut.hole.size


################ End data modification

################ Plotting tree

tempo.gg.name <- "gg.indiv.plot."
tempo.gg.count <- 0
assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::ggplot(
    data = obs3,
    mapping = ggplot2::aes(x = x, y = Freq, fill = Germline_Clone_Name), 
    color = "white"
))
bar_width = 1
assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::geom_col(color = "white", size = 1.5, width = bar_width)) # size is size of the separation
# assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::geom_text(
#     ggplot2::aes(label = Freq), 
#     position = ggplot2::position_stack(vjust = 0.5)
# ))
assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::xlim(c(- bar_width / 2 - donut.hole.size, bar_width))) # must be centered on x = 2
assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::ylim(c(0, max(cumsum(obs3$Freq)))))
assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::annotate(
    geom = "text", 
    x = - bar_width / 2 - donut.hole.size, 
    y = 0, 
    label = sum(obs3$Freq), 
    size = 15
))
assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::coord_polar(theta = "y", direction = -1, start = 0, clip = "on"))
if( ! is.null(donut.colors)){
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::scale_fill_brewer(palette = donut.colors))
}
assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::theme_void())
assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::theme(legend.title=element_blank()))


tempo.log <- obs[ , 1] == obs[ , 2]
if( ! all(tempo.log)){ # warning: I can do that because I have duplicated the first column of the dataset in the second column, in order to change the name in the first column with metadata
    tempo <- rev(cumsum(rev(obs3$Freq)))
    obs4 <- data.frame(obs3, text_y = tempo - (tempo - c(tempo[-1], 0)) / 2)
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggrepel::geom_text_repel(
        data = obs4, 
        mapping = ggplot2::aes(
            x = x, 
            y = text_y, 
            label = labels
        ), 
        size = 3, 
        force_pull = 100, 
        nudge_x = bar_width / 2 + (bar_width - bar_width / 2) / 2, # add nudge_x to the center of the bar
        show.legend = FALSE
    ))
}


if( ! is.null(donut_limit_legend)){
    if(sum(obs3$Prop >= donut_limit_legend) == 0){
        tempo.cat <- paste0("ERROR IN donut.R for kind ", kind, ":\n\nTHE donut_limit_legend PARAMETER VALUE (", donut_limit_legend, ") IS TOO HIGH FOR THE PROPORTIONS IN THE DONUT PLOT:\n", paste0(obs3$Prop, collapse = "\n"))
        stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) 
    }else{
        assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::scale_fill_discrete(
            breaks = as.character(obs3$Germline_Clone_Name[obs3$Prop >= donut_limit_legend])
        ))
    }
}


tempo.plot <- eval(parse(text = paste(paste0(tempo.gg.name, 1:tempo.gg.count), collapse = " + ")))
tempo.title <- paste0(
    "Kind of sequences: ", kind, " (see the corresponding ", ifelse(kind == "all", "productive_seq.tsv", ifelse(kind == "tree", "seq_for_trees.tsv", "=====ERROR=====")), " file)\n",
    "Chain: ", chain, "\n",
    "The total number of sequences is indicated in the center of the donut.\nSee the donut_stats.tsv file for the stats of the donut."
)
if( ! is.null(donut_limit_legend)){
    if(sum(obs3$Prop >= donut_limit_legend) < length(obs3$Prop)){
        tempo.title <- paste0(tempo.title, "\nClasses have been removed from the legend (donut_limit_legend = ", round(donut_limit_legend, 2), "). See the donut_stats.tsv file")
    }
}

title.grob <- grid::textGrob(
    label = tempo.title,
    x = grid::unit(0, "lines"), 
    y = grid::unit(0, "lines"),
    hjust = 0,
    vjust = 0,
    gp = grid::gpar(fontsize = 7)
)
pdf(NULL)
tempo.plot <- suppressMessages(suppressWarnings(gridExtra::arrangeGrob(tempo.plot, top = title.grob, left = " ", right = " "))) # , left = " ", right = " " : trick to add margins in the plot. padding =  unit(0.5, "inch") is for top margin above the title


################ end Plotting tree

################ saving plots

ggplot2::ggsave(filename = paste0(kind, "_donutchart.png"), plot = tempo.plot, device = "png", path = ".", width = 5, height = 5, units = "in", dpi = 300)
ggplot2::ggsave(filename = paste0(kind, "_donutchart.svg"), plot = tempo.plot, device = "svg", path = ".", width = 5, height = 5, units = "in", dpi = 300)
ggplot2::ggsave(filename = paste0(kind, "_donutchart.pdf"), plot = tempo.plot, device = "pdf", path = ".", width = 5, height = 5, units = "in", dpi = 300)


################ end saving plots

################ data saving

write.table(obs2, file = paste0("./", kind, "_donut_stats.tsv"), row.names = FALSE, sep = "\t")

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
    tempo.cat <- paste0("IN donut.R OF THE NEXFLOW EXECUTION:\n\n", warn)
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
