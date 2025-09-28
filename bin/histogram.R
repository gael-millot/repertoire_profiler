#!/usr/bin/Rscript
#########################################################################
##                                                                     ##
##     histogram.R                                                     ##
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
if(version$version.string != "R version 4.3.2 (2023-10-31)"){
    stop(paste0("\n\n================\n\nERROR IN histogram.R\n", version$version.string, " IS NOT THE 4.3.2 RECOMMANDED\n\n================\n\n"))
}
# other initializations
erase.objects = TRUE # write TRUE to erase all the existing objects in R before starting the algorithm and FALSE otherwise. Beginners should use TRUE
if(erase.objects == TRUE){
    rm(list = ls(all.names = TRUE))
    erase.objects = TRUE
}
erase.graphs = TRUE # write TRUE to erase all the graphic windows in R before starting the algorithm and FALSE otherwise
script <- "histogram"


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
        stop(paste0("\n\n================\n\nERROR IN histogram.R\nTHE args OBJECT HAS NA\n\n================\n\n"), call. = FALSE)
    }
    tempo.arg.names <- c(
        "file.name", 
        "clone_model", 
        "clone_normalize", 
        "clone_distance", 
        "cute", 
        "log"
    ) # objects names exactly in the same order as in the bash code and recovered in args. Here only one, because only the path of the config file to indicate after the histogram.R script execution
    if(length(args) != length(tempo.arg.names)){
        stop(paste0("\n\n================\n\nERROR IN histogram.R\nTHE NUMBER OF ELEMENTS IN args (", length(args),") IS DIFFERENT FROM THE NUMBER OF ELEMENTS IN tempo.arg.names (", length(tempo.arg.names),")\nargs:", paste0(args, collapse = ","), "\ntempo.arg.names:", paste0(tempo.arg.names, collapse = ","), "\n\n================\n\n"), call. = FALSE)
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

# setwd("C:/Users/gael/Documents/Git_projects/repertoire_profiler/work/72/b43d4e26736d47baef36ed10330319")
# file.name <- "nearest_distance.tsv"
# clone_model <- "ham"
# clone_normalize <- "len"
# clone_distance <- "0.15"
# cute = "https://gitlab.pasteur.fr/gmillot/cute_little_R_functions/-/raw/v12.3.0/cute_little_R_functions.R"
# log = "histogram.log"





################################ end Test

################################ Recording of the initial parameters


param.list <- c(
    "erase.objects", 
    "erase.graphs", 
    "script", 
    "run.way",
    "tempo.arg.names", 
    if(run.way == "SCRIPT"){"command"}, 
    "file.name", 
    "clone_model", 
    "clone_normalize", 
    "clone_distance", 
    "cute", 
    "log"
)
if(any(duplicated(param.list))){
    stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 1 IN histogram.R\nTHE param.list OBJECT CONTAINS DUPLICATED ELEMENTS:\n", paste(param.list[duplicated(param.list)], collapse = " "), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE) # message for developers
}
if(erase.objects == TRUE){
    created.object.control <- ls()[ ! ls() %in% "param.list"]
    if( ! (all(created.object.control %in% param.list) & all(param.list %in% created.object.control))){
        stop(paste0("\n\n================\n\nINTERNAL CODE ERROR 2 IN histogram.R\nINCONSISTENCIES BETWEEN THE ARGUMENTS USED AND THE PARAMETERS REQUIRED IN THE EXECUTABLE CODE FILE\nTHE ARGUMENTS NOT PRESENT IN THE EXECUTABLE FILE (histogram.R) ARE:\n", paste(created.object.control[ ! created.object.control %in% param.list], collapse = " "), "\nTHE PARAMETERS OF THE EXECUTABLE FILE (histogram.R) NOT PRESENT IN THE ARGUMENTS ARE:\n", paste(param.list[ ! param.list %in% created.object.control], collapse = " "), "\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n================\n\n"), call. = FALSE) # message for developers
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
    stop(paste0("\n\n============\n\nERROR IN histogram.R\ncute PARAMETER MUST BE LENGTH 1: ", paste(cute, collapse = " "), "\n\n============\n\n"), call. = FALSE)
}else if(grepl(x = cute, pattern = "^http")){
    tempo.try <- try(suppressWarnings(suppressMessages(source(cute, local = .GlobalEnv))), silent = TRUE)
    if(any(grepl(x = tempo.try, pattern = "^[Ee]rror"))){
        stop(paste0("\n\n============\n\nERROR IN histogram.R\nHTTP INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else if( ! grepl(x = cute, pattern = "^http")){
    if( ! file.exists(cute)){
        stop(paste0("\n\n============\n\nERROR IN histogram.R\nFILE INDICATED IN THE cute PARAMETER DOES NOT EXISTS: ", cute, "\n\n============\n\n"), call. = FALSE)
    }else{
        source(cute, local = .GlobalEnv) # source the fun_ functions used below
    }
}else{
    tempo.cat <- paste0("\n\n================\n\nINTERNAL CODE ERROR 3 IN histogram.R: CODE HAS TO BE MODIFIED\nPLEASE, SEND AN ISSUE AT https://gitlab.pasteur.fr/gmillot/repertoire_profiler OR REPORT AT gael.millot@pasteur.fr\n\n============\n\n")
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
    tempo.cat <- paste0("ERROR IN histogram.R\nREQUIRED cute FUNCTION", ifelse(length(tempo) > 1, "S ARE", " IS"), " MISSING IN THE R ENVIRONMENT:\n", paste0(tempo, collapse = "()\n"))
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end required function checking


################ local function: package import


# R Packages required
req.package.list <- c(
    "lubridate", 
    "ggplot2", 
    "shazam"
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
# management of NULL arguments, WARNING: only for histogram.R because NULL is "NULL" in the nextflow.config file
tempo.arg <-c(
    "file.name", 
    "clone_model", 
    "clone_normalize", 
    "clone_distance", 
    "log"
)
tempo.log <- sapply(lapply(tempo.arg, FUN = get, env = sys.nframe(), inherit = FALSE), FUN = is.null)
if(any(tempo.log) == TRUE){# normally no NA with is.null()
    tempo.cat <- paste0("ERROR IN histogram.R:\n", ifelse(sum(tempo.log, na.rm = TRUE) > 1, "THESE ARGUMENTS\n", "THIS ARGUMENT\n"), paste0(tempo.arg[tempo.log], collapse = "\n"),"\nCANNOT BE NULL")
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n"), call. = FALSE) # == in stop() to be able to add several messages between ==
}
# end management of NULL arguments, WARNING: only for histogram.R because NULL is "NULL" in the nextflow.config file
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
tempo <- fun_check(data = clone_model, options = c("ham", "aa", "hh_s1f", "hh_s5f", "mk_rs1nf", "mk_rs5nf", "m1n_compat", "hs1f_compat"), length = 1) ; eval(ee)
tempo <- fun_check(data = clone_normalize, options = c("len", "none"), length = 1) ; eval(ee)
if(length(clone_distance) != 1 & any(grepl(clone_distance, pattern = "^(1)|(0)|(0\\.[0-9]*)$"))){ # positive prop
    tempo.cat <- paste0("ERROR IN histogram.R:\nTHE clone_distance PARAMETER MUST BE A SINGLE POSITIVE PROPORTION\nHERE IT IS: \n", paste0(clone_distance, collapse = " "))
    text.check2 <- c(text.check2, tempo.cat)
    arg.check2 <- c(arg.check2, TRUE)
}else{
    clone_distance <- as.numeric(clone_distance)
    tempo <- fun_check(data = clone_distance, class = "vector", prop = TRUE, neg.values = FALSE, length = 1) ; eval(ee)
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


################ Ignition


fun_report(data = paste0("\n\n################################################################ histogram PROCESS\n\n"), output = log, path = "./", overwrite = FALSE)
ini.date <- Sys.time()
ini.time <- as.numeric(ini.date) # time of process begin, converted into seconds
fun_report(data = paste0("\n\n################################ RUNNING DATE AND STARTING TIME\n\n"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0(ini.date, "\n\n"), output = log, path = "./", overwrite = FALSE)
fun_report(data = paste0("\n\n################################ RUNNING\n\n"), output = log, path = "./", overwrite = FALSE)


################ End ignition


################ Graphical parameter initialization


################ End graphical parameter initialization


################ Data import

db <- read.table(file.name, header = TRUE, sep = "\t", comment.char = "")


################ End Data import

################ Plotting




if(all(is.na(db$dist_nearest))){
    cat("\n\nNO DISTANCE HISTOGRAM PLOTTED: shazam::distToNearest() FUNCTION RETURNED ONLY NA (SEE THE dist_nearest COLUMN O THE nearest_distance.tsv FILE)")
    # no need to use pdf(NULL) with fun_gg_empty_graph()
    tempo.plot <- fun_gg_empty_graph(text = "NO DISTANCE HISTOGRAM PLOTTED\nshazam::distToNearest() FUNCTION RETURNED ONLY NA\nSEE THE dist_nearest COLUMN OF THE nearest_distance.tsv FILE", text.size = 3)
    tempo.name <- "sequence_distance" # do not use "seq_distance" !! It creates a replacement of the seq_distance.pdf file
}else{
    tempo.gg.name <- "gg.indiv.plot."
    tempo.gg.count <- 0
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::ggplot(
        data = subset(db, ! is.na(dist_nearest)),
        mapping = ggplot2::aes(x = dist_nearest), 
    ))
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::geom_histogram(color="white", binwidth=0.02))
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::xlab(paste0(clone_model, " distance")))
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::ylab("Count"))
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::geom_vline(xintercept = clone_distance, color = "firebrick", linetype = 2))
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::theme_bw())
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::theme(
        plot.margin = ggplot2::margin(t = 0.25, l = 0.1, b = 0.1, r = 0.1, unit = "in"),
        axis.line.y.right = ggplot2::element_line(color = NA), 
        axis.line.x.top = ggplot2::element_line(color = NA)
    ))
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::scale_x_continuous(breaks=seq(0, 1, 0.1)))
    tempo.title <- paste0("Manual method | Distance manually set to ", clone_distance)
    title.grob <- grid::textGrob(
        label = tempo.title,
        x = grid::unit(0, "lines"), 
        y = grid::unit(0, "lines"),
        hjust = 0,
        vjust = 0,
        gp = grid::gpar(fontsize = 6)
    )
    # end title
    pdf(NULL)
    tempo.plot <- suppressMessages(suppressWarnings(gridExtra::arrangeGrob(grobs = list(eval(parse(text = paste(paste0(tempo.gg.name, 1:tempo.gg.count), collapse = " + ")))), ncol = 1, widths = 1, top = title.grob, left = " ", right = " "))) # , left = " ", right = " " : trick to add margins in the plot. padding =  unit(0.5, "inch") is for top margin above the title
    dev.off()
    tempo.name <- "seq_distance_manual"

# Plotting of ancillary histograms
    # Automated threshold detection via smoothed density
    tempo.get <- fun_get_message("shazam::findThreshold(db$dist_nearest, method='density')", kind = "error")
    if(is.null(tempo.get)){
        output <- shazam::findThreshold(db$dist_nearest, method="density")
        save(list = ls(), file = paste0("./all.RData"))
        if( ! is.null(output)){
            png(filename = "seq_distance_smoothed_auto.png", width = 5, height = 5, units = "in", res = 300)
            plot(output, title = paste0("Density Method | threshold estimated: ", round(output@threshold, 3)))
            dev.off()
            svg(filename = "seq_distance_smoothed_auto.svg", width = 5, height = 5)
            plot(output, title = paste0("Density Method | threshold estimated: ", round(output@threshold, 3)))
            dev.off()
            pdf(file = "seq_distance_smoothed_auto.pdf", width = 5, height = 5)
            plot(output, title = paste0("Density Method | threshold estimated: ", round(output@threshold, 3)))
            dev.off()
        }else{
            # no need to use pdf(NULL) with fun_gg_empty_graph()
            tempo.plot2 <- fun_gg_empty_graph(text = "NO SMOOTHED DENSITY DISTANCE HISTOGRAM PLOTTED\nshazam::findThreshold FUNCTION RETURNED NULL", text.size = 3)
            ggplot2::ggsave(filename = "seq_distance_smoothed_auto.png", plot = tempo.plot2, device = "png", path = ".", width = 5, height = 5, units = "in", dpi = 300)
            ggplot2::ggsave(filename = "seq_distance_smoothed_auto.svg", plot = tempo.plot2, device = "svg", path = ".", width = 5, height = 5, units = "in", dpi = 300)
            ggplot2::ggsave(filename = "seq_distance_smoothed_auto.pdf", plot = tempo.plot2, device = "pdf", path = ".", width = 5, height = 5, units = "in", dpi = 300)
        }
    }else{
    # no need to use pdf(NULL) with fun_gg_empty_graph()
        tempo.plot2 <- fun_gg_empty_graph(text = "NO SMOOTHED DENSITY DISTANCE HISTOGRAM PLOTTED\nshazam::findThreshold FUNCTION RETURNED NULL", text.size = 3)
        ggplot2::ggsave(filename = "seq_distance_smoothed_auto.png", plot = tempo.plot2, device = "png", path = ".", width = 5, height = 5, units = "in", dpi = 300)
        ggplot2::ggsave(filename = "seq_distance_smoothed_auto.svg", plot = tempo.plot2, device = "svg", path = ".", width = 5, height = 5, units = "in", dpi = 300)
        ggplot2::ggsave(filename = "seq_distance_smoothed_auto.pdf", plot = tempo.plot2, device = "pdf", path = ".", width = 5, height = 5, units = "in", dpi = 300)
    }

    # Automated threshold detection via a mixture model
    tempo.get <- fun_get_message("shazam::findThreshold(db$dist_nearest, method='gmm', model='gamma-gamma'))", kind = "error")
    if(is.null(tempo.get)){
        output2 <- shazam::findThreshold(db$dist_nearest, method="gmm", model="gamma-gamma")
        save(list = ls(), file = paste0("./all2.RData"))
        if( ! is.null(output2)){
            png(filename = "seq_distance_mixture_auto.png", width = 5, height = 5, units = "in", res = 300)
            plot(output2, title = paste0("GMM Method: gamma-gamma Method | threshold estimated: ", round(output2@threshold, 3)), binwidth = 0.02)
            dev.off()
            svg(filename = "seq_distance_mixture_auto.svg", width = 5, height = 5)
            plot(output2, title = paste0("GMM Method: gamma-gamma Method | threshold estimated: ", round(output2@threshold, 3)), binwidth = 0.02)
            dev.off()
            pdf(file = "seq_distance_mixture_auto.pdf", width = 5, height = 5)
            plot(output2, title = paste0("GMM Method: gamma-gamma Method | threshold estimated: ", round(output2@threshold, 3)), binwidth = 0.02)
            dev.off()
        }else{
            # no need to use pdf(NULL) with fun_gg_empty_graph()
            tempo.plot3 <- fun_gg_empty_graph(text = "NO MIXTURE MODEL DISTANCE HISTOGRAM PLOTTED\nshazam::findThreshold FUNCTION RETURNED NULL", text.size = 3)
            ggplot2::ggsave(filename = "seq_distance_mixture_auto.png", plot = tempo.plot3, device = "png", path = ".", width = 5, height = 5, units = "in", dpi = 300)
            ggplot2::ggsave(filename = "seq_distance_mixture_auto.svg", plot = tempo.plot3, device = "svg", path = ".", width = 5, height = 5, units = "in", dpi = 300)
            ggplot2::ggsave(filename = "seq_distance_mixture_auto.pdf", plot = tempo.plot3, device = "pdf", path = ".", width = 5, height = 5, units = "in", dpi = 300)
        }
    }else{
    # no need to use pdf(NULL) with fun_gg_empty_graph()
        tempo.plot3 <- fun_gg_empty_graph(text = "NO MIXTURE MODEL DISTANCE HISTOGRAM PLOTTED\nshazam::findThreshold FUNCTION RETURNED NULL", text.size = 3)
        ggplot2::ggsave(filename = "seq_distance_mixture_auto.png", plot = tempo.plot3, device = "png", path = ".", width = 5, height = 5, units = "in", dpi = 300)
        ggplot2::ggsave(filename = "seq_distance_mixture_auto.svg", plot = tempo.plot3, device = "svg", path = ".", width = 5, height = 5, units = "in", dpi = 300)
        ggplot2::ggsave(filename = "seq_distance_mixture_auto.pdf", plot = tempo.plot3, device = "pdf", path = ".", width = 5, height = 5, units = "in", dpi = 300)
    }
    # end Plotting of ancillary histograms
}
ggplot2::ggsave(filename = paste0(tempo.name, ".png"), plot = tempo.plot, device = "png", path = ".", width = 5, height = 5, units = "in", dpi = 300)
ggplot2::ggsave(filename = paste0(tempo.name, ".svg"), plot = tempo.plot, device = "svg", path = ".", width = 5, height = 5, units = "in", dpi = 300)
ggplot2::ggsave(filename = paste0(tempo.name, ".pdf"), plot = tempo.plot, device = "pdf", path = ".", width = 5, height = 5, units = "in", dpi = 300)




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
    tempo.cat <- paste0("IN histogram.R OF THE NEXFLOW EXECUTION:\n\n", warn)
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
