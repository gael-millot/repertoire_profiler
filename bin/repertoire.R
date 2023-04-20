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



fun_gg_heatmap2 <- function(
    data1, 
    x = NULL, 
    y = NULL, 
    z, 
    label.size = 12, 
    color.low = "white", 
    color.high = "blue", 
    title = "", 
    title.size = 12
){
    # AIM
    # Plot a ggplot2 heatmap using contingency data
    # For ggplot2 specifications, see: https://ggplot2.tidyverse.org/articles/ggplot2-specs.html
    # WARNINGS
    # Size arguments (hole.text.size, border.size, title.text.size and annotation.size) are in mm. See Hadley comment in https://stackoverflow.com/questions/17311917/ggplot2-the-unit-of-size. See also http://sape.inf.usi.ch/quick-reference/ggplot2/size). Unit object are not accepted, but conversion can be used (e.g., grid::convertUnit(grid::unit(0.2, "inches"), "mm", valueOnly = TRUE))
    # ARGUMENTS
    # data1: a dataframe compatible with ggplot2 or a matrix of contingency
    # x: single character string of the data1 column name of the x-axis heatmap if data1 is a data.frame. Ignored if data1 is a matrix. Write NULL if data1 is a data.frame with a single categorical column and that the returned heatmap is vertical (single colum)
    # x: single character string of the data1 column name of the y-axis heatmap if data1 is a data.frame. Ignored if data1 is a matrix. Write NULL if data1 is a data.frame with a single categorical column and that the returned heatmap is horizontal (single row)
    # z: single character string of the data1 column name of the heatmap frequencie if data1 is a data.frame. Ignored if data1 is a matrix
    # label.size: single positive numeric value of the x-axis and y-axis font size in mm
    # color.low: a single character string or integer of the color corresponding to the lower value in the gradient color. Colors can be color names (see ?colors() in R), hexadecimal color codes, or integers (according to the ggplot2 palette)
    # color.high: as the color.lower argument but for the higher value in the gradient color
    # title: single character string of the graph title
    # title.size: single numeric value of the title font size in mm

    # add: character string allowing to add more ggplot2 features (dots, lines, themes, facet, etc.). Ignored if NULL
    # WARNING: (1) the string must start with "+", (2) the string must finish with ")" and (3) each function must be preceded by "ggplot2::". Example: "+ ggplot2::coord_flip() + ggplot2::theme_bw()"
    # If the character string contains the "ggplot2::theme" string, then the article argument of fun_gg_donut() (see above) is ignored with a warning. In addition, some arguments can be overwritten, like x.angle (check all the arguments)
    # Handle the add argument with caution since added functions can create conflicts with the preexisting internal ggplot2 functions
    # WARNING: the call of objects inside the quotes of add can lead to an error if the name of these objects are some of the fun_gg_donut() arguments. Indeed, the function will use the internal argument instead of the global environment object. Example article <- "a" in the working environment and add = '+ ggplot2::ggtitle(article)'. The risk here is to have TRUE as title. To solve this, use add = '+ ggplot2::ggtitle(get("article", envir = .GlobalEnv))'
    # return: logical (either TRUE or FALSE). Return the graph parameters?
    # return.ggplot: logical (either TRUE or FALSE). Return the ggplot object in the output list? Ignored if return argument is FALSE. WARNING: always assign the fun_gg_donut() function (e.g., a <- fun_gg_donut()) into something if the return.ggplot argument is TRUE, otherwise, double plotting is performed. See $ggplot in the RETURN section below for more details
    # return.gtable: logical (either TRUE or FALSE). Return the full graph (main, title and legend) as a gtable of grobs in the output list? See $gtable in the RETURN section below for more details
    # plot: logical (either TRUE or FALSE). Plot the graphic? If FALSE and return argument is TRUE, graphical parameters and associated warnings are provided without plotting
    # warn.print: logical (either TRUE or FALSE). Print warnings at the end of the execution? ? If FALSE, warning messages are never printed, but can still be recovered in the returned list. Some of the warning messages (those delivered by the internal ggplot2 functions) are not apparent when using the argument plot = FALSE
    # lib.path: vector of character strings indicating the absolute path of the required packages (see below). if NULL, the function will use the R library default folders
    # RETURN
    # a heatmap  if plot argument is TRUE
    # a list of the graph info if return argument is TRUE:
    # $data: the initial data with modifications and with graphic information added
    # $removed.row.nb: a list of the removed rows numbers in data frames (because of NA). NULL if no row removed
    # $removed.rows: a list of the removed rows in data frames (because of NA). NULL if no row removed
    # $plot.data
    # $panel: the variable names used for the panels (NULL if no panels). WARNING: NA can be present according to ggplot2 upgrade to v3.3.0
    # $axes: the x-axis and y-axis info
    # $warn: the warning messages. Use cat() for proper display. NULL if no warning. WARNING: warning messages delivered by the internal ggplot2 functions are not apparent when using the argument plot = FALSE
    # $ggplot: ggplot object that can be used for reprint (use print($ggplot) or update (use $ggplot + ggplot2::...). NULL if return.ggplot argument is FALSE. Warning: the legend is not in $ggplot as it is in a separated grob (use $gtable to get it). Of note, a non-null $ggplot in the output list is sometimes annoying as the manipulation of this list prints the plot
    # $gtable: gtable object that can be used for reprint (use gridExtra::grid.arrange(...$ggplot) or with additionnal grobs (see the grob decomposition in the examples). Contrary to $ggplot, a non-NULL $gtable in the output list is not annoying as the manipulation of this list does not print the plot
    # REQUIRED PACKAGES
    # ggplot2
    # gridExtra
    # grid
    # lemon (in case of use in the add argument)
    # REQUIRED FUNCTIONS FROM THE cute PACKAGE
    # fun_gg_palette()
    # fun_gg_get_legend()
    # fun_pack()
    # fun_check()
    # EXAMPLES
    # obs1 <- data.frame(X = "A", Var1 = c("TUUT", "WIIM", "BIP", "WROUM"), Count = c(20,15,0,1), stringsAsFactors = TRUE) ; fun_gg_heatmap(data1 = obs1, x = "X", y = "Var1", z = "Count", label.size = 12, size.min = 1, color.low = "white", color.high = "blue", title = tempo.title, title.size = 12)
    # DEBUGGING
    # obs1 <- data.frame(X = "A", Var1 = c("TUUT", "WIIM", "BIP", "WROUM"), Count = c(20,15,0,1), stringsAsFactors = TRUE) ; data1 = obs1 ; x = "X" ; y = "Var1" ; z = "Count" ; label.size = 12 ; size.min = 1 ; color.low = "white" ; color.high = "blue" ; title = tempo.title ; title.size = 12 ; add = NULL ; return = TRUE ; return.ggplot = FALSE ; return.gtable = TRUE ; plot = TRUE ; warn.print = FALSE ; lib.path = NULL
    # function name
    tempo.gg.name <- "gg.indiv.plot."
    tempo.gg.count <- 0


    if(class(data1) == "data.frame" & (is.null(x) | is.null(y)) & ncol(data1) != 2){
        stop(paste0("\n\n================\n\nERROR IN repertoire.R:\ndata1 ARGUMENT MUST BE A TWO COLUMN DATA FRAME IF data1 ARGUMENT IS A DATA FRAME AND IF x OR y ARGUMENT IS NULL\n\n================\n\n"), call. = FALSE)
    }
    if(class(data1) == "data.frame" & is.null(x) & ncol(data1) == 2){
        data1 <- data.frame(data1, X = "A")
        assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::ggplot(
            data = data1,
            ggplot2::aes_string(x = "X", y = y, fill= z) # Var1 with first capital letter because converted by table()
        ))
    }else if(class(data1) == "data.frame" & is.null(y) & ncol(data1) == 2){
        data1 <- data.frame(data1, Y = "A")
        assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::ggplot(
            data = data1,
            ggplot2::aes_string(x = x, y = "Y", fill= z) # Var1 with first capital letter because converted by table()
        ))
    }else if(class(data1) == "data.frame" & ncol(data1) == 3){
        assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::ggplot(
            data = data1,
             ggplot2::aes_string(x = x, y = y, fill= z) # Var1 with first capital letter because converted by table()
        ))
    }else{
        stop(paste0("\n\n================\n\nINTERNAL CODE ERROR IN repertoire.R\n\n================\n\n"), call. = FALSE)
    }
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::geom_tile())
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::scale_fill_gradient2(
        low = color.low, 
        high = color.high,
        breaks = seq(0, max(data1[ , z], na.rm = TRUE), length.out = 5),
        limits = c(0, max(data1[ , z], na.rm = TRUE))
    ))
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::coord_fixed(ratio = 1))
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::scale_y_discrete(
        expand = c(0, 0)
    ))
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::scale_x_discrete(
        expand = c(0, 0)
    ))
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::theme(
        axis.text.x =  if(is.null(x)){ggplot2::element_blank()}else{ggplot2::element_text(size = label.size, angle = 90)},
        axis.title.x =  ggplot2::element_blank(),
        axis.ticks.x = if(is.null(x)){ggplot2::element_line(size = NA)}else{ggplot2::element_line(size = 0.1)},
        axis.text.y = if(is.null(y)){ggplot2::element_blank()}else{ggplot2::element_text(size = label.size)},
        axis.title.y =  ggplot2::element_blank(),
        axis.ticks.y = if(is.null(y)){ggplot2::element_line(size = NA)}else{ggplot2::element_line(size = 0.1)},
        panel.border = element_rect(linetype = "solid", fill = NA, size = 0.1),
        panel.background = element_rect(fill = "white", colour = "white"),
        legend.title = ggplot2::element_text(size = 12), 
        legend.text = ggplot2::element_text(size = 8)
    ))
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::guides(
        fill = ggplot2::guide_colourbar(
            ticks = FALSE
        )
    ))
    bef.final.plot <- ggplot2::ggplot_build(eval(parse(text = paste(paste0(tempo.gg.name, 1:tempo.gg.count), collapse = " + "))))
    legend.final <- fun_gg_get_legend(ggplot_built = bef.final.plot) # get legend
    assign(paste0(tempo.gg.name, tempo.gg.count <- tempo.gg.count + 1), ggplot2::guides(fill = "none", color = "none", size = "none")) # inactivate the initial legend
    title.grob <- grid::textGrob(
        label = title,
        x = grid::unit(0, "lines"), 
        y = grid::unit(0, "lines"),
        hjust = 0,
        vjust = 0,
        gp = grid::gpar(fontsize = 12)
    )
    # end title
    pdf(NULL)
    output <- suppressMessages(suppressWarnings(gridExtra::arrangeGrob(grobs = list(eval(parse(text = paste(paste0(tempo.gg.name, 1:tempo.gg.count), collapse = " + "))), legend.final), ncol=2, widths=c(1, 0.5), top = title.grob, left = " ", right = " "))) # , left = " ", right = " " : trick to add margins in the plot. padding =  unit(0.5, "inch") is for top margin above the title
    suppressMessages(dev.off())
    return(output)
}

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


################ data modification, plotting and saving


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
        for(i4 in c("non-zero", "all")){
            if(sum(tempo.table, na.rm = TRUE) > 0){
                if(i4 == "non-zero"){
                    tempo.table.gg <- tempo.table[tempo.table > 0]
                }else{
                    tempo.table.gg <- tempo.table
                }
                tempo.table.gg <- as.data.frame(tempo.table.gg)
                names(tempo.table.gg)[names(tempo.table.gg) == "Freq"] <- "Count"
                output.name <- names(alleles)[tempo.pos]
                tempo.title <- paste0(
                    "Locus: ", names(alleles)[tempo.pos], "\n",
                    "Alleles: ", i4
                )
                label.size <- -5/132 * nrow(tempo.table.gg) + 1081/66 # use https://www.wolframalpha.com/widgets/view.jsp?id=f995c9aeb1565edd78adb37d2993d66
                final.plot <- fun_gg_heatmap2(
                    data1 = tempo.table.gg,
                    x = NULL,
                    y = "Var1",
                    z = "Count",
                    label.size = ifelse(label.size <= 0, 1, label.size),
                    color.low = "white",
                    color.high = "blue",
                    title = tempo.title,
                    title.size = 12
                )
            }else{
                # no need to use pdf(NULL) with fun_gg_empty_graph()
                final.plot <- fun_gg_empty_graph(text = "NO GRAPH PLOTTED FOR ", output.name, "\nNO ALLELE DETECTED", text.size = 3)
            }
            ggplot2::ggsave(filename = paste0(output.name, "_", i4, ".png"), plot = final.plot, device = "png", path = ".", width = 4, height = 10, units = "in", dpi = 300) # do not modify width and height. Otherwise impair axis.text.y, axis.ticks.y, panel.border sizes
            ggplot2::ggsave(filename = paste0(output.name, "_", i4, ".svg"), plot = final.plot, device = "svg", path = ".", width = 4, height = 10, units = "in", dpi = 300)
            ggplot2::ggsave(filename = paste0(output.name, "_", i4, ".pdf"), plot = final.plot, device = "pdf", path = ".", width = 4, height = 10, units = "in", dpi = 300)
        }
        # end plot
    }
}


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
            tempo.table2 <- table(df[names(df) %in% c(var1[i0], var1[i1])])
            output.name <- paste0(names(alleles)[tempo.pos1], "_x_", names(alleles)[tempo.pos2])
            write.table(tempo.table2, file = paste0("./rep_", output.name, ".tsv"), row.names = TRUE, col.names = NA, sep = "\t") # separate repertoires
            # plot
            for(i4 in c("non-zero", "all")){
                if(sum(tempo.table2, na.rm = TRUE) > 0){
                    if(i4 == "non-zero"){
                        tempo.log <- apply(tempo.table2, 1, sum, na.rm = TRUE) > 0
                        tempo.table3 <- tempo.table2[tempo.log, ]
                        tempo.log <- apply(tempo.table3, 2, sum, na.rm = TRUE) > 0
                        tempo.table3 <- tempo.table3[ , tempo.log]
                        tempo.table.gg <- tempo.table3
                    }else{
                        tempo.table.gg <- tempo.table2
                    }
                    tempo.table.gg <- as.data.frame(tempo.table.gg)
                    names(tempo.table.gg)[names(tempo.table.gg) == "Freq"] <- "Count"
                    tempo.title <- paste0(
                        "Locus: ", output.name, "\n",
                        "Alleles: ", i4
                    )
                    label.size <- -5/132 * nrow(tempo.table.gg) + 1081/66 # use https://www.wolframalpha.com/widgets/view.jsp?id=f995c9aeb1565edd78adb37d2993d66
                    final.plot <- fun_gg_heatmap2(
                        data1 = tempo.table.gg,
                        x = var1[2],
                        y = var1[1],
                        z = "Count",
                        label.size = ifelse(label.size <= 0, 1, label.size), 
                        color.low = "white",
                        color.high = "blue",
                        title = tempo.title,
                        title.size = 12
                    )
                }else{
                    # no need to use pdf(NULL) with fun_gg_empty_graph()
                    final.plot <- fun_gg_empty_graph(text = "NO GRAPH PLOTTED FOR ", output.name, "\nNO ALLELE DETECTED", text.size = 3)
                }
                ggplot2::ggsave(filename = paste0(output.name, "_", i4, ".png"), plot = final.plot, device = "png", path = ".", width = 4, height = 10, units = "in", dpi = 300) # do not modify width and height. Otherwise impair axis.text.y, axis.ticks.y, panel.border sizes
                ggplot2::ggsave(filename = paste0(output.name, "_", i4, ".svg"), plot = final.plot, device = "svg", path = ".", width = 4, height = 10, units = "in", dpi = 300)
                ggplot2::ggsave(filename = paste0(output.name, "_", i4, ".pdf"), plot = final.plot, device = "pdf", path = ".", width = 4, height = 10, units = "in", dpi = 300)
            }
            # end plot
        }
    }
}

tempo <- qpdf::pdf_combine(input = list.files(path = ".", pattern = ".pdf$"), output = "./repertoire.pdf") # assignation to prevent a returned element


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
