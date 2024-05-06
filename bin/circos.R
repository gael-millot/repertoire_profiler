tsv.path <- "C:\\Users\\gmillot\\Desktop\\mouse2_circos_data.tsv" # "C:\\Users\\gmillot\\Documents\\Hub projects\\20210914 Alice Dejoux 19463\\dataset\\24_03_18_circo plot VH-VL appariesv_meta_data.xlsx"
out.path <- "C:\\Users\\gmillot\\Desktop"
pdf.name <- "20240502_unpaired_HL_mouse2"
Name <- "Name" # column of the sequence names
col1 <- "VJ_VH" # first column name
col2 <- "VJ_VL" # second column name
col1.selection <- c("IGHV1-55_IGHJ2", "IGHV3-6_IGHJ3") # single vector of characters indicating the values to select inside col1 for the circos representation. Write NULL if not required
col2.selection <- c("IGKV4-68_IGKJ5", "IGKV12-89_IGKJ4") # single vector of characters indicating the values to select inside col2 for the circos representation. Write NULL if not required
coupled.selection <- FALSE # single logical value. Couple col1.selection and col2.selection? If TRUE, col1.selection and col2.selection must be of same length: element 1 of col1.selection is parired with 1 col2.selection, etc. If FALSE, select the values without pairing. paired signifie que seuls les Ig qui ont le même couple H et L que les fonctionnels ont été gardés. Cela veut dire que si IP1 est H1 pour H et L1 pour L, j'ai récup tous les Ig H1-L1. Unpaired signifie que seuls les Ig qui ont un H des H fonctionnels avec l'un des L fonctionnels ont été gardés. Cela signifie que si IP1 est H1 pour H et L1 pour L, et IP2 est H2 et L2, j'ai récup tous les IG H1-L1, H1-L2, H2-L1 et H2, L2.
metadata <- c("mAb_name", "Kd_nM") # single vector of characters indicating the name of the column with metadata
scale.size <- 0.7 # size of the numbers
label.size <- 0.7 # size of the label text
label.adj <- c(-0.1, 0.5) # label text adjustment
palette.name <- "fun_gg_palette" # palette of color. See http://127.0.0.1:22643/library/grDevices/html/palettes.html
kind <- "std" # for gg_palette only. Either "std" for standard gg colors, "dark" for darkened gg colors, or "light" for pastel gg colors
alpha <- 0.5 # transparency
margins <- 1.5 # in inches
windows.size <- 7 #in inches


#output:
# note message "Note: 1 point is out of plotting region in sector 'IGHV1-18_IGHJ2', track '1'". Does not compromise the plotting in fact. See https://github.com/jokergoo/circlize/blob/master/R/global.R. "points.overflow.warning" Since each cell is in fact not a real plotting region but only an ordinary rectangle, it does not eliminate points that are plotted out of the region. So if some points are out of the plotting region, circlize would continue drawing the points and printing warnings. In some cases, draw something out of the plotting region is useful, such as draw some legend or text. Set this value to FALSE to turn off the warnings. By canvas.xlim and canvas.ylim(described in the same page above), you can set up your canvas in order to avoid, or just ignore the warning.



 suppressPackageStartupMessages(library(saferGraph))
 suppressPackageStartupMessages(library(circlize))
 suppressPackageStartupMessages(library(rJava)) # dependency of xlsx. Warning: requires the install of java 64.bit, not 32, using https://www.java.com/fr/download/manual.jsp. Otherwise, error message: le chargement du package ou de l'espace de noms a échoué pour ‘xlsx’ :.onLoad a échoué dans loadNamespace() pour 'rJava', détails
 suppressPackageStartupMessages(library(xlsx))
set.seed(999)


fun_gg_palette <- function(n, kind = kind, alpha = alpha){
# AIM
# provide colors used by ggplot2
# the interest is to use another single color that is not the red one used by default
# ARGUMENTS
# n: number of groups on the graph
# kind: either "std" for standard gg colors, "dark" for darkened gg colors, or "light" for pastel gg colors
# RETURN
# the vector of hexadecimal colors
# EXAMPLES
# the ggplot2 palette when asking for 7 different colors
# plot(1:7, pch = 16, cex = 5, col = fun_gg_palette(n = 7))
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = if(kind == "std"){65}else if(kind == "dark"){35}else if(kind == "light"){85}, c = 100, alpha = alpha)[1:n]
}


open2(pdf = TRUE, pdf.path = out.path, pdf.name = pdf.name, width = windows.size, height = windows.size)
prior_plot(
    param.reinitial = TRUE,
    xlog.scale = FALSE,
    ylog.scale = FALSE,
    remove.label = TRUE,
    remove.x.axis = TRUE,
    remove.y.axis = TRUE,
    std.x.range = TRUE,
    std.y.range = TRUE,
    down.space = margins,
    left.space = margins,
    up.space = margins,
    right.space = margins,
    orient = 1,
    dist.legend = 3.5,
    tick.length = 0.5,
    box.type = "n",
    amplif.label = 1,
    amplif.axis = 1,
    display.extend = TRUE,
    return.par = FALSE
)


obs1 <- read.table(tsv.path, header = TRUE, sep = "\t")
if(coupled.selection == FALSE){
    if(( ! is.null(col1.selection)) & ! is.null(col2.selection)){
        if( ! all(col1.selection %in% obs1[ , col1]) ){
            stop("ERROR 4")
        }
        if( ! all(col2.selection %in% obs1[ , col2]) ){
            stop("ERROR 5")
        }
        obs1 <- obs1[obs1[ , col1] %in%  col1.selection & obs1[ , col2] %in%  col2.selection, ]
    }else{
        if( ! is.null(col1.selection)){
            if( ! all(col1.selection %in% obs1[ , col1]) ){
                stop("ERROR 1")
            }else{
                obs1 <- obs1[obs1[ , col1] %in%  col1.selection, ]
            }
        }else if( ! is.null(col2.selection)){
            if( ! all(col2.selection %in% obs1[ , col2]) ){
                stop("ERROR 2")
            }else{
                obs1 <- obs1[obs1[ , col2] %in%  col2.selection, ]
            }
        }
    }
}else if(coupled.selection == TRUE){
    if(( ! is.null(col1.selection)) & ! is.null(col2.selection)){
        if(length(col1.selection) != length(col2.selection)){
            stop("ERROR 3")
        }
        if( ! all(col1.selection %in% obs1[ , col1]) ){
            stop("ERROR 4")
        }
        if( ! all(col2.selection %in% obs1[ , col2]) ){
            stop("ERROR 5")
        }
        tempo_df <- data.frame()
        for(i2 in 1:length(col1.selection)){
                    tempo_df <- rbind(tempo_df, obs1[obs1[ , col1] ==  col1.selection[i2] & obs1[ , col2] ==  col2.selection[i2], ])
        }
        obs1 <- tempo_df
    }else{
        stop("ERROR 6")
    }
}
obs2 <- table(obs1[ , c(col1, col2)])
obs3 <- as.data.frame(obs2)
if(palette.name == "fun_gg_palette"){
    grid.col <- setNames(
        eval(parse(text = paste0(palette.name, "(n = ", length(unlist(dimnames(obs2))), ", kind = \"", kind, "\", alpha = ", alpha, ")"))), 
        union(rownames(obs2), colnames(obs2))
    )
}else{
    grid.col <- setNames(
        eval(parse(text = paste0(palette.name, "(n = ", length(unlist(dimnames(obs2))), ", alpha = ", alpha, ")"))), 
        union(rownames(obs2), colnames(obs2))
    )
}
# https://jokergoo.github.io/circlize_book/book/the-chorddiagram-function.html
# http://127.0.0.1:22643/library/circlize/html/chordDiagram.html
circos.par(start.degree = 270) # H on left and K on right
chordDiagram( # circos plot of a table
    obs2, # table. The circos connect rows versus columns with stripes with width depending on the values inside the table. Cann also be a matrix or a data frame (pass all argument to chordDiagramFromMatrix or chordDiagramFromDataFrame). If it is in the form of a matrix, it should be an adjacency matrix. If it is in the form of a data frame, it should be an adjacency list
    annotationTrack = "grid", # name track and axis over track removed. Only the track "color of sector". Thus only 1 track here
    # preAllocateTracks = 1, # a new track added. Thus a total of 2 tracks
    grid.col = grid.col, 
    link.lwd = 1,
    link.lty = 1,
    link.border = "black",
    directional = 2
)

# https://jokergoo.github.io/circlize_book/book/circular-layout.html#panel-fun
circos.trackPlotRegion( # deal with features inside a track
    track.index = 1, # track number, here the track "color of sector"
    panel.fun = function(x, y) { # function for each sector
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        sector.name = get.cell.meta.data("sector.index")
        # Draw sector names
        circos.text(mean(xlim), ylim[1] + 2, sector.name, facing = "clockwise", niceFacing = TRUE, adj = label.adj, cex = label.size) # message because ylim[1] + 2 is to big
        # Annotate with metadata

        # Draw the axis
        circos.axis(h = "top", labels.cex = scale.size, major.tick.length = 0.2, sector.index = sector.name) # add axis over track 1. track.index = 1 removed because already on track 1. But can be used to specified another track while we have in the panel.fun of  track.index = 1 
    }, 
    bg.border = NA
)
circos.clear()
close2()
cat(paste0("\n\nNumber of H-L sequences: ", nrow(obs1)), "\n\n")

