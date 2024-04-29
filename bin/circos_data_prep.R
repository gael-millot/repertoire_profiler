erase.objects = TRUE # write TRUE to erase all the existing objects in R before starting the algorithm and FALSE otherwise. Beginners should use TRUE
if(erase.objects == TRUE){
    rm(list = ls(all.names = TRUE))
    erase.objects = TRUE
}
erase.graphs = TRUE # write TRUE to erase all the graphic windows in R before starting the algorithm and FALSE otherwise
script <- "circos_data_prep"

VH.path <- "X:/ROCURONIUM PROJECT/01 Primary data/04.Repertoire analysis/SORT1/SORT1 Seq-original/RESULT/repertoire_profiler_1708127835/all_passed_seq.tsv" # Must be VH
VL.path <- "X:/ROCURONIUM PROJECT/01 Primary data/04.Repertoire analysis/SORT1/SORT1 Seq-original/RESULT/repertoire_profiler_1708128192/all_passed_seq.tsv" # must be VL
merging_colums_names <- c("initial_sequence_id") # single vector of character strings of the names of the column used to merge the 2 files. Of note, if _VH or _VL in the columns names, they are removed.
metadata_colums_names <- c("sequence_id") # single character string of the name of the metadata column in the two files. Write NULL if no metadata.
kept_colums_names <- c("germline_v_call",	"germline_j_call") # single vector of character strings of the names of the column of VH.path file that correspond to V (first position) and J. Of note, Same names must be for VH and VL.



df1 <- read.table(VH.path, sep = "\t", header = TRUE)
df2 <- read.table(VL.path, sep = "\t", header = TRUE)
if( ! all(merging_colums_names %in% names(df1))){
    stop("ERROR 1")
}
if( ! all(merging_colums_names %in% names(df2))){
    stop("ERROR 2")
}
if( ! all(kept_colums_names %in% names(df1))){
    stop("ERROR 3")
}
if( ! all(kept_colums_names %in% names(df2))){
    stop("ERROR 4")
}
if(any(duplicated(df1[ , merging_colums_names]))){
    stop("ERROR 5")
}
if(any(duplicated(df2[ , merging_colums_names]))){
    stop("ERROR 6")
}

# removal of _VH and _VL in these columns
df1[ , merging_colums_names] <- gsub(df1[ , merging_colums_names], pattern = "_VH", replacement = "")
df2[ , merging_colums_names] <- gsub(df2[ , merging_colums_names], pattern = "_VL", replacement = "")
df1[ , metadata_colums_names] <- gsub(df1[ , metadata_colums_names], pattern = "_VH", replacement = "")
df2[ , metadata_colums_names] <- gsub(df2[ , metadata_colums_names], pattern = "_VL", replacement = "")
# end removal of _VH and _VL in these columns

# columns rename before merging
renamed_colums_names <- c("V_allele", "J_allele") # single vector of character strings of the names of VH.path file to keep, beyond metadata. Of note, Same names must be for VH and VL.
names(df1)[match(kept_colums_names, names(df1))] <- paste0(renamed_colums_names, "_VH")
names(df1)[match(metadata_colums_names, names(df1))] <- "Metadata_VH"
names(df2)[match(kept_colums_names, names(df2))] <- paste0(renamed_colums_names, "_VL")
names(df2)[match(metadata_colums_names, names(df2))] <- "Metadata_VL"
# end columns rename before merging

# selection of the columns of interest
df1 <- df1[ , c(merging_colums_names,"Metadata_VH", paste0(renamed_colums_names, "_VH"))]
df2 <- df2[ , c(merging_colums_names, "Metadata_VL", paste0(renamed_colums_names, "_VL"))]
# end selection of the columns of interest

# merging
tempo.log1 <- is.na(match(df1[ , merging_colums_names], df2[ , merging_colums_names]))
tempo.log2 <- is.na(match(df2[ , merging_colums_names], df1[ , merging_colums_names]))
if(any(tempo.log1)){
    cat(paste0("\n\nWarning: SOME SEQUENCES FROM VH ARE REMOVED BECAUSE NOT IN THE OTHER VL FILE:\n", paste(df1[ , merging_colums_names][tempo.log1], collapse = "\n"), "\n\n"))
}
if(any(tempo.log2)){
    cat(paste0("\n\nWarning: SOME SEQUENCES FROM VL ARE REMOVED BECAUSE NOT IN THE OTHER VH FILE:\n", paste(df2[ , merging_colums_names][tempo.log2], collapse = "\n"), "\n\n"))
}
df1 <- df1[ ! tempo.log1, ]
df2 <- df2[ ! tempo.log2, ]
nrow(df1)
nrow(df2)
identical(sort(df1[ , merging_colums_names]), sort(df2[ , merging_colums_names]))
output_all_seq <- merge(df1, df2, by = c("initial_sequence_id"), sort = FALSE) # rd for random. Send the coord of the boxes into the coord data.frame of the dots (in the column x.y). WARNING: by = c("PANEL", "group") without fill column because PANEL & group columns are enough as only one value of x column per group number in box.coord. Thus, no need to consider fill column
nrow(output_all_seq)
if(nrow(output_all_seq) != nrow(df1)){
    tempo.cat <- paste0("INTERNAL CODE ERROR. THE merge() FUNCTION DID NOT RETURN A CORRECT output DATA FRAME. CODE HAS TO BE MODIFIED")
    stop(paste0("\n\n================\n\n", tempo.cat, "\n\n================\n\n", call. = FALSE)) # == in stop() to be able to add several messages between ==
}
# end merging

# new columns wo * in alleles
output_all_seq <- data.frame(
    output_all_seq, 
    V_VH = (output_all_seq$V_allele_VH)
)
# end new columns wo * in alleles


# new fusion columns
output_all_seq <- data.frame(
    output_all_seq, )
# end new fusion columns

# metadata in output_all_seq
if( ! all(output_all_seq$Metadata_VH == output_all_seq$Metadata_VL)){
    stop("ERROR 7")
}
output_all_seq <- output_all_seq[ , names(output_all_seq) != "Metadata_VL"] # removal of Metadata_VL
names(output_all_seq)[names(output_all_seq) == "Metadata_VH"] <- "Metadata" # remane Metadata_VH into Metadata
tempo.log <- is.na(match(output_all_seq[ , merging_colums_names], output_all_seq$Metadata)) # detection of the Metadata names
output_all_seq$Metadata[ ! tempo.log] <- NA # keep only the Metadata names in the Metadata column
# end metadata in output_all_seq

# export output_all_seq

# end export output_all_seq
