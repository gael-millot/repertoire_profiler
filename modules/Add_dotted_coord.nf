// add corrdinates for IMGT dotted alignment sequences
process Add_dotted_coord {
    label 'r_ext'
    cache 'true'

    input:
    path trimtranslate_ch // parallelization expected

    output:
    path "add_dotted_coord.tsv", emit: add_dotted_coord_ch // productive file with column sequence_alignment_aa added
    path "add_dotted_coord.log", emit: add_dotted_coord_log_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    Rscript -e '
        df <- read.table("${trimtranslate_ch}", sep = "\\t", header = TRUE)
        if( ! "sequence_alignment_with_gaps" %in% names(df)){
            stop(paste0("\\n\\n============\\n\\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\\n\\nsequence_alignment_with_gaps IS NOT A COLUMN NAME OF trimtranslate.tsv\\n\\n============\\n\\n"), call. = FALSE)
        }
        map_gapped <- function(seq, ungapped_pos) {
        # AIM
        # shift vdjc or fwr/cdr position provided in the .tsv file, according to the dots (...) inserted in the aligned sequences, to have the correct positions of the features in the aligned sequences in the html file.
        # ARGUMENTS
        # seq: aligned sequence (with dots)
        # ungapped_pos: initial positions in the unaligned sequence (without hyphens)
        # RETURN
        # shifted position
        # EXAMPLES
        # map_gapped("atg..t...g", c(3,4,5)) # initial sequence is atgtg, input coords are third g, fourth t and fifth g
        # DEBUGGING
        # seq = "atg..t...g" ; ungapped_pos = c(3,4,5)
            # Remove gaps but keep index positions of non-gaps
            aligned_chars <- strsplit(seq, "")[[1]]
            non_gap_indices <- which( ! (aligned_chars == "." | aligned_chars == "-"))
            dots_nb <- sum(aligned_chars == "." | aligned_chars == "-", na.rm = TRUE)
            # Return the gapped index corresponding to the ungapped position
            if(length(non_gap_indices) < ungapped_pos){
                return(ungapped_pos + dots_nb)
            }else{
                return(non_gap_indices[ungapped_pos])
            }
        }
        features <- c("v", "d", "j", "c", "fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")
        tempo_seq <- sub(pattern = "-+\$", replacement = "", x = df\$sequence_alignment_with_gaps) # removal of hyphens at the end of the seq only # returns NA if df\$sequence_alignment_with_gaps is NA
        for(i2 in features){
            start_coord <- df[ , paste0(i2, "_alignment_start")]
            end_coord <- df[ , paste0(i2, "_alignment_end")]
            # shifed coordinates due to hyphens in the aligned seq
            if(( ! is.na(start_coord)) & ! is.na(tempo_seq)){
                start_coord <- map_gapped(seq = tempo_seq, ungapped_pos = start_coord)
            }else{
                start_coord <- NA
            }
            if(( ! is.na(end_coord)) & ! is.na(tempo_seq)){
                end_coord <- map_gapped(seq = tempo_seq, ungapped_pos = end_coord)
            }else{
                end_coord <- NA
            }
            df <- data.frame(df, TEMPO_NAME_START = start_coord, TEMPO_NAME_END = end_coord)
            names(df)[names(df) == "TEMPO_NAME_START"] <- paste0(i2, "_alignment_with_gaps_start")
            names(df)[names(df) == "TEMPO_NAME_END"] <- paste0(i2, "_alignment_with_gaps_end")
        }
        write.table(df, file = "./add_dotted_coord.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\\t")
    ' |& tee -a add_dotted_coord.log
    """
}