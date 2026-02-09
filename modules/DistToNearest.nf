
// Compares cdr3 différences from sequences with same cdr3 length, same v and same j
// Number of substitutions normalized by the length of the cdr3
// NA means cdr3 sequences are exactly the same for the group
process DistToNearest {
    label 'immcantation'
    publishDir path: "${out_path}/tsv", mode: 'copy', pattern: "{dist_ignored.tsv}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{distToNearest.log}", overwrite: false
    cache 'true'

    input:
    path trimtranslate_ch2 // no parallelization
    val clone_model
    val clone_normalize

    output:
    path "nearest_distance.tsv", emit: distToNearest_ch
    path "dist_ignored.tsv", emit: dist_ignored_ch
    path "distToNearest.log"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    Rscript -e '
         # WEIRD stuf: if db alone is returned, and if distToNearest_ch is used for the clone_assignment process and followings, everything is fine. But if db3 is returned with db3 <- data.frame(db, dist_nearest = db2\$dist_nearest) or db3 <- data.frame(db, d = db2\$dist_nearest) or data.frame(db, d = db\$sequence_id) or db3 <- data.frame(db, d = as.numeric(db2\$dist_nearest)) or db3 <- data.frame(db[1:3], d = db\$sequence_id, db[4:length(db)]), the get_germ_tree process cannot make trees, while the wanted_seq.tsv seem identical at the end, between the use of db or db3, except that the clone_id order is not the same
        db <- read.table("${trimtranslate_ch2}", header = TRUE, sep = "\\t")
        if(all(is.na(db\$junction)) || all(is.na(db\$locus))){
            not_selected_df <- db[ -(1:nrow(db)), ]
            db2 <- data.frame(db, dist_nearest = NA)
        }else{
            db2 <- shazam::distToNearest(db, sequenceColumn = "junction", locusColumn = "locus", model = "${clone_model}", normalize = "${clone_normalize}", nproc = 1)
            missing_line_nb <- NULL
            tempo_log <- is.na(match(db\$sequence_id, db2\$sequence_id))
            if(any(tempo_log)){
                missing_line_nb <- which(tempo_log)
            }
            if(is.null(missing_line_nb)){
                not_selected_df <- db[ -(1:nrow(db)), ]
            }else{
                not_selected_df <- db[missing_line_nb, ]
                tempo_db <- as.data.frame(matrix(data = NA, nrow = nrow(not_selected_df), ncol = ncol(db2)))
                names(tempo_db) <- names(db2)
                tempo_db\$sequence_id <- db\$sequence_id[missing_line_nb]
                db2 <- rbind(db2, tempo_db)
            }
        }
        write.table(not_selected_df, file = paste0("./dist_ignored.tsv"), row.names = FALSE, col.names = TRUE, sep = "\\t", quote = FALSE)
        write.table(db2, file = paste0("./nearest_distance.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\\t")
    ' |& tee -a distToNearest.log
    """
}
