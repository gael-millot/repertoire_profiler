
// Add the dist_to_nearest values and gene columns in the seq_name_replacement file, thereby creating the wanted_seq.tsv file
// The gene columns are created with the allele columns, minus the allele part
// + isotype class (4 first characters of the c_gene column)
process data_assembly {
    label 'immcantation'
    publishDir path: "${out_path}/tsv", mode: 'copy', pattern: "{wanted_seq.tsv}", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "{data_assembly.log}", overwrite: false
    cache 'true'

    input:
    path seq_name_replacement_ch2 // no parallelization
    path distToNearest_ch

    output:
    path "wanted_seq.tsv", emit: wanted_ch
    path "data_assembly.log"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    Rscript -e '
        db <- read.table("${seq_name_replacement_ch2}", header = TRUE, sep = "\\t")
        dtn <- read.table("${distToNearest_ch}", header = TRUE, sep = "\\t")


        # replace TRUE and FALSE of the read.table() conversion by the initial T and F
        # tempo.log <- sapply(db, FUN = function(x){class(x) == "logical"})
        # db[tempo.log] <- lapply(db[tempo.log], FUN = function(x){substr(as.character(x), 1, 1)})
        # end replace TRUE and FALSE of the read.table() conversion by the initial T and F


        if( ! any(c("sequence_id", "initial_sequence_id") %in% names(db))){
            tempo.cat <- paste0("\\n\\n========\\n\\nINTERNAL ERROR IN THE NEXTFLOW EXECUTION OF THE data_assembly PROCESS\\ndb SHOULD HAVE \\"sequence_id\\" AND ALSO POTENTIALLY \\"initial_sequence_id\\" AS COLUMN NAME. HERE:\\nNAMES: ", paste(names(db), collapse = " "), "\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n")
            stop(tempo.cat)
        }else{
            tempo.col.name <- ifelse("initial_sequence_id" %in% names(db), "initial_sequence_id", "sequence_id")
        }
        if(all(c("sequence_id", "initial_sequence_id") %in% names(db))){
            if(all(db\$sequence_id == db\$initial_sequence_id)){
                tempo.cat <- paste0("\\n\\n\\n\\nWARNING IN THE data_assembly PROCESS OF NEXTFLOW\\nTHE meta_path AND meta_name_replacement PARAMETERS ARE NOT \\"NULL\\" BUT NO SEQUENCE NAMES HAVE BEEN REPLACED WHEN USING THE meta_name_replacement COLUMN\\n\\n\\n\\n")
                cat(tempo.cat)
            }
        }
        if(any(is.na(match(db[, tempo.col.name], dtn\$sequence_id)))){
            tempo.cat <- paste0("\\n\\n========\\n\\nINTERNAL ERROR IN THE NEXTFLOW EXECUTION OF THE data_assembly PROCESS\\nNO NA SHOULD APPEAR AT THAT STAGE WITH match()\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n")
            stop(tempo.cat)
        }
        if(sum(db[, tempo.col.name] %in% dtn\$sequence_id, na.rm = TRUE) != nrow(db)){
            tempo.cat <- paste0("\\n\\n========\\n\\nINTERNAL ERROR IN THE NEXTFLOW EXECUTION OF THE data_assembly PROCESS\\nsum(db[, tempo.col.name] %in% dtn\$sequence_id, na.rm = TRUE) SHOULD BE EQUAL TO nrow(db)\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n")
            stop(tempo.cat)
        }
        if(nrow(dtn) < nrow(db)){
            tempo.cat <- paste0("\\n\\n========\\n\\nINTERNAL ERROR IN THE NEXTFLOW EXECUTION OF THE data_assembly PROCESS\\ndtn CANNOT HAVE LESS ROWS THAN db. HERE:\\ndb: ", nrow(db), "\\ndtn: ", nrow(dtn), "\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n")
            stop(tempo.cat)
        }
        # selection of rows and ordering of dtn
        tempo.dtn <- dtn[dtn\$sequence_id %in% db[, tempo.col.name], ] # nb rows of dtn reduced to the one of db
        if(nrow(tempo.dtn) != nrow(db)){
            tempo.cat <- paste0("\\n\\n========\\n\\nINTERNAL ERROR IN THE NEXTFLOW EXECUTION OF THE data_assembly PROCESS\\ntempo.dtn SHOULD HAVE THE SAME NUMBER OF ROWS AS db. HERE:\\ndb: ", nrow(db), "\\ntempo.dtn: ", nrow(tempo.dtn), "\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n")
            stop(tempo.cat)
        }
        tempo.dtn <- tempo.dtn[match(db[, tempo.col.name], tempo.dtn\$sequence_id), ]
        if( ! all(db[, tempo.col.name] == tempo.dtn\$sequence_id)){
            tempo.cat <- paste0("\\n\\n========\\n\\nINTERNAL ERROR IN THE NEXTFLOW EXECUTION OF THE data_assembly PROCESS\\n", tempo.col.name, " COLUMNS OF db AND sequence_id COLUMN OF dtn SHOULD BE IDENTICAL. HERE THEY ARE\\ndb\$", tempo.col.name, ":\\n", paste0(db[, tempo.col.name], collapse = ","), "\\ndtn\$sequence_id:\\n", paste0(dtn\$sequence_id, collapse = ","), "\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n")
            stop(tempo.cat)
        }
        db3 <- data.frame(db, dist_nearest = tempo.dtn\$dist_nearest)
        if( ! "c_call" %in% names(db)){
            tempo.cat <- paste0("\\n\\n========\\n\\nINTERNAL ERROR IN THE NEXTFLOW EXECUTION OF THE data_assembly PROCESS\\ndb SHOULD HAVE \\"c_call\\" AS COLUMN NAME. HERE:\\nNAMES: ", paste(names(db), collapse = " "), "\\n\\nPLEASE, REPORT AN ISSUE HERE https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues OR AT gael.millot<AT>pasteur.fr.\\n\\n========\\n\\n")
            stop(tempo.cat)
        }

        #### remove allele info

        tempo_v <- strsplit(db3\$v_call, ",")
        sub_v <- sapply(X = tempo_v, FUN = function(x){y <- sub(pattern = "\\\\*.*", replacement = "", x = x) ; paste0(unique(y), collapse = ",")})
        tempo_j <- strsplit(db3\$j_call, ",")
        sub_j <- sapply(X = tempo_j, FUN = function(x){y <- sub(pattern = "\\\\*.*", replacement = "", x = x) ; paste0(unique(y), collapse = ",")})
        # The c_call column must be handled with more care because it can have empty values
        subclass <- sapply(db3\$c_call, function(x) {
            if (is.na(x) || x == "") { # The c_call column can have empty values if the sequences are too short
                return(NA)  # Si NA ou "", on affecte NA
            } else {
                # Traitement normal: suppression des suffixes après *
                y <- sub(pattern = "\\\\*.*", replacement = "", x = unlist(strsplit(x, ",")))
                return(paste0(unique(y), collapse = ","))
            }
        })
        class <- sapply(X = subclass, FUN = function(x) {
            if (is.na(x)) {
                return(NA)  # Si x est NA, on garde NA
            } else {
                y <- sub(pattern = "\\\\*.*", replacement = "", x = x) 
                y <- substr(y, 1, 4)
                return(paste0(unique(y), collapse = ","))
            }
        })
        db4 <- data.frame(db3, v_gene = sub_v, j_gene = sub_j, isotype_class = class, c_gene = subclass)

        #### end remove allele info

        write.table(db4, file = paste0("./wanted_seq.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\\t")
    ' |& tee -a data_assembly.log
    """
}