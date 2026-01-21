process Clone_id_count {
    label 'immcantation'
    publishDir path: "${out_path}/tsv", mode: 'copy', pattern: "clone_id_count.tsv", overwrite: false
    publishDir path: "${out_path}/reports", mode: 'copy', pattern: "clone_id_count.log", overwrite: false
    cache 'true'

    input:
    path clone_assigned_seq_ch // no parallelization


    output:
    path "clone_id_count.tsv"
    path "clone_id_count.log"

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    Rscript -e '
        obs <- read.table("./clone_assigned_seq.tsv", sep = "\\t", header = TRUE)
        clone_id_count <- as.data.frame(table(obs\$clone_id))
        names(clone_id_count) <- c("clone_id", "count")
        write.table(clone_id_count, file = "clone_id_count.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\\t")
    ' |& tee -a clone_id_count.log
    """
}
