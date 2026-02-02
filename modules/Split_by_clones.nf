

process Split_by_clones { // split the file into multiple files according to the clone_id column
    label 'immcantation'
    cache 'true'

    input:
    tuple path(clone_ch), val(nlines) // no parallelization

    when:
    nlines > 1 // over 1 means "more than header only"

    output:
    path "*_wanted_seq_clone-pass.tsv", emit: clone_split_ch // multiple files -> parall expected

    script:
    """
    #!/bin/bash -ue
    set -o pipefail
    FILENAME=\$(basename -- ${clone_ch}) # recover a file name without path
    cp -Lr ${clone_ch} "./TEMPO.tsv" # to have the hard file, not the symlink, because modifications will be performed inside
    chmod 777 TEMPO.tsv
    rm \$FILENAME # remove the initial file to avoid to send it into the channel
    cp -rp TEMPO.tsv "\$FILENAME" # -p for preserve permissions
    rm TEMPO.tsv

    # Determine the position of the clone_id column
    clone_id_col=\$(Rscript -e '
        args <- commandArgs(trailingOnly = TRUE)
        filename <- args[1]

        # Read the header of the file
        header <- read.table(filename, sep = "\\t", header = TRUE, nrows = 1)

        # Find the position of the clone_id column
        clone_id_col <- which(names(header) == "clone_id")

        if (length(clone_id_col) == 0) {
          cat("ERROR: clone_id column not found in the input file\\n")
          quit(status = 1)
        } else {
          cat(clone_id_col, "\\n")
        }
    ' \$FILENAME)

    if [ -z "\$clone_id_col" ]; then
        echo "ERROR IN THE SPLIT_BY_CLONES PROCESS: clone_id column not found in the input file"
        exit 1
    fi

    awk -v var1=\$FILENAME -v col=\$clone_id_col -F "\\t" '{
        print \$col
        if(NR == 1){
            header=\$0 ; next
        }else{
            clone_id=\$col
            if(system( "[ -f " \$col"_"var1 " ] " ) > 0){ # test if a file exists
                print header > \$col"_"var1
            }
            print \$0 > \$col"_"var1
        }
    }' \$FILENAME
    """
}

