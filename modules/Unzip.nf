

process Unzip {
    label 'unzip'
    cache 'true'

    input:
    path zip
    val sample_path

    output:
    path "*", emit: unzip_ch

    script:
    """
    #!/bin/bash -ue
    set -o pipefail # return the exit code of the first nonzero command in the pipeline, not just the last one (important when using tee)
    FILENAME_INI=\$(basename -- "${sample_path}")
    FILE_EXTENSION="\${FILENAME_INI##*.}"
    FILENAME="\${FILENAME_INI%.*}"
    if [[ ! "\${FILE_EXTENSION}" =~ zip ]] ; then
        echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nTHE FILE EXTENSION MUST BE \\".zip\\" AND NOT \${FILENAME_INI}\\n\\n========\\n\\n"
        exit 1
    else
        # hidden files (.DS_Store, ._* AppleDouble files of macOS, etc.) and __MACOSX folder (added by macOS in zip files) are neither listed nor extracted
        # [^/] and [!/] below: names starting by ./ (zip made by bsdtar) are kept
        SINGLE_FILE_NAME=\$(unzip -Z1 "${zip}" | awk '!/^__MACOSX\\// && !/(^|\\/)\\.[^\\/]/')
        if [[ -z "\${SINGLE_FILE_NAME}" ]] ; then
            echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nTHE ZIP FILE CONTAINS NO FILE (HIDDEN FILES AND __MACOSX FOLDER IGNORED):\\n\${FILENAME_INI}\\n\\n========\\n\\n"
            exit 1
        elif [[ "\${SINGLE_FILE_NAME}" == "\$(unzip -Z1 "${zip}")" ]] ; then # nothing to exclude
            unzip "${zip}"
        else
            unzip "${zip}" -x '__MACOSX/*' '.[!/]*' '*/.[!/]*' # a pattern that matches nothing only prints "caution: excluded filename not matched" on stderr
        fi
    fi
    rm \$FILENAME_INI # remove the zipped file 
    if [ -d "\$FILENAME" ]; then # if the unzipped file is a directory
        for file in "\$FILENAME"/* ; do # all the elements of the folder (hidden files are not extracted)
            # Check if the file is a regular file and not a directory
            if [ -f "\$file" ] ; then
                # Check if the file does not match the extensions .fasta or .fa
                if [[ ! "\$file" =~ \\.(fasta|fa|fas|fna|txt|seq|faa)\$ ]] ; then
                    echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nALL THE UNZIPPED FILE EXTENSIONS MUST BE fasta, fa, fas, fna, txt, seq OR faa.\\nHERE A FILE IS:\\n\$file.\\n\\n========\\n\\n"
                    exit 1
                fi
            else
                echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nTHE UNZIPPED FOLDER MUST CONTAIN ONLY FILES, WITHOUT SUBFOLDER.\\nHERE AN ELEMENT IS NOT A FILE:\\n\$file\\n\\n========\\n\\n"
                exit 1
            fi
        done
        cp "\$FILENAME"/* .
        rm -r \$FILENAME
    else
        if [[ ! "\${SINGLE_FILE_NAME}" =~ \\.(fasta|fa|fas|fna|txt|seq|faa)\$ ]] ; then
            echo -e "\\n\\n========\\n\\nERROR IN NEXTFLOW EXECUTION\\n\\nALL THE UNZIPPED FILE EXTENSIONS MUST BE fasta, fa, fas, fna, txt, seq OR faa.\\nHERE A FILE IS:\\n\${SINGLE_FILE_NAME}.\\n\\n========\\n\\n"
            exit 1
        fi
    fi
    """
}