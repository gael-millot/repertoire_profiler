#!/usr/bin/env bash

#########################################################################
##                                                                     ##
##     TrimTranslate.sh                                                  ##
##                                                                     ##
##     Gael A. Millot                                                  ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################




select_ch=${1}


mkdir productive_nuc
mkdir productive_aa
if (( $(cat ${select_ch} | wc -l ) > 1 )) ; then
    SEQ="sequence" # name of the column containing the initial sequences
    SEQ_ALIGN="sequence_alignment" # name of the column containing the aligned sequence by igblast
    # make fasta files of the filtered sequences (only productive ones because this process is called after the productive filtering)
    awk -v var1=${SEQ} -v var2=${SEQ_ALIGN} 'BEGIN{FS="\t" ; ORS="\n" ; OFS="\t"}{
        if(NR==1){
            COL_NAME="FALSE"
            COL_SEQ="FALSE"
            COL_SEQ_ALIGN="FALSE"
            for(i4=1; i4<=NF; i4++){
                if($i4=="sequence_id"){COL_NAME=i4}
                if($i4==var1){COL_SEQ=i4}
                if($i4==var2){COL_SEQ_ALIGN=i4}
            }
            if(COL_NAME=="FALSE"){
                print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nNO sequence_id COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                exit 1
            }
            if(COL_SEQ=="FALSE"){
                print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nNO "var1" COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                exit 1
            }
            if(COL_SEQ_ALIGN=="FALSE"){
                print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nNO "var2" COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                exit 1
            }
        }else{
            if($COL_SEQ!~/^[-NATGC.]*$/){
                print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\n"var1" COLUMN NAME OF THE INPUT FILE MUST BE A NUCLEOTIDE SEQUENCE\nHERE IT MIGHT BE MADE OF AMINO ACIDS:\n"$COL_SEQ"\n\n========\n\n"
                exit 1
            }
            gsub(/\./, "", $COL_SEQ_ALIGN) # Remove dots from the sequence
            if($COL_SEQ_ALIGN!~/^[NATGC]*$/){
                print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\n"var2" COLUMN NAME OF THE INPUT FILE MUST BE A NUCLEOTIDE SEQUENCE\nHERE IT MIGHT BE MADE OF AMINO ACIDS:\n"$COL_SEQ_ALIGN"\n\n========\n\n"
                exit 1
            }
            print ">"$COL_NAME"\n"$COL_SEQ > $COL_NAME"_ini.fasta" # make fasta files of the nuc sequences
            print ">"$COL_NAME"\n"$COL_SEQ_ALIGN > $COL_NAME"_align.fasta" # make fasta files of the nuc sequences
            print $COL_NAME > "COL_NAME.txt" # name of the line in a file
        }
    }' ${select_ch} |& tee -a trimtranslate.log
    # end make fasta files of the productive sequences

    # triming the initial nucleotide sequence (corresponding to the leader peptide before the FWR1 region)
        # This assumes sequence B appears exactly once in sequence A. If there are multiple matches, it will use the last occurrence found.
        # Get sequence B as string
        COL_NAME=$(cat COL_NAME.txt) # recover the name of the line
        SEQ_B=$(seqkit seq -s "$COL_NAME"_align.fasta)
        # Find start position of B in A
        START_POS=$(seqkit locate -p "$SEQ_B" "$COL_NAME"_ini.fasta | tail -n 1 | cut -f5)
        # Trim A from the start position onwards
        if (( $START_POS > 1 )) ; then
            TRIM=TRUE
        else
            TRIM=FALSE # FALSE means no tream
        fi
        seqkit subseq -w 0 -r ${START_POS}:-1 "$COL_NAME"_ini.fasta > ./productive_nuc/"$COL_NAME"_trim.fasta
        TRIM_SEQ=$(sed -n '2p' ./productive_nuc/"$COL_NAME"_trim.fasta)
    # end triming the initial nucleotide sequence (corresponding to the leader peptide before the FWR1 region)
    # translation into aa (for a potential second round of analysis) since analysis is performed at the nuc level
    FILENAME=$(basename -- ./productive_nuc/*.*) # recover a file name without path. Here a single file
    # translate fasta files
    seqkit translate -T 1 -f 1 --allow-unknown-codon ./productive_nuc/${COL_NAME}_trim.fasta > ./productive_aa/${FILENAME}_tempo |& tee -a trimtranslate.log
        # no trim, no translate unknown code to 'X'
        # -T 1 : human genetic code
        # -f 1 : only the first frame is translated
        # --allow-unknown-codon : convert unknown codon (for instance ...) to X
    awk 'BEGIN{ORS=""}{if($0~/^>.*/){if(NR>1){print "\n"} ; print $0"\n"} else {print $0 ; next}}END{print "\n"}' productive_aa/${FILENAME}_tempo > productive_aa/${FILENAME} |& tee -a trimtranslate.log # remove \n
    rm productive_aa/${FILENAME}_tempo |& tee -a trimtranslate.log
    # end translate fasta files
    # assemble name and productive_aa seq
    awk '{
        lineKind=(NR-1)%2 ; 
        if(lineKind==0){
            gsub("> *", "", $0)
            print $0 >> "name.txt"
        }
        if(lineKind==1){
            print $0 >> "seq.txt"
        }
    }' productive_aa/${FILENAME} |& tee -a trimtranslate.log
    paste --delimiters='\t' name.txt seq.txt > aa.tsv |& tee -a trimtranslate.log
    # end assemble name and aa seq
    # add the aa seq into the trimtranslate.tsv
    awk -v var1=${SEQ} -v var2=${TRIM_SEQ} -v var3=${TRIM} 'BEGIN{FS="\t" ; ORS="" ; OFS=""}
        FNR==NR{ # means that work only on the first file
            a[$1] = $1 # a is name (sequence_id)
            b[$1] = $2 # b is sequence
            next
        }{ # mean that works only for the second file
            if(FNR==1){
                print $0"\tsequence_ini\tis_sequence_trimmed\tsequence_aa\n" > "trimtranslate.tsv" # header added to trimtranslate.tsv
                for(i4=1; i4<=NF; i4++){
                    if($i4=="sequence_id"){COL_NAME=i4}
                    if($i4==var1){COL_SEQ=i4}
                }
                # no need the recheck as above because already done above
            }else{
                if($COL_NAME in a){
                    for(i5=1; i5<=NF; i5++){ # instead of print $0 (to replace the initial sequence in the sequence column by the trimmed sequence
                        if($i5!=var1){
                            print $i5 > "trimtranslate.tsv"
                            if(i5!=FN){print "\t" > "trimtranslate.tsv"} # because of BEGIN
                        }else{
                            print var2 > "trimtranslate.tsv" # trimmed sequence replace the initial sequence in the SEQ column (i.e., "sequence" column)
                            if(i5!=FN){print "\t" > "trimtranslate.tsv"} # because of BEGIN
                        }
                    }
                    print $COL_SEQ"\t"var3"\t"b[$COL_NAME]"\n" > "trimtranslate.tsv"
                }
            }
        }
    ' aa.tsv ${select_ch} |& tee -a trimtranslate.log
    # end add the aa seq into the trimtranslate.tsv
    sed -i '1i sequence_id\tsequence_aa' aa.tsv |& tee -a trimtranslate.log # header added to aa.tsv
    # echo -e "sequence_id\tsequence_alignment\n" | cat aa.tsv > caca.tsv 
fi



