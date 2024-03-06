#!/usr/bin/env bash

#########################################################################
##                                                                     ##
##     translation.sh                                                  ##
##                                                                     ##
##     Gael A. Millot                                                  ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################




aligned_seq=${1}
igblast_aa=${2}


if [[ "${igblast_aa}" == "false" ]] ; then
    # translation into aa (for a potential second round of analysis) since analysis is performed at the nuc level
    mkdir aa
    # translate fasta files
    for i2 in $(ls aligned_seq) ; do
        seqkit translate -T 1 -f 1 --allow-unknown-codon aligned_seq/$i2 > aa/${i2}_tempo |& tee -a translation.log
            # no trim, no translate unknown code to 'X'
            # -T 1 : human genetic code
            # -f 1 : only the first frame is translated
            # --allow-unknown-codon : convert unknown codon (for instance ...) to X
        awk 'BEGIN{ORS=""}{if($0~/^>.*/){if(NR>1){print "\n"} ; print $0"\n"} else {print $0 ; next}}END{print "\n"}' aa/${i2}_tempo > aa/${i2} |& tee -a translation.log # remove \n
        rm aa/${i2}_tempo |& tee -a translation.log
    done
    # end translate fasta files
    # assemble name and aa seq
    for i2 in $(ls aa) ; do
        awk '{
            lineKind=(NR-1)%2 ; 
            if(lineKind==0){
                gsub("> *", "", $0)
                print $0 >> "name.txt"
            }
            if(lineKind==1){
                print $0 >> "seq.txt"
            }
        }' aa/${i2} |& tee -a translation.log
    done
    paste --delimiters='\t' name.txt seq.txt > aa.tsv |& tee -a translation.log
    # end assemble name and aa seq
    # add the aa seq into the translation.tsv
    awk '
        FNR==NR{ # means that work only on the first file
            var1=$1
            a[var1] = var1 # a is name
            b[var1] = $2 # b is sequence
            next
        }{ # mean that works only for the second file
            if(FNR==1){
                print $0"\tsequence_alignment_aa" > "translation.tsv"
                for(i4=1; i4<=NF; i4++){
                    if($i4=="sequence_id"){COL_NAME=i4}
                }
                # no need the recheck as above because already done above
            }else{
                if($COL_NAME in a){print $0"\t"b[$COL_NAME] > "translation.tsv"}
            }
        }
    ' aa.tsv ${select_ch} |& tee -a translation.log
    # end add the aa seq into the translation.tsv
    sed -i '1i sequence_id\tsequence_alignment' aa.tsv |& tee -a translation.log # header added
    # echo -e "sequence_id\tsequence_alignment\n" | cat aa.tsv > caca.tsv 
else
    cat ${select_ch} > translation.tsv |& tee -a translation.log
fi



