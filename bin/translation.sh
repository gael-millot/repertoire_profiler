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



select_ch=${1}
aligned_seq=${2}
igblast_aa=${3}


if [[ "${igblast_aa}" == "false" ]] ; then
    # translation into aa (for a potential second round of analysis) since analysis is performed at the nuc level
    mkdir aa
    # translate fasta files
        seqkit translate -T 1 -f 1 --allow-unknown-codon ${aligned_seq} > aa/${aligned_seq}_tempo |& tee -a translation.log
            # no trim, no translate unknown code to 'X'
            # -T 1 : human genetic code
            # -f 1 : only the first frame is translated
            # --allow-unknown-codon : convert unknown codon (for instance ...) to X
        awk 'BEGIN{ORS=""}{if($0~/^>.*/){if(NR>1){print "\n"} ; print $0"\n"} else {print $0 ; next}}END{print "\n"}' aa/${aligned_seq}_tempo > aa/${i2} |& tee -a translation.log # remove \n
        rm aa/${aligned_seq}_tempo |& tee -a translation.log
    # end translate fasta files
    # assemble name and aa seq
    awk '{
        lineKind=(NR-1)%2 ; 
        if(lineKind==0){
            gsub("> *", "", $0)
            print $0 >> "name.txt"
        }
        if(lineKind==1){
            print $0 >> "seq.txt"
        }
    }' aa/${aligned_seq} |& tee -a translation.log
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



