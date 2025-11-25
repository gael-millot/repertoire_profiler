#!/usr/bin/env bash

#########################################################################
##                                                                     ##
##     TrimTranslate.sh                                                ##
##                                                                     ##
##     Gael A. Millot                                                  ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################




select_ch=${1}

set -o pipefail

mkdir -p productive_nuc/trimmed
mkdir -p productive_nuc/removed
mkdir -p productive_nuc/query
mkdir -p productive_nuc/align
mkdir -p productive_nuc/align_with_gaps

mkdir -p productive_aa/trimmed
mkdir -p productive_aa/igblast
mkdir -p productive_aa/query
mkdir -p productive_aa/align

if (( $(cat ${select_ch} | wc -l ) > 1 )) ; then
    SEQ="sequence" # name of the column containing the initial sequences
    SEQ_ALIGN="sequence_alignment" # name of the column containing the aligned sequence by igblast
    SEQ_ALIGN_GAP="sequence_alignment_with_gaps" # name of the column containing the aligned sequence by igblast
    GERM_ALIGN_GAP="germline_alignment_with_gaps" # name of the column containing the aligned sequence by igblast
    SEQ_AA="sequence_aa"
    SEQ_ALIGN_AA="sequence_alignment_aa"
    # make fasta files of the filtered sequences (only productive ones because this process is called after the productive filtering)
    awk -v var1=${SEQ} -v var2=${SEQ_ALIGN} -v var3=${SEQ_ALIGN_GAP} -v var4=${GERM_ALIGN_GAP} -v var5=${SEQ_AA} -v var6=${SEQ_ALIGN_AA} 'BEGIN{FS="\t" ; ORS="\n" ; OFS="\t"}{
        if(NR==1){
            NAME_COL_NB="FALSE"
            SEQ_COL_NB="FALSE"
            SEQ_ALIGN_COL_NB="FALSE"
            SEQ_ALIGN_GAP_COL_NB="FALSE"
            GERM_ALIGN_GAP_COL_NB="FALSE"
            SEQ_AA_COL_NB="FALSE"
            SEQ_ALIGN_AA_COL_NB="FALSE"
            for(i4=1; i4<=NF; i4++){
                if($i4=="sequence_id"){NAME_COL_NB=i4}
                if($i4==var1){SEQ_COL_NB=i4}
                if($i4==var2){SEQ_ALIGN_COL_NB=i4}
                if($i4==var3){SEQ_ALIGN_GAP_COL_NB=i4}
                if($i4==var4){GERM_ALIGN_GAP_COL_NB=i4}
                if($i4==var5){SEQ_AA_COL_NB=i4}
                if($i4==var6){SEQ_ALIGN_AA_COL_NB=i4}
            }
            if(NAME_COL_NB=="FALSE"){
                print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nFOR "$NAME_COL_NB", NO sequence_id COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                exit 1
            }
            if(SEQ_COL_NB=="FALSE"){
                print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nFOR "$NAME_COL_NB", NO "var1" COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                exit 1
            }
            if(SEQ_ALIGN_COL_NB=="FALSE"){
                print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nFOR "$NAME_COL_NB", NO "var2" COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                exit 1
            }
            if(SEQ_ALIGN_GAP_COL_NB=="FALSE"){
                print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nFOR "$NAME_COL_NB", NO "var3" COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                exit 1
            }
            if(GERM_ALIGN_GAP_COL_NB=="FALSE"){
                print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nFOR "$NAME_COL_NB", NO "var4" COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                exit 1
            }
            if(SEQ_AA_COL_NB=="FALSE"){
                print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nFOR "$NAME_COL_NB", NO "var5" COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                exit 1
            }
            if(SEQ_ALIGN_AA_COL_NB=="FALSE"){
                print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nFOR "$NAME_COL_NB", NO "var6" COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                exit 1
            }
        }else{
            if($SEQ_COL_NB!~/^[NATGC]*$/){
                print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nFOR "$NAME_COL_NB", "var1" COLUMN NAME OF THE INPUT FILE MUST BE A NUCLEOTIDE SEQUENCE WITHOUT . OR - \nHERE IT MAY BE AN ALIGNED SEQUENCE OR A AMINO ACIDS SEQUENCE:\n"$SEQ_COL_NB"\n\n========\n\n"
                exit 1
            }
            if($SEQ_ALIGN_COL_NB!~/^[-NATGC.]*$/){
                print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nFOR "$NAME_COL_NB", "var2" COLUMN NAME OF THE INPUT FILE MUST BE A NUCLEOTIDE SEQUENCE\nHERE IT MIGHT BE MADE OF AMINO ACIDS:\n"$SEQ_COL_NB"\n\n========\n\n"
                exit 1
            }
            gsub(/\./, "", $SEQ_ALIGN_COL_NB) # Remove dots from the sequence
            if($SEQ_ALIGN_COL_NB!~/^[-NATGC]*$/){
                print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nFOR "$NAME_COL_NB", "var2" COLUMN NAME OF THE INPUT FILE MUST BE A NUCLEOTIDE SEQUENCE\nHERE IT MIGHT BE MADE OF AMINO ACIDS:\n"$SEQ_ALIGN_COL_NB"\n\n========\n\n"
                exit 1
            }
            if($SEQ_ALIGN_GAP_COL_NB!~/^[-NATGC.]*$/){
                print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nFOR "$NAME_COL_NB", "var3" COLUMN NAME OF THE INPUT FILE MUST BE A NUCLEOTIDE SEQUENCE WITH POTENTIAL IMGT GAPS (DOTS)\nHERE IT MIGHT BE MADE OF AMINO ACIDS:\n"$SEQ_COL_NB"\n\n========\n\n"
                exit 1
            }
            TEMPO_SEQ_ALIGN_GAP=$SEQ_ALIGN_GAP_COL_NB
            gsub(/\./, "", TEMPO_SEQ_ALIGN_GAP)
            if(TEMPO_SEQ_ALIGN_GAP!~/^[-NATGC]*$/){
                print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nFOR "$NAME_COL_NB", "var3" COLUMN NAME OF THE INPUT FILE MUST BE A NUCLEOTIDE SEQUENCE WITH POTENTIAL IMGT GAPS (DOTS)\nHERE IT MIGHT BE MADE OF AMINO ACIDS:\n"TEMPO_SEQ_ALIGN_GAP"\n\n========\n\n"
                exit 1
            }
            if($GERM_ALIGN_GAP_COL_NB!~/^[-NATGC.]*$/){
                print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nFOR "$NAME_COL_NB", "var4" COLUMN NAME OF THE INPUT FILE MUST BE A NUCLEOTIDE SEQUENCE WITH POTENTIAL IMGT GAPS (DOTS)\nHERE IT MIGHT BE MADE OF AMINO ACIDS:\n"$SEQ_COL_NB"\n\n========\n\n"
                exit 1
            }
            print ">"$NAME_COL_NB"\n"$SEQ_COL_NB > "./productive_nuc/query/"$NAME_COL_NB"_query.fasta" # make fasta files of the query nuc sequences
            print ">"$NAME_COL_NB"\n"$SEQ_ALIGN_COL_NB > "./productive_nuc/align/"$NAME_COL_NB"_align.fasta" # make fasta files of the aligned nuc sequences without dots. This sequence will help to know if trimming must be done. Because igblast only returns aligned VDJ seq (without leader seq).
            print ">"$NAME_COL_NB"\n"$SEQ_ALIGN_GAP_COL_NB > "./productive_nuc/align_with_gaps/"$NAME_COL_NB"_seq_align_with_gaps.fasta"
            # print ">"$NAME_COL_NB"\n"$GERM_ALIGN_GAP_COL_NB > "./productive_nuc/align_with_gaps/"$NAME_COL_NB"_germ_align_with_gaps.fasta"
            print ">"$NAME_COL_NB"\n"$SEQ_AA_COL_NB > "./productive_aa/igblast/aa_"$NAME_COL_NB".fasta" # make fasta files of the sequence_aa column
            print ">"$NAME_COL_NB"\n"$SEQ_ALIGN_AA_COL_NB > "./productive_aa/align/aa_"$NAME_COL_NB"_align.fasta"
            print $NAME_COL_NB > "NAME.txt" # name of the line in a file
        }
    }' ${select_ch} |& tee -a trimtranslate.log
    # end make fasta files of the productive sequences

    # triming the initial nucleotide sequence (corresponding to the leader peptide before the FWR1 region)
        # This assumes sequence B appears exactly once in sequence A. If there are multiple matches, it will use the last occurrence found.
        # Get sequence B as string
        NAME=$(cat NAME.txt) # recover the name of the line
        SEQ_B=$(seqkit seq -s ./productive_nuc/align/"${NAME}"_align.fasta | sed 's/-//g') #call seq and remove ---
        # Find start position of B in A
        START_POS=$(seqkit locate -p "$SEQ_B" ./productive_nuc/query/"${NAME}"_query.fasta | tail -n 1 | cut -f5) # determine if trimming has been done by igblast when returning aligned VDJ seq
        if ! [[ "$START_POS" =~ ^-?[0-9]+$ ]]; then
            print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nSTART_POS MUST BE AN INTEGER OF POSITION.\nHERE IT IS:\n"$START_POS"\n\n========\n\n"
            exit 1
        fi
        # Trim A from the start position onwards
        if (( $START_POS > 1 )) ; then
            TRIM_LOG=TRUE
        else
            TRIM_LOG=FALSE # FALSE means no tream
        fi
        seqkit subseq -w 0 -r ${START_POS}:-1 ./productive_nuc/query/"${NAME}"_query.fasta > ./productive_nuc/trimmed/"${NAME}"_trim.fasta
        # -w 0 → disables line wrapping so the sequence is on one line.
        # -r ${START_POS}:-1 → specifies the range: from position ${START_POS} up to the last base (-1 means the end).
        TRIM_SEQ=$(sed -n '2p' ./productive_nuc/trimmed/"${NAME}"_trim.fasta) # prints only the second line of the file.
        seqkit subseq -w 0 -r 1:$(( ${START_POS} - 1 )) ./productive_nuc/query/"${NAME}"_query.fasta > ./productive_nuc/removed/"${NAME}"_removed.fasta
        # -w 0 → disables line wrapping so the sequence is on one line.
        # -r ${START_POS}:-1 → specifies the range: from position ${START_POS} up to the last base (-1 means the end).
        REMOVED_SEQ=$(sed -n '2p' ./productive_nuc/removed/"${NAME}"_removed.fasta) # prints only the second line of the file.
    # end triming the initial nucleotide sequence (corresponding to the leader peptide before the FWR1 region)
    # translation into aa (for a potential second round of analysis) since analysis is performed at the nuc level
    # translate fasta files
    seqkit translate -T 1 -f 1 -w 0 --allow-unknown-codon ./productive_nuc/trimmed/"${NAME}"_trim.fasta > ./productive_aa/trimmed/aa_"${NAME}"_trim.fasta |& tee -a trimtranslate.log
    seqkit translate -T 1 -f 1 -w 0 --allow-unknown-codon ./productive_nuc/query/"${NAME}"_query.fasta > ./productive_aa/query/aa_"${NAME}"_query.fasta |& tee -a trimtranslate.log
        # no trim, no translate unknown code to 'X'
        # -T 1 : human genetic code
        # -f 1 : only the first frame is translated
        # --allow-unknown-codon : convert unknown codon (for instance ...) to X
    # end translate fasta files
    INI_AA_SEQ=$(sed -n '2p' ./productive_aa/query/aa_"${NAME}"_query.fasta)
    TRIM_AA_SEQ=$(sed -n '2p' ./productive_aa/trimmed/aa_"${NAME}"_trim.fasta)
    # assemble name and seq into aa.tsv 
    # awk '{
    #     lineKind=(NR-1)%2 ; 
    #     if(lineKind==0){
    #         gsub("> *", "", $0)
    #         print $0 >> "name.txt"
    #     }
    #     if(lineKind==1){
    #         print $0 >> "seq.txt"
    #     }
    # }' productive_aa/aa_"${NAME}"_trim.fasta |& tee -a trimtranslate.log
    # paste --delimiters='\t' name.txt seq.txt > aa.tsv |& tee -a trimtranslate.log
    # sed -i '1i sequence_id\tsequence_aa' aa.tsv |& tee -a trimtranslate.log # header added to aa.tsv
    # assemble name and seq into aa.tsv 
    # add the aa seq into the trimtranslate.tsv
    awk -v var1=${SEQ} -v var2=${TRIM_SEQ} -v var3=${TRIM_LOG} -v var4=${REMOVED_SEQ} -v var5=${TRIM_AA_SEQ} -v var6=${INI_AA_SEQ} -v var7=${SEQ_AA} 'BEGIN{FS="\t" ; ORS="" ; OFS=""}
        {
            if(FNR==1){ # first line
                gsub(/\r/, "") # remove CR
                print $0"\ttrimmed_sequence\tis_sequence_trimmed\tremoved_sequence\tremoved_sequence_length\ttrimmed_sequence_aa\tquery_sequence_aa\talign_seq_identical\taa_identical\tsequence_aa_stop\tsequence_alignment_aa_stop\tgermline_alignment_aa_stop\ttrimmed_sequence_aa_stop\tquery_sequence_aa_stop\tfwr1_alignment_start\tfwr1_alignment_end\tcdr1_alignment_start\tcdr1_alignment_end\tfwr2_alignment_start\tfwr2_alignment_end\tcdr2_alignment_start\tcdr2_alignment_end\tfwr3_alignment_start\tfwr3_alignment_end\tcdr3_alignment_start\tcdr3_alignment_end\tfwr4_alignment_start\tfwr4_alignment_end\n" > "trimtranslate.tsv"
                # header added to trimtranslate.tsv
                NAME_COL_NB=i4="FALSE"
                SEQ_COL_NB="FALSE"
                SEQ_AA_COL_NB="FALSE"
                fwr1_start_NB="FALSE"
                fwr1_end_NB="FALSE"
                cdr1_start_NB="FALSE"
                cdr1_end_NB="FALSE"
                fwr2_start_NB="FALSE"
                fwr2_end_NB="FALSE"
                cdr2_start_NB="FALSE"
                cdr2_end_NB="FALSE"
                fwr3_start_NB="FALSE"
                fwr3_end_NB="FALSE"
                cdr3_start_NB="FALSE"
                cdr3_end_NB="FALSE"
                fwr4_start_NB="FALSE"
                fwr4_end_NB="FALSE"
                v_alignment_start_NB="FALSE"
                SEQ_ALIGN_AA_COL_NB="FALSE"
                GERM_ALIGN_AA_COL_NB="FALSE"
                SEQ_ALIGN_COL_NB="FALSE"
                SEQ_ALIGN_GAP_COL_NB="FALSE"
                for(i4=1; i4<=NF; i4++){
                    if($i4=="sequence_id"){NAME_COL_NB=i4}
                    if($i4==var1){SEQ_COL_NB=i4}
                    if($i4==var7){SEQ_AA_COL_NB=i4}
                    if($i4=="fwr1_start"){fwr1_start_NB=i4}
                    if($i4=="fwr1_end"){fwr1_end_NB=i4}
                    if($i4=="cdr1_start"){cdr1_start_NB=i4}
                    if($i4=="cdr1_end"){cdr1_end_NB=i4}
                    if($i4=="fwr2_start"){fwr2_start_NB=i4}
                    if($i4=="fwr2_end"){fwr2_end_NB=i4}
                    if($i4=="cdr2_start"){cdr2_start_NB=i4}
                    if($i4=="cdr2_end"){cdr2_end_NB=i4}
                    if($i4=="fwr3_start"){fwr3_start_NB=i4}
                    if($i4=="fwr3_end"){fwr3_end_NB=i4}
                    if($i4=="cdr3_start"){cdr3_start_NB=i4}
                    if($i4=="cdr3_end"){cdr3_end_NB=i4}
                    if($i4=="fwr4_start"){fwr4_start_NB=i4}
                    if($i4=="fwr4_end"){fwr4_end_NB=i4}
                    if($i4=="v_alignment_start"){v_alignment_start_NB=i4}
                    if($i4=="sequence_alignment_aa"){SEQ_ALIGN_AA_COL_NB=i4}
                    if($i4=="germline_alignment_aa"){GERM_ALIGN_AA_COL_NB=i4}
                    if($i4=="sequence_alignment"){SEQ_ALIGN_COL_NB=i4}
                    if($i4=="sequence_alignment_with_gaps"){SEQ_ALIGN_GAP_COL_NB=i4}
                }
                if(SEQ_COL_NB=="FALSE"){
                    print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nNO "var1" COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                    exit 1
                }
                if(SEQ_AA_COL_NB=="FALSE"){
                    print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nNO "var7" COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                    exit 1
                }
                if(fwr1_start_NB=="FALSE"){
                    print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nNO fwr1_start COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                    exit 1
                }
                if(fwr1_end_NB=="FALSE"){
                    print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nNO fwr1_end COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                    exit 1
                }
                if(cdr1_start_NB=="FALSE"){
                    print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nNO cdr1_start COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                    exit 1
                }
                if(cdr1_end_NB=="FALSE"){
                    print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nNO cdr1_end COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                    exit 1
                }
                if(fwr2_start_NB=="FALSE"){
                    print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nNO fwr2_start COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                    exit 1
                }
                if(fwr2_end_NB=="FALSE"){
                    print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nNO fwr2_end COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                    exit 1
                }
                if(cdr2_start_NB=="FALSE"){
                    print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nNO cdr2_start COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                    exit 1
                }
                if(cdr2_end_NB=="FALSE"){
                    print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nNO cdr2_end COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                    exit 1
                }
                if(fwr3_start_NB=="FALSE"){
                    print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nNO fwr3_start COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                    exit 1
                }
                if(fwr3_end_NB=="FALSE"){
                    print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nNO fwr3_end COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                    exit 1
                }
                if(cdr3_start_NB=="FALSE"){
                    print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nNO cdr3_start COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                    exit 1
                }
                if(cdr3_end_NB=="FALSE"){
                    print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nNO cdr3_end COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                    exit 1
                }
                if(fwr4_start_NB=="FALSE"){
                    print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nNO fwr4_start COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                    exit 1
                }
                if(fwr4_end_NB=="FALSE"){
                    print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nNO fwr4_end COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                    exit 1
                }
                if(v_alignment_start_NB=="FALSE"){
                    print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nNO v_alignment_start COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                    exit 1
                }
                if(SEQ_ALIGN_AA_COL_NB=="FALSE"){
                    print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nNO sequence_alignment_aa COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                    exit 1
                }
                if(GERM_ALIGN_AA_COL_NB=="FALSE"){
                    print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nNO germline_alignment_aa COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                    exit 1
                }
                if(SEQ_ALIGN_COL_NB=="FALSE"){
                    print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nFOR "$NAME_COL_NB", NO "var2" COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                    exit 1
                }
                if(SEQ_ALIGN_GAP_COL_NB=="FALSE"){
                    print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nFOR "$NAME_COL_NB", NO "var3" COLUMN NAME FOUND IN THE INPUT FILE\n\n========\n\n"
                    exit 1
                }
                # no need the recheck as above because already done above
            }else{
                gsub(/\r/, "") # remove CR
                # inactivated because now "sequence" column is the query one and "trimmed_sequence" is the trimmed one
                # for(i5=1; i5<=NF; i5++){ # instead of print $0 (to replace the initial sequence in the sequence column by the trimmed sequence
                #     if(i5!=SEQ_COL_NB){
                #         print $i5 > "trimtranslate.tsv"
                #         if(i5!=FN){print "\t" > "trimtranslate.tsv"} # because of BEGIN
                #     }else{
                #         OLD_SEQ=$i5
                #         print var2 > "trimtranslate.tsv" # trimmed sequence replace the initial sequence in the SEQ column (i.e., "sequence" column)
                #         if(i5!=FN){print "\t" > "trimtranslate.tsv"} # because of BEGIN
                #     }
                # }
                REMOVED_SEQ_LENGTH=length(var4)
                if($v_alignment_start_NB!="NA"){
                    if($v_alignment_start_NB!=1){
                        print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nv_alignment_start SHOULD BE 1 OR NA AND NOT "$v_alignment_start_NB" FOR "$NAME_COL_NB"\n\n========\n\n"
                        exit 1
                    }
                }
                if($fwr1_start_NB!="NA"){
                    DIFF=REMOVED_SEQ_LENGTH+1-$fwr1_start_NB
                    if(DIFF!=0){
                        print "\n\n========\n\nERROR IN NEXTFLOW EXECUTION OF THE TrimTranslate PROCESS\n\nremoved_sequence_length + 1 - fwr1_start SHOULD BE EQUAL TO ZERO AND NOT "DIFF" FOR "$NAME_COL_NB"\n\n========\n\n"
                        exit 1
                    }
                }
                if($SEQ_AA_COL_NB==var5){
                    IDENTICAL="TRUE"
                }else{
                    IDENTICAL="FALSE"
                }
                if($SEQ_AA_COL_NB ~ /\*/){
                    SEQ_AA_STOP="TRUE"
                }else{
                    SEQ_AA_STOP="FALSE"
                }
                if($SEQ_ALIGN_AA_COL_NB ~ /\*/){
                    SEQ_ALIGN_AA_STOP="TRUE"
                }else{
                    SEQ_ALIGN_AA_STOP="FALSE"
                }
                if($GERM_ALIGN_AA_COL_NB ~ /\*/){
                    GERM_ALIGN_AA_STOP="TRUE"
                }else{
                    GERM_ALIGN_AA_STOP="FALSE"
                }
                if(var5 ~ /\*/){
                    TRIM_AA_SEQ_STOP="TRUE"
                }else{
                    TRIM_AA_SEQ_STOP="FALSE"
                }
                if(var6 ~ /\*/){
                    INI_AA_SEQ_STOP="TRUE"
                }else{
                    INI_AA_SEQ_STOP="FALSE"
                }
                TEMPO_SEQ_ALIGN_GAP=$SEQ_ALIGN_GAP_COL_NB
                gsub(/\./, "", TEMPO_SEQ_ALIGN_GAP)
                if(TEMPO_SEQ_ALIGN_GAP==$SEQ_ALIGN_COL_NB){
                    IDENT_SEQ_ALIGN="TRUE"
                }else{
                    IDENT_SEQ_ALIGN="FALSE"
                }
                if($fwr1_start_NB!="NA"){fwr1_start2=$fwr1_start_NB-REMOVED_SEQ_LENGTH}else{fwr1_start2="NA"}
                if($fwr1_end_NB!="NA"){fwr1_end2=$fwr1_end_NB-REMOVED_SEQ_LENGTH}else{fwr1_end2="NA"}
                if($cdr1_start_NB!="NA"){cdr1_start2=$cdr1_start_NB-REMOVED_SEQ_LENGTH}else{cdr1_start2="NA"}
                if($cdr1_end_NB!="NA"){cdr1_end2=$cdr1_end_NB-REMOVED_SEQ_LENGTH}else{cdr1_end2="NA"}
                if($fwr2_start_NB!="NA"){fwr2_start2=$fwr2_start_NB-REMOVED_SEQ_LENGTH}else{fwr2_start2="NA"}
                if($fwr2_end_NB!="NA"){fwr2_end2=$fwr2_end_NB-REMOVED_SEQ_LENGTH}else{fwr2_end2="NA"}
                if($cdr2_start_NB!="NA"){cdr2_start2=$cdr2_start_NB-REMOVED_SEQ_LENGTH}else{cdr2_start2="NA"}
                if($cdr2_end_NB!="NA"){cdr2_end2=$cdr2_end_NB-REMOVED_SEQ_LENGTH}else{cdr2_end2="NA"}
                if($fwr3_start_NB!="NA"){fwr3_start2=$fwr3_start_NB-REMOVED_SEQ_LENGTH}else{fwr3_start2="NA"}
                if($fwr3_end_NB!="NA"){fwr3_end2=$fwr3_end_NB-REMOVED_SEQ_LENGTH}else{fwr3_end2="NA"}
                if($cdr3_start_NB!="NA"){cdr3_start2=$cdr3_start_NB-REMOVED_SEQ_LENGTH}else{cdr3_start2="NA"}
                if($cdr3_end_NB!="NA"){cdr3_end2=$cdr3_end_NB-REMOVED_SEQ_LENGTH}else{cdr3_end2="NA"}
                if($fwr4_start_NB!="NA"){fwr4_start2=$fwr4_start_NB-REMOVED_SEQ_LENGTH}else{fwr4_start2="NA"}
                if($fwr4_end_NB!="NA"){fwr4_end2=$fwr4_end_NB-REMOVED_SEQ_LENGTH}else{fwr4_end2="NA"}
                print $0"\t"var2"\t"var3"\t"var4"\t"REMOVED_SEQ_LENGTH"\t"var5"\t"var6"\t"IDENT_SEQ_ALIGN"\t"IDENTICAL"\t"SEQ_AA_STOP"\t"SEQ_ALIGN_AA_STOP"\t"GERM_ALIGN_AA_STOP"\t"TRIM_AA_SEQ_STOP"\t"INI_AA_SEQ_STOP"\t"fwr1_start2"\t"fwr1_end2"\t"cdr1_start2"\t"cdr1_end2"\t"fwr2_start2"\t"fwr2_end2"\t"cdr2_start2"\t"cdr2_end2"\t"fwr3_start2"\t"fwr3_end2"\t"cdr3_start2"\t"cdr3_end2"\t"fwr4_start2"\t"fwr4_end2"\n" > "trimtranslate.tsv"
            }
        }
    ' ${select_ch} 
    # end add the aa seq into the trimtranslate.tsv
fi



