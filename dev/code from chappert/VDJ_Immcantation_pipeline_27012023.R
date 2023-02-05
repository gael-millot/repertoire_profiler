#'VDJ analysis of anti-RBD Memory B cells in naive vaccinated patients (+/- Omicron):
#'Dedicated to culture sequences only. 

#######################################################
## Part1: Generation of aligned VDJ repertoire database
#######################################################
#'from sequences to post-ChangeO db 

library(Biostrings)
library(data.table)
library(plyr)
library(alakazam)
library(shazam)
library(openxlsx)

setwd("~/YourFolder")
folder <- "Analysis_AllDonors_XXXX"
if(isFALSE(dir.exists(folder))){
  dir.create(folder)
}

## Create FASTA file with all sequences:
########################################
memB_seq <- read.xlsx("All_seq_Naives_Donors_XXXX.xlsx", 1)
key_attributes <- c("cell_id", "donor_id", "donor_id_anno", "time_point", "cell_type", "sort_specificity", "assay")
table(memB_seq$donor_id, memB_seq$time_point)
memB_seq$cell_id[duplicated(memB_seq$cell_id)] #check for duplicated cell_id (issue later at the MakeDb step)

seq <- memB_seq$sequence
names(seq) <- memB_seq$cell_id
dna = DNAStringSet(seq)
fasta_file <- paste0(folder,"/All_donors_sequences.fasta")
writeXStringSet(dna, fasta_file)

rm(dna)

## Generate AIRR Rearrangement data (ChangeO/IgBlast):
######################################################
## !!! to be run in Terminal:
cd ~/YourFolder/Analysis_AllDonors_XXXX

#' Step1: IgBlast for VDJ aligment (AssignGenes/MakeDb):
AssignGenes.py igblast -s All_donors_sequences.fasta -b ~/share/igblast --organism human --loci ig --format blast
MakeDb.py igblast -i All_donors_sequences_igblast.fmt7 -s All_donors_sequences.fasta -r ~/share/germlines/imgt/human/vdj/imgt_human_*.fasta --extended

#' Step2: Blast for Ig_C gene alignment (blastn):
blastn -query All_donors_sequences.fasta -task blastn -out All_donors_Ig_C_genes_results.csv -outfmt "10 qseqid sseqid score ppos qlen qstart sstart length" -max_target_seqs 1 -db ~/YourFolder/Ig_C_seq_Hg38.fasta 
#' !! end Terminal

## Add c_call and additional attributes to culture_airr.tsv table:
######################################################
culture_airr <- read.table(paste0(folder, "/All_donors_sequences_igblast_db-pass.tsv"), sep = '\t', header = TRUE)
culture_airr$cell_id <- culture_airr$sequence_id

#failed sequences:
memB_seq$cell_id[!memB_seq$cell_id %in% culture_airr$cell_id]
write.csv(memB_seq[!memB_seq$cell_id %in% culture_airr$cell_id,], paste0(folder, "/Failed_IgBlast-MakeDb_sequences.csv"))

#add isotype information:
c_results <- read.csv(paste0(folder, "/All_donors_Ig_C_genes_results.csv"), 
                      col.names = c("sequence_id", "c_call_top_match", "c_call_alignment_score", "p_pos", "sequence_length", "sequence_align_start", "IGH_align_start", "c_call_alignment_length"),
                      header = FALSE)
c_results <- c_results[c_results$IGH_align_start == 1,] #remove match within the Ig C gene sequence (should start at first base)
table(duplicated(c_results$sequence_id)) #check if there still is any duplicated match (should not occur as Blast only return the top match)

c_infos <- read.csv("~/YourFolder/Ig_C_seq_Hg38_infos.csv")
c_results$c_call <- c_infos$Isotype[match(c_results$c_call_top_match, c_infos$ENST_id)]
c_results <- c_results[c_results$c_call_alignment_length >= 35,] #keep only long enough match for c_call (<35bp not enough to distinguish IGHG1 and IGHG2)

culture_airr <- merge(culture_airr, 
                      c_results[, names(c_results) %in% c("sequence_id", "c_call", "c_call_alignment_score", "c_call_alignment_length")], 
                      by = "sequence_id", all.x = TRUE)

# add additional attributes:
culture_airr <- merge(culture_airr, memB_seq[, names(memB_seq) %in% key_attributes], 
                      by = "cell_id", all.x = TRUE)

write.table(culture_airr, file = paste0(folder, "/All_donors_sequences_igblast_db-pass_extended.tsv"), row.names = FALSE, sep = "\t")

## calculating distance to nearest:
###################################

# Calculates nearest neighbor distances:
culture_airr <- distToNearest(culture_airr, cellIdColumn="cell_id", sequenceColumn="junction", locusColumn="locus", model="ham", normalize="len", nproc=1)

#Plot results:
pdf(paste0(folder, "/All_donors_sequences_distance_to_nearest.pdf"))
ggplot(subset(culture_airr, !is.na(dist_nearest)), aes(x=dist_nearest)) + theme_bw() + 
  xlab("Hamming distance") + ylab("Count") +
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=0.15, color="firebrick", linetype=2) 
dev.off()

culture_airr$dist_change0 <- 0.15

write.table(culture_airr, file = paste0(folder, "/All_donors_sequences_full.tsv"), row.names = FALSE, sep = "\t")

rm(c_results, c_infos, memB_seq, culture_airr)

## define clones (ChangeO):
###########################
#' !!! to be run in Terminal:
#' !!! DefineClones will change collumns names and remove all Capitals letters...
#' Step 1: Define clones based on heavy chain only:
#' Step 2: Add masked germline-dmask sequence (useful for SHM calculations)
DefineClones.py -d All_donors_sequences_full.tsv --act set --model ham --norm len --dist 0.15
CreateGermlines.py -d All_donors_sequences_full_clone-pass.tsv -g dmask --cloned -r ~/share/germlines/imgt/human/vdj/imgt_human_*.fasta
#' !! end Terminal

## calculate Mutational load:
###############################
#' Clonal assignment and germline sequences reconstruction should have been performed 
#' using the DefineClone.py and CreateGerlines.py in ChangeO
#' A "germline_alignment_d_mask" collumn should be present. 
#' If germline sequences reconstruction has been performed after clonal assignment,
#' a single germline_alignment_d_mask" consensus sequence should be present for each clone.

VDJ_db <- read.table(paste0(folder, "/All_donors_sequences_full_clone-pass_germ-pass.tsv"), sep = '\t', header = TRUE)

# Calculate R and S mutation counts
VDJ_db <- observedMutations(VDJ_db, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=IMGT_V,
                            frequency=FALSE, 
                            nproc=1)
# Calculate combined R and S mutation counts
VDJ_db <- observedMutations(VDJ_db, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=IMGT_V,
                            frequency=FALSE, 
                            combine=TRUE,
                            nproc=1)

# Calculate R and S mutation frequencies
VDJ_db <- observedMutations(VDJ_db, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=IMGT_V,
                            frequency=TRUE, 
                            nproc=1)

# Calculate combined R and S mutation frequencies
VDJ_db <- observedMutations(VDJ_db, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=IMGT_V,
                            frequency=TRUE, 
                            combine=TRUE,
                            nproc=1)

VDJ_db <- arrange(VDJ_db, clone_id)

write.table(VDJ_db, file = paste0(folder, "/All_donors_sequences_full_clone-pass_germ-pass_shm-pass.tsv"), row.names = FALSE, sep = "\t")
write.csv(VDJ_db, paste0(folder, "/All_donors_sequences_full_clone-pass_germ-pass_shm-pass.csv"))

interesting_collumns <- c("clone_id",	"sequence_id", "sequence",	
                          "donor_id", "donor_id_anno", "time_point", "cell_type", "sort_specificity",	"assay",
                          "v_call",	"d_call",	"j_call",	"cdr3", "junction", "junction_aa", "junction_length",	"c_call",	
                          "dist_nearest",	"dist_change0",	
                          "mu_count", "mu_count_cdr_r",	"mu_count_cdr_s",	"mu_count_fwr_r",	"mu_count_fwr_s",
                          "c_call_alignment_length")
setdiff(interesting_collumns, colnames(VDJ_db))

##recap table by clones:
library(openxlsx)
OUT <- createWorkbook()
addWorksheet(OUT, "All_sequences_non_productive")
writeData(OUT, sheet = "All_sequences_non_productive", x = VDJ_db[VDJ_db$productive == FALSE, interesting_collumns], colNames = TRUE, rowNames = FALSE)
addWorksheet(OUT, "All_sequences_productive")
writeData(OUT, sheet = "All_sequences_productive", x = VDJ_db[VDJ_db$productive == TRUE, interesting_collumns], colNames = TRUE, rowNames = FALSE)
donors <- levels(as.factor(VDJ_db$donor_id))
for(i in (1:length(donors))){
  donor <- donors[[i]]
  db_prod <- VDJ_db[VDJ_db$productive == TRUE & VDJ_db$donor_id == donor,]
  clones <- db_prod[duplicated(db_prod$clone_id),]$clone_id
  clones <- unique(clones)
  addWorksheet(OUT, paste0("Clones_", donor))
  writeData(OUT, sheet = paste0("Clones_", donor), x = db_prod[db_prod$clone_id %in% clones, interesting_collumns], colNames = TRUE, rowNames = FALSE)
  addWorksheet(OUT, paste0("Uniques_", donor))
  writeData(OUT, sheet = paste0("Uniques_", donor), x = db_prod[!(db_prod$clone_id %in% clones), interesting_collumns], colNames = TRUE, rowNames = FALSE)
  db_prod <- arrange(db_prod, time_point, cell_type)
  addWorksheet(OUT, paste0("All_seq_", donor))
  writeData(OUT, sheet = paste0("All_seq_", donor), x = db_prod[, interesting_collumns], colNames = TRUE, rowNames = FALSE)
}
saveWorkbook(OUT, file = paste0(folder, "/All_sequences_0.15_recap.xlsx"), overwrite = TRUE)

rm(db_prod)


#######
#Trees:
#######
#'Run in Docker image or install igphyml via direct compiling (faster than Docker image)
library(dowser)
library(ggtree)
setwd(paste0(analysis_folder, "/", integration_folder))
VDJ_db <- read.table("All_donors_sequences_full_clone-pass_germ-pass_shm-pass.tsv", sep = '\t', header = TRUE)

if(isFALSE(dir.exists("EGC&GC/Trees"))){
  dir.create("EGC&GC/Trees")
}

##IgPhyML trees
db <- VDJ_db[VDJ_db$productive == TRUE & !is.na(VDJ_db$productive),]

donors <- levels(as.factor(db$donor_id))
for (i in seq_along(donors)){
  AIRR <- db[db$donor_id == donors[i],]
  clones <- formatClones(AIRR, traits=c("time_point"), text_fields=c("cell_id"), mask_char = ".", max_mask = 18, minseq=5)
  if(length(clones$clone_id) != 0){
    trees <- getTrees(clones, build="igphyml", exec="~/igphyml/src/igphyml", nproc=4)
    save(trees, file=paste0("EGC&GC/Trees/", donors[i], "_igphyml_trees.Rdata"))
    plots <- plotTrees(trees, tips="numeric_time_point", tipsize="collapse_count") 
    treesToPDF(plots, file=paste0("EGC&GC/Trees/", donors[i], "_trees.pdf"), nrow=2, ncol=2)
  } else {warning(paste0(donors[i], " doesn't have clones of size >= 5 observed during that time frame"))}
}

##AA CDR3 logos:
#' run seperatly from Docker image
library(ggseqlogo)
plots_folder <- "CDR3_logo_post"
if(isFALSE(dir.exists(plots_folder))){
  dir.create(plots_folder)
}

clones <- formatClones(db, traits=c("numeric_time_point"), text_fields=c("cell_id"), mask_char = ".", max_mask = 18, minseq=5)
clone_ids <- clones$clone_id
for(j in seq_along(clone_ids)){
  junction_aa = db[db$donor_clone_id==clone_ids[j],]$junction_aa
  l = db[db$donor_clone_id==clone_ids[j],]$junction_length
  g = ggseqlogo(junction_aa, seq_type = "aa", method = "prob")
  ggsave(g, filename = paste0(plots_folder, "/", clone_ids[j],"CDR3_logo.pdf"), width=(2/8+l[1]/8), height=2)
}





