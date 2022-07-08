---

  # Paper : C. innocum and its diversity
  # Figure 6: toxins, virulence factors + blast output
  # Author: DB

---


## Input files

### Figure 6A: Folder: figure_6
#    - pathofact_pred.tsv
#    - tox_lib.tsv
#    - cluster_groups2.txt

### Figure 6B: Folder: figure_6
#    - vir_facs_orf_id_14501.tsv
#    - tox_orf_id_14501.tsv
#    - blast_Ref_Cinnocuum_14501.txt; blast_amr_mge_Ref_Cinnocuum_14501.txt; 
#    - Toxin_gene_library_Ref_Cinnocuum_14501_dvf_report.tsv
#    - Ref_Cinnocuum_14501.gtf

## Output files

### Figure 6B - vir_facs_orf_id_14501.tsv; tox_orf_id_14501.tsv; amr_pathofact_orf_id.tsv; path_vir2.txt; path_amr.txt; path_tox.txt; gi.csv


rm(list = ls())
dev.off()

# Libraries

library(stringr)
library(stringi)
library(readxl)
library(readr)
library(tidyverse)
library(ggplot2)
library(readr)
library(tidyr)
library(plyr)
library(RColorBrewer)
library(dplyr)

# Toxins: use toxin data from Pathofact
# data input, using for loop to catch all the files 
# with all predictions first

#file_list_tox <- read_tsv(file = "~/Box Sync/Disha/Projects/carb_utilization/cinnocuum/Pathofact_reports_all/sample2")

#file_list_tox <- as.vector(file_list_tox$list)

#list.data.tox<-list()

#filepath_tox = "~/Box Sync/Disha/Projects/carb_utilization/cinnocuum/Pathofact_reports_all/PathoFact_"

#for(i in (1:length(file_list_tox))){
#  #print(i)
#  list.data.tox[[i]] <- read_tsv(paste0(filepath_tox, file_list_tox[i], '_predictions.tsv'))
#}

# combining all the predictions.tsv from a list to one large data frame and simultaneously adding their species name

#CM01_52_S208 <- list.data.tox[[1]]
#CM01_52_S208$genome_name <- "CM1_52_S208"


#pathofact_pred <- CM01_52_S208

#for (i in 2:113) {
#  pathofact_pred <- rbind.fill(pathofact_pred, list.data.tox[[i]])
#  pathofact_pred <- mutate(pathofact_pred, genome_name = replace(genome_name, is.na(genome_name), file_list_tox[[i]]))
#}

#write_tsv(pathofact_pred, file = "~/Box Sync/SeekatzLab/Grants/Paper1/github/Data/figure_6/pathofact_pred.tsv")

# all toxin gene libraries are next

#list.data.tox.lib<-list()

#filepath.tox.lib = "~/Box Sync/Disha/Projects/carb_utilization/cinnocuum/Pathofact_reports_all/Toxin_gene_library_"

#for(i in (1:length(file_list_tox))){
#  #print(i)
#  list.data.tox.lib[[i]] <- read_tsv(paste0(filepath.tox.lib, file_list_tox[i], '_report.tsv'))
#}

# combining all the overview.txt from a list to one large data frame and simultaneously adding their species name

#CM1_52_S208.tox.lib <- list.data.tox.lib[[1]]
#CM1_52_S208.tox.lib$genome_name <- "CM1_52_S208"

## This peice of code I tried out binds all the cmit together and I have the species name associated with it; this will save me time in the future

#tox.lib <- CM1_52_S208.tox.lib

#for (i in 2:113) {
#  tox.lib <- rbind.fill(tox.lib, list.data.tox.lib[[i]])
#  tox.lib <- mutate(tox.lib, genome_name = replace(genome_name, is.na(genome_name), file_list_tox[[i]]))
#}

#write_tsv(tox.lib, file = "~/Box Sync/SeekatzLab/Grants/Paper1/github/Data/figure_6/tox_lib.tsv")

## Read in the following files

pathofact_pred <- read_tsv(file = "~/Box Sync/SeekatzLab/Grants/Paper1/github/Data/figure_6/tox_lib.tsv")
tox.lib <- read_tsv(file = "~/Box Sync/SeekatzLab/Grants/Paper1/github/Data/figure_6/tox_lib.tsv")

pathofact_pred$newgenome_name <- sub("^([^_]*_[^_]*)_.*$", "\\1", pathofact_pred$genome_name)
pathofact_pred$newgenome_name[pathofact_pred$genome_name == "Ref_Cinnocuum_14501_dvf"] <- "Cinnocuum_14501"
pathofact_pred$newgenome_name[pathofact_pred$genome_name == "Ref_Cinnocuum_LCLUMC_dvf"] <- "Cinnocuum_LCLUMC"
pathofact_pred$newgenome_name[pathofact_pred$genome_name == "Ref_Cinnocuum_2959_dvf"] <- "Cinnocuum_2959"
pathofact_pred$newgenome_name[pathofact_pred$genome_name == "Ref_Cinnocuum_I46_dvf"] <- "Cinnocuum_I46"

# merging the 2 db without removing anything; secreted toxin and non secreted toxin

merge_path <- inner_join(pathofact_pred, tox.lib[,-2], by = c("genome_name" = "genome_name", "ORF_ID" = "ORF_ID"))

# lost the other "-" because the toxin gene library removed it
# separating out secreted and non secreted toxin

secreted_toxins <- subset(merge_path, merge_path$Toxin_confidence_level == "1: Secreted Toxin")
non_secreted_toxins <- subset(merge_path, merge_path$Toxin_confidence_level == "2: Non-secreted Toxin")

# clade coloring

cluster_group <- read_tsv(file = "~/Box Sync/Disha/Projects/carb_utilization/cinnocuum/cluster_groups2.txt")
cluster_group1 <- cluster_group %>%  arrange(clades)

cluster_group1$color[cluster_group1$clades == "clade_1"] <- "#E8A419"
cluster_group1$color[cluster_group1$clades == "clade_2"] <- "#9FC095"
cluster_group1$color[cluster_group1$clades == "clade_3"] <- "#3B99B1"
cluster_group1$color[cluster_group1$clades == "clade_4"] <- "#F5191C"

# plotting

clust_sec_tox <- inner_join(secreted_toxins, cluster_group1, by = c("newgenome_name" = "isolate")) %>% arrange(clades)

clust_non_sec_tox <- inner_join(non_secreted_toxins, cluster_group1, by = c("newgenome_name" = "isolate")) %>% arrange(clades) %>% distinct()

cutoff_clust_sec_tox <- subset(clust_sec_tox, clust_sec_tox$Score >=50)
cutoff_clust_non_sec_tox <- subset(clust_non_sec_tox, clust_non_sec_tox$Score >=50)

ccnst_sum <- aggregate(cutoff_clust_non_sec_tox$Score, list(cutoff_clust_non_sec_tox$clades, cutoff_clust_non_sec_tox$NAME), FUN=mean)

ggplot(cutoff_clust_sec_tox, aes(x=newgenome_name, y=NAME, fill=Score)) + geom_tile() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(limits = unique(cutoff_clust_sec_tox$newgenome_name)) +
  scale_fill_gradient2(low = "yellow", high= "navy") +
  #coord_fixed() + 
  ggtitle("Secreted toxins")

ggplot(cutoff_clust_non_sec_tox, aes(x=newgenome_name, y=NAME, fill=Score)) + geom_tile() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(limits = unique(cutoff_clust_non_sec_tox$newgenome_name)) +
  scale_fill_gradient2(low = "yellow", high= "navy") +
  coord_fixed() + 
  ggtitle("Non-secreted toxins")

ggplot(ccnst_sum, aes(x=Group.1, y=Group.2, fill=x)) + geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  #scale_x_discrete(limits = unique(cutoff_clust_non_sec_tox$newgenome_name)) +
  scale_fill_gradient2(low = "white", high= "navy") +
  coord_fixed() + 
  ggtitle("Non-secreted toxins")


# new data set for secreted toxins with secreted virulence factors

secreted_vir_tox <- subset(merge_path, merge_path$Toxin_confidence_level == "1: Secreted Toxin" & merge_path$Virulence_confidence_level == "1: Secreted Virulence factor")
m.secreted_vir_tox <- inner_join(secreted_vir_tox, cluster_group1, by = c("newgenome_name" = "isolate")) %>% arrange(clades)

ggplot(m.secreted_vir_tox, aes(x=clades, y=NAME, fill=Score)) + geom_tile()

m.number_sec_vt <- group_by(clust_non_sec_tox, NAME, clades, newgenome_name) %>% tally()

# final figure for toxins with highest category, mid category and non secreted toxins of selected toxins

tox_list <- c('tcdAB', 'entD', 'entB', 'tlyC', 'Phage_holin_4', 'holin_tox_secr', 'YcfA', 'ToxN_toxin', 'RelB', 'RelE_StbE', 'hlyIII', 'hcnC', 'hcnB', 'GGGtGRT')

finals_tox <- merge_path[merge_path$NAME %in% tox_list,] %>% distinct()
m.finals_tox <- inner_join(finals_tox, cluster_group1, by = c("newgenome_name" = "isolate")) %>% arrange(factor(NAME, levels = tox_list))

# plotting Figure 6A:
ggplot(m.finals_tox, aes(x=clades, y=NAME, fill=Score)) +geom_tile() +
  theme_bw() +
  scale_y_discrete(limits = unique(m.finals_tox$NAME)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient(high = "#000080", low = "#E78100") +
  coord_fixed()


# Figure 6B: Circos plot for C. innocuum 14501


## subset pathofact_pred to Ref Cinnocuum 14501 to get just the 14501
## make orf_id list to blast against the genome to get coordinates on the genome

### virulence factors
cinn14501_vir <- subset(pathofact_pred, pathofact_pred$newgenome_name == "Cinnocuum_14501")
cinn14501_vir2 <- subset(cinn14501_vir, cinn14501_vir$Virulence_confidence_level == "1: Secreted Virulence factor" | cinn14501_vir$Virulence_confidence_level == "2: Non-secreted Virulence factor")

vir_df <- as.data.frame(cinn14501_vir2$ORF_ID) %>% distinct()
colnames(vir_df) <- "orf_id"

write_tsv(vir_df, file = "~/Box Sync/SeekatzLab/Grants/Paper1/github/Data/figure_6/vir_facs_orf_id_14501.tsv")

### toxin

cinn14501_tox2 <- subset(cinn14501_vir, cinn14501_vir$Toxin_confidence_level == "1: Secreted Toxin" | cinn14501_vir$Toxin_confidence_level == "2: Non-secreted Toxin")

tox_df <- as.data.frame(cinn14501_tox2$ORF_ID) %>% distinct()
colnames(tox_df) <- "orf_id"

write_tsv(tox_df, file = "~/Box Sync/SeekatzLab/Grants/Paper1/github/Data/figure_6/tox_orf_id_14501.tsv") 

### antimicrobial resistance genes

## creating one bing df for all amr_mge_predictions.tsv
#filepath_amr = "~/Box Sync/Disha/Projects/carb_utilization/cinnocuum/Pathofact_reports_all/AMR_MGE_prediction_"

#list.data.amr <- list()
#for(i in (1:length(file_list_tox))){
#  #print(i)
#  list.data.amr[[i]] <- read_tsv(paste0(filepath_amr, file_list_tox[i], '_report.tsv'))
#}

# combining all the predictions.tsv from a list to one large data frame and simultaneously adding their species name

#all.amr.mge <- list.data.amr[[1]]
#all.amr.mge$genome_name <- "CM1_52_S208"


#for (i in 2:113) {
#  all.amr.mge <- rbind.fill(all.amr.mge, list.data.amr[[i]])
#  all.amr.mge <- mutate(all.amr.mge, genome_name = replace(genome_name, is.na(genome_name), file_list_tox[[i]]))
#}

#amr_mge2 <- subset(all.amr.mge, !(all.amr.mge$ARG == "-"))

#amr_mge2$newgenome_name <- sub("^([^_]*_[^_]*)_.*$", "\\1", amr_mge2$genome_name)
#amr_mge2$newgenome_name[amr_mge2$genome_name == "Ref_Cinnocuum_14501_dvf"] <- "Cinnocuum_14501"
#amr_mge2$newgenome_name[amr_mge2$genome_name == "Ref_Cinnocuum_LCLUMC_dvf"] <- "Cinnocuum_LCLUMC"
#amr_mge2$newgenome_name[amr_mge2$genome_name == "Ref_Cinnocuum_2959_dvf"] <- "Cinnocuum_2959"
#amr_mge2$newgenome_name[amr_mge2$genome_name == "Ref_Cinnocuum_I46_dvf"] <- "Cinnocuum_I46"

#write_tsv(amr_mge2, file = "~/Box Sync/SeekatzLab/Grants/Paper1/github/Data/figure_6/Antibiotic_resistance.tsv")

## generating the blast file from Ref_cinnocuum_14501_dvf_sl.faa to find the genes pathofact in the .faa is talking about; done on Palmetto
amr_mge2 <- read_tsv(file = "/Users/disha/Desktop/gbk_cinnocuum/coordinates/Antibiotic_resistance.tsv")
amr_mge2 <- subset(amr_mge2, amr_mge2$newgenome_name == "Cinnocuum_14501")
amr_orf_id <- amr_mge2 %>% select(ORF_ID) %>% distinct()

write_tsv(amr_orf_id, file = "/Users/disha/Desktop/gbk_cinnocuum/amr_pathofact_orf_id.tsv")
#### I need this file for blast against the Ref_Cinnocuum_14501.faa


## compare it to coordinates and produce file for circos to use

## for virulence factors after blast
vir_14501 <- read_tsv(file = "/Users/disha/Desktop/gbk_cinnocuum/blast_Ref_Cinnocuum_14501.txt", col_names = c("Query_ID", "gtf_id", "Percentage_of_identical_matches", "Alignment_length", "Number_of_mismatches", "Numberofgapopenings", "Startofalignmentinquery", "Endofalignmentinquery", "Startofalignmentinsubject", "Endofalignmentinsubject", "Expected_value", "Bit_score"))
final_vir <- subset(vir_14501, vir_14501$Percentage_of_identical_matches >= 99.00)
final_vir2 <- final_vir %>% select(Query_ID, gtf_id)

### for coordinates of virulence factors from pathofact
gtf_14501 <- read_tsv(file = "/Users/disha/Desktop/gbk_cinnocuum/coordinates/gtf_files/Ref_Cinnocuum_14501.gtf", col_names = c("name", "software", "gene_type", "start", "end", "blah1", "orientation", "blah2", "gene_id"))
gtf_14501<- separate(gtf_14501, col = gene_id, into = c("blah3", "gtf_id"), sep = " ")

final_vir3 <- inner_join(final_vir2, gtf_14501, by = "gtf_id")

write_tsv(final_vir3, file = "/Users/disha/Desktop/gbk_cinnocuum/data/path_vir2.txt") # <------ Goes into circos

## after blast antibiotic resistance

amr_14501 <- read_tsv(file = "/Users/disha/Desktop/gbk_cinnocuum/blast_amr_mge_Ref_Cinnocuum_14501.txt", col_names = c("Query_ID", "gtf_id", "Percentage_of_identical_matches", "Alignment_length", "Number_of_mismatches", "Numberofgapopenings", "Startofalignmentinquery", "Endofalignmentinquery", "Startofalignmentinsubject", "Endofalignmentinsubject", "Expected_value", "Bit_score"))
final_amr <- subset(amr_14501, amr_14501$Percentage_of_identical_matches >= 99.00)
final_amr2 <- final_amr %>% select(Query_ID, gtf_id)


final_amr3 <- inner_join(final_amr2, gtf_14501, by = "gtf_id")
# this goes into circos
write_tsv(final_amr3, file = "/Users/disha/Desktop/gbk_cinnocuum/data/path_amr.txt")

## after blast toxins

tox_14501 <- read_tsv(file = "/Users/disha/Desktop/gbk_cinnocuum/blast_tox_Ref_Cinnocuum_14501.txt", col_names = c("Query_ID", "gtf_id", "Percentage_of_identical_matches", "Alignment_length", "Number_of_mismatches", "Numberofgapopenings", "Startofalignmentinquery", "Endofalignmentinquery", "Startofalignmentinsubject", "Endofalignmentinsubject", "Expected_value", "Bit_score"))
final_tox <- subset(tox_14501, tox_14501$Percentage_of_identical_matches >= 99.00)
final_tox2 <- final_tox %>% select(Query_ID, gtf_id)
final_tox3 <- inner_join(final_tox2, gtf_14501, by = "gtf_id")

# this goes into circos
write_tsv(final_tox3, file = "/Users/disha/Desktop/gbk_cinnocuum/data/path_tox.txt") 


## for genomic islands
cinn14501_gi <- read_excel(path = "/Users/disha/Desktop/gbk_cinnocuum/genomic_islands/Cinn14501.xls")

cinn14501_gi2 <- cinn14501_gi %>% select(`Island start`, `Island end`) %>% distinct()

#this goes into circos after modifications

write_csv(cinn14501_gi2, file = "/Users/disha/Desktop/gbk_cinnocuum/genomic_islands/gi.csv")

## mods to output files before circos
#### remove the extra columns in Excel, put in chr1 on the first column and fill_color for the last column forr all output files before running through circos
### with given output files + circos config set to user preference visualization


# Figure S5: Antimicrobial resistance across clades
# antibiotic/antimicrobial peptide resistance

m.amr_mge2 <- inner_join(amr_mge2, cluster_group1, by = c('newgenome_name' = 'isolate'))
#plot
ggplot(m.amr_mge2, aes(x= clades, y = AMR_category)) + geom_tile() + coord_fixed() + theme_bw() 



