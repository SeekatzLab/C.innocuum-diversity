---

  # Paper: Cinnocuum and its diversity
  # Author: DB
  # Figure 4: Cazyme analysis
  # Distribution of cazymes per clade and heat map per strain showing the no. of genes present

---

# Figure 4A: Total number of cazymes and their distribution per clade

# Figue 4B: Heatmap with clustering trees on strains

# Input files:
#   - all_dbcan_cinn.txt
#   - cluster_groups2.txt

# Output files: figure plots saved

# Libraries used in the entire script

library(stringr)
library(stringi)
library(readxl)
library(tidyverse)
library(ggplot2)
library(readr)
library(gplots)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(plyr)
library(pheatmap)

# Figure 4A: cazyme distribution

## Making the full cazyme df for all species

### sample id file from output dbcan folder
#file_list <- read_tsv(file = "~/Box Sync/Disha/Projects/carb_utilization/cinnocuum/output_dbcan_latest/sample2")
### file list as a vector
#file_list <- as.vector(file_list$list)
### create an empty list which will be populated with all the overview.txt
#list.data<-list()
### file path to overview.txt from run_dbcan
#filepath = "~/Box Sync/Disha/Projects/carb_utilization/cinnocuum/output_dbcan_latest/dbcan_"

#for(i in (1:length(file_list))){
  #print(i)
#  list.data[[i]] <- read_tsv(paste0(filepath, file_list[i], '/overview.txt'))
#}

### combining all the overview.txt from a list to one large data frame and simultaneously adding their species name

# CM06_79_S180 <- list.data[[1]]
# CM06_79_S180$genome_name <- "CM06_79_S180"

### bind all the separate txt files from the list together and have the species name associated with it; 

# all_dbcan_cinn <- CM06_79_S180

# for (i in 2:113) {
#  all_dbcan_cinn <- rbind.fill(all_dbcan_cinn, list.data[[i]])
#  all_dbcan_cinn <- mutate(all_dbcan_cinn, genome_name = replace(genome_name, is.na(genome_name), file_list[[i]]))
# }

### writing the entire df as txt file, that will be used as input for the rest
# write_csv(all_dbcan_cinn, file = "~/Box Sync/Disha/Projects/carb_utilization/cinnocuum/output_dbcan_latest/all_dbcan_cinn.txt")

 all_dbcan_cinn <- read_csv(file = "~/Box Sync/Disha/Projects/carb_utilization/cinnocuum/output_dbcan_latest/all_dbcan_cinn.txt")

### modifying all_dbcan_cinn; shortening the genome name for easier use later on
all_dbcan_cinn$newgenome_name <- sub("^([^_]*_[^_]*)_.*$", "\\1", all_dbcan_cinn$genome_name)
all_dbcan_cinn$newgenome_name <- gsub("\\..*","",all_dbcan_cinn$newgenome_name)
all_dbcan_cinn$newgenome_name[all_dbcan_cinn$genome_name == "Ref_Cinnocuum_14501"] <- "Cinnocuum_14501"
all_dbcan_cinn$newgenome_name[all_dbcan_cinn$genome_name == "Ref_Cinnocuum_LCLUMC"] <- "Cinnocuum_LCLUMC"
all_dbcan_cinn$newgenome_name[all_dbcan_cinn$genome_name == "Ref_Cinnocuum_2959"] <- "Cinnocuum_2959"
all_dbcan_cinn$newgenome_name[all_dbcan_cinn$genome_name == "Ref_Cinnocuum_I46"] <- "Cinnocuum_I46"

### modifying the all_dbcan_cinn to separate out the rows that have 2 families in the same row separated by '+'
all_dbcan_cinn1 <- separate_rows(all_dbcan_cinn, HMMER, Hotpep, DIAMOND, sep = "\\+.*")
### Selecting only those that have more than one tool : various papers do this in their cazy db section
all_dbcan_cinn2 <- subset(all_dbcan_cinn1, `#ofTools` > 1) 

### removing the blank spaces and the '-' from the 3 tools hmmer, hotpep and diamond so it doesnot interfere in the coalesce stage
all_dbcan_cinn2$HMMER <- replace(all_dbcan_cinn2$HMMER, grepl("^\\s*$", all_dbcan_cinn2$HMMER) == TRUE, NA)
all_dbcan_cinn2$HMMER <- replace(all_dbcan_cinn2$Hotpep, grepl("^\\s*$", all_dbcan_cinn2$Hotpep) == TRUE, NA)
all_dbcan_cinn2$DIAMOND <- replace(all_dbcan_cinn2$DIAMOND, grepl("^\\s*$", all_dbcan_cinn2$DIAMOND) == TRUE, NA)
all_dbcan_cinn2$HMMER <- replace(all_dbcan_cinn2$HMMER, grepl("^\\-*$", all_dbcan_cinn2$HMMER) == TRUE, NA)
all_dbcan_cinn2$Hotpep <- replace(all_dbcan_cinn2$Hotpep, grepl("^\\-*$", all_dbcan_cinn2$Hotpep) == TRUE, NA)
all_dbcan_cinn2$DIAMOND <- replace(all_dbcan_cinn2$DIAMOND, grepl("^\\-*$", all_dbcan_cinn2$DIAMOND) == TRUE, NA)

### coalesce the three tools together to generate 1 consensus between three tools under cazyfam

super_dbcan <- all_dbcan_cinn2
super_dbcan <- mutate(super_dbcan, cazyfam= coalesce(HMMER, Hotpep, DIAMOND))

### extract the first 2 letters of the cazyfam to identify the type of cazy

super_dbcan$cazyfaminit <- substr(super_dbcan$cazyfam, 1, 2)

### assigning category to cazyfam

super_dbcan$cazyfamcategory[super_dbcan$cazyfaminit=="GT"] <- "GlycosylTransferases"
super_dbcan$cazyfamcategory[super_dbcan$cazyfaminit=="GH"] <- "Glycoside Hydrolases"
super_dbcan$cazyfamcategory[super_dbcan$cazyfaminit=="CB"] <- "Carbohydrate-Binding Modules"
super_dbcan$cazyfamcategory[super_dbcan$cazyfaminit=="CE"] <- "Carbohydrate Esterases"
super_dbcan$cazyfamcategory[super_dbcan$cazyfaminit=="PL"] <- "Polysaccharide Lyases"
super_dbcan$cazyfamcategory[super_dbcan$cazyfaminit=="AA"] <- "Auxiliary Activities"

### separating out all bracketed values associated with cazy fam and subfam so only fam and subfam can be counted at the end

super_dbcan2 <- super_dbcan
super_dbcan2$cazyfamid <- sub("\\(.*","", super_dbcan2$cazyfam)

### cluster grouping

cluster_group <- read_tsv(file = "~/Box Sync/Disha/Projects/carb_utilization/cinnocuum/cluster_groups2.txt")
cluster_group1 <- cluster_group %>%  arrange(clades)

cluster_group1$color[cluster_group1$clades == "clade_1"] <- "#E8A419"
cluster_group1$color[cluster_group1$clades == "clade_2"] <- "#9FC095"
cluster_group1$color[cluster_group1$clades == "clade_3"] <- "#3B99B1"
cluster_group1$color[cluster_group1$clades == "clade_4"] <- "#F5191C"

### selecting non NA rows from cazyfamid
super_dbcan3 <- subset(super_dbcan2, !is.na(super_dbcan2$cazyfamid))
mycols2 <- as.vector(cluster_group1$color)

### merge the all_meta and super_dbcan5
merge_meta <- inner_join(super_dbcan3, cluster_group1, by = c("newgenome_name" = "isolate"))

### now the entire file is ready to be counted and plotted; ordering all tables for correct coloring
summary_dbcan <- merge_meta %>% group_by(newgenome_name, cazyfamid, cazyfamcategory, color, clades) %>% tally() %>% arrange(clades)


## Figure 4A plot

ggplot(summary_dbcan, aes(y=cazyfamid, x=n, fill = clades)) + geom_col() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 3), axis.text.x = element_text(size = 3), panel.background = element_rect(colour = "black", size=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 3), legend.title = element_text(size = 3)) +
  scale_fill_manual(breaks = c("clade_1", "clade_2", "clade_3", "clade_4"), values = c("#E8A419", "#9FC095", "#3B99B1", "#F5191C"))


# Figure 4B: heatmap

### modifying the summary_dbcan into numric matrix
sum_mat <- summary_dbcan[,c(1,2,6)]

sum_mat1 <- sum_mat %>% pivot_wider(names_from = cazyfamid, values_from = n, values_fill = 0)
sum_mat2 <- sum_mat1[,-1]
rownames(sum_mat2) <- sum_mat1$newgenome_name

sum_mat3 <- t(sum_mat2)

### order the y axis for 4B in the sum_mat3 matrix
ordered.list <- as.character(sort(unique(summary_dbcan$cazyfamid)))
sum_mat3_ordered <- sum_mat3[ order(match(rownames(sum_mat3), rev(ordered.list))), ]

### for annotation_col
cluster_group2 <- cluster_group1[,c(1,3)] %>% arrange(clades, isolate)
cluster_group3 <- cluster_group2[,-1]
rownames(cluster_group3) <- colnames(sum_mat3)
cluster_group3 <- as.data.frame(cluster_group3)

### for annotation_colors
ann_colors <- list(clades = c(clade_1 = "#E8A419", clade_2 = "#9FC095", clade_3 = "#3B99B1", clade_4 = "#F5191C"))

## for heatmap colors/breaks:
cols <- brewer.pal(6, "Blues")
# to visualize colors:
# barplot(c(1:6), col=brewer.pal(6, "Blues"))
heat_colors <- c("#9E9E9E", cols[2:6])
heat_breaks <- c(0, 0.99, 1.99, 4.99, 9.99, 19.99, 30)


## Figure 4B plot
pheatmap(sum_mat3_ordered, cluster_rows = FALSE, annotation_col= cluster_group3, annotation_colors = ann_colors, 
	cutree_cols = 3, border_color = "black", breaks = heat_breaks, color = heat_colors, 
	legend_breaks = heat_breaks, clustering_method = "median", fontsize = 6)





