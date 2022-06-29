---

  # Paper: Cinnocuum and its diversity
  # Author: DB
  # Figure 3: core genome tree, PCA plot , cog % in pangenome, functional enrichment

---


# Fig 3A: core genome tree- from graphlan; Not mentioned in this script; see figure_3_metabolism

# Fig 3B: PCA plot no loadings

# Fig 3C: percentage of COG assignments in core, softshell, accessory-singletons

# Fig 3D: functional enrichment

# Input files:
  ## Fig 3B: full_cinnocuum.tsv; cluster_prokka.txt; cog-20.def.txt
  ## Fig 3C: NEWCINNWGS_gene_clusters_summary.txt
  ## Fig 3D: enriched.txt, functl_enr_fin.txt, list_enriched_mods.txt 


# Output files: plot files, list_enriched_mods.txt


# Libraries

library(stringr)
library(stringi)
library(readxl)
library(tidyverse)
library(ggplot2)
library(readr)
library(tidyr)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(philentropy)
library(vegan)
library(ade4)


# Fig 3B

### Creating a gene presence absence from prokka.tsv so we can have good gene name

prokka_full_tsv <- read_tsv(file = "~/Box Sync/Disha/Projects/carb_utilization/cinnocuum/full_cinnocuum.tsv")
prokka_full_tsv <- separate(prokka_full_tsv, col = locus_tag, into = c("prokka_id", "locus_no"), sep = "_")

prokka_meta <- read_tsv(file = "~/Box Sync/Disha/Projects/carb_utilization/cinnocuum/cluster_prokka.txt")

cog_def <- read_tsv(file = "~/Box Sync/Disha/Projects/carb_utilization/cog-20.def.txt")
cog_def <- cog_def[, 1:5]

cluster_group1 <- prokka_meta %>%  arrange(clades)

cluster_group1$color[cluster_group1$clades == "clade_1"] <- "#E8A419"
cluster_group1$color[cluster_group1$clades == "clade_2"] <- "#9FC095"
cluster_group1$color[cluster_group1$clades == "clade_3"] <- "#3B99B1"
cluster_group1$color[cluster_group1$clades == "clade_4"] <- "#F5191C"

merge_prokka <- inner_join(prokka_full_tsv, cluster_group1, by = 'prokka_id') %>% arrange(clades)

merge_prokka$prab <- 1

merge_prab <- select(merge_prokka, species, product, prab) %>% distinct()
merge_prab1 <- merge_prab %>% pivot_wider(names_from = species, values_from = prab, values_fill = 0)

merge_prab2 <- merge_prab1 
rownames(merge_prab2) <- merge_prab1$product
merge_prab3 <- merge_prab2[,-1]

rownames(merge_prab3) <- merge_prab2$product
merge_prab4 <- t(merge_prab3)

## create pca plot
set.seed(1234)


cog.pca <- dudi.pca(merge_prab4, scannf = FALSE, center = TRUE, scale = FALSE)


cog.eig <- cog.pca$eig
prop_var <- cog.eig/sum(cog.eig) * 100


pca.axes.cog <- cog.pca$li
pca.axes.cog$species <- rownames(pca.axes.cog)
pca.axes.cog1 <- inner_join(pca.axes.cog, cluster_group1, by = 'species') %>% arrange(clades)

pca.axes.cog1$isolate_id <- "0"
pca.axes.cog1$isolate_id[grep("^CM", pca.axes.cog1$species)] <- "this_study"
pca.axes.cog1$isolate_id[pca.axes.cog1$isolate_id == "0"] <- "other"

ggplot(pca.axes.cog1, aes(x= Axis1, y= Axis2, fill = clades, color = isolate_id)) +geom_point(shape =21) +
  #geom_segment(data = cog.loading.def, aes(x = 0, y = 0, xend = (CS1*20), yend = (CS2*20)), arrow = arrow(length = unit(1/2, "picas")), color = "blue") +
  #annotate("text", x = (cog.loading.def$CS1*20), y = (cog.loading.def$CS2*20),
  #         label = cog.loading.def$coggene, repel = TRUE) +
  theme_classic() +
  scale_fill_manual(breaks = pca.axes.cog1$clades, values = pca.axes.cog1$color) +
  scale_color_manual(breaks = pca.axes.cog1$isolate_id, values = c("black", rgb(0, 0, 0, alpha=0))) + 
  theme(text = element_text(size = 5), panel.background = element_rect(colour = "black", size=0.3), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "Axis1 (33.22%)", y = "Axis2 (15.49%)") +
  scale_x_continuous(limits = c(-5,15)) +
  scale_y_continuous(limits = c(-5,15)) +
  coord_fixed() +
  ggtitle("PCA of gene presence and absence COG")


# Fig 3C

## soft-core = soft-shell core

anvio_pan_summary <- read_tsv(file = "~/Box Sync/Disha/Projects/carb_utilization/cinnocuum/NEWCINNWGS_gene_clusters_summary.txt")

anvio_pan_summary1 <-  anvio_pan_summary %>% select(-c(aa_sequence))
anvio_pan_summary1$bin_name[is.na(anvio_pan_summary1$bin_name)] <- "core"

pan_sum <- anvio_pan_summary1 
pan_sum$COG20_CATEGORY_ACC <- gsub("\\|.*", "", pan_sum$COG20_CATEGORY_ACC)
pan_sum$COG20_CATEGORY <- gsub("\\|.*", "", pan_sum$COG20_CATEGORY)

## to analyze the main categories in core, soft shell and accessory
## energy_C <- pan_sum %>% select(bin_name, COG20_CATEGORY_ACC, COG20_FUNCTION, COG20_PATHWAY) %>% subset(COG20_CATEGORY_ACC == "C") %>% distinct()

## carb_G <- pan_sum %>% select(bin_name, COG20_CATEGORY_ACC, COG20_FUNCTION, COG20_PATHWAY, COG20_FUNCTION_ACC) %>% subset(COG20_CATEGORY_ACC == "G") %>% distinct()

## efhip <- pan_sum %>% select(bin_name, COG20_CATEGORY_ACC, COG20_FUNCTION, COG20_PATHWAY) %>% subset(COG20_CATEGORY_ACC == "E" | COG20_CATEGORY_ACC == "F" | COG20_CATEGORY_ACC == "H" | COG20_CATEGORY_ACC == "I" | COG20_CATEGORY_ACC == "P") %>% distinct()

## mob_x <- pan_sum %>% select(bin_name, COG20_CATEGORY_ACC, COG20_FUNCTION, COG20_PATHWAY) %>% subset(COG20_CATEGORY_ACC == "X") %>% distinct()

summ_anvio <- pan_sum %>% group_by(bin_name, COG20_CATEGORY_ACC, COG20_CATEGORY) %>% tally()

summ_anvio1 <- summ_anvio
summ_anvio1$COG20_CATEGORY_ACC[is.na(summ_anvio1$COG20_CATEGORY_ACC)] <- "Not_Assigned"
summ_anvio1$COG20_CATEGORY[is.na(summ_anvio1$COG20_CATEGORY)] <- "Not_Assigned"

summ_anvio1$newcog <- summ_anvio1$COG20_CATEGORY_ACC


#write_tsv(summ_anvio, "~/Box Sync/Disha/Projects/carb_utilization/cinnocuum/summ_anvio.txt")

## calculating percentages from summ_anvio2

core_summary <- subset(summ_anvio1, summ_anvio1$bin_name == "core")
#core_summary <- core_summary %>% group_by(bin_name, COG20_CATEGORY, newcog) %>% tally()
core_summary$perc <- round((core_summary$n/sum(core_summary$n)) *100, digits = 6)

accsin_summary <- subset(summ_anvio1, summ_anvio1$bin_name == "accessory_singletons")
#accsin_summary <- accsin_summary %>% group_by(bin_name, newcog) %>% tally()
accsin_summary$perc <- round((accsin_summary$n/sum(accsin_summary$n)) *100, digits = 6)

soft_summary <- subset(summ_anvio1, summ_anvio1$bin_name == "soft-core")
#soft_summary <- soft_summary %>% group_by(bin_name, newcog) %>% tally()
soft_summary$perc <- round((soft_summary$n/sum(soft_summary$n)) *100, digits = 6)

final_summary <- rbind(core_summary, accsin_summary, soft_summary)

#write_tsv(final_summary, "~/Box Sync/SeekatzLab/Grants/Paper1/FIGURES_2/full_summary.txt")

## making the color palette and sorting the data

color_pal <- as.data.frame(unique(final_summary$newcog))
colnames(color_pal) <- "newcog"
color_pal$color[color_pal$newcog == "R" | color_pal$newcog == "S" | color_pal$newcog == "W" | color_pal$newcog == "Not_Assigned"] <- "grey67"
color_pal$color[color_pal$newcog == "J" | color_pal$newcog == "K" | color_pal$newcog == "L"] <- "#A71B4B"
color_pal$color[color_pal$newcog == "D" | color_pal$newcog == "V" | color_pal$newcog == "T" | color_pal$newcog == "M" | color_pal$newcog == "N" | color_pal$newcog == "O" | color_pal$newcog == "U"] <- "#D04939"
color_pal$color[color_pal$newcog == "C"] <- "#EB7803"
color_pal$color[color_pal$newcog == "G"] <- "#F9BC53"
color_pal$color[color_pal$newcog == "E"] <- "#FEF1A6"
color_pal$color[color_pal$newcog == "F"] <- "#E2F8B5"
color_pal$color[color_pal$newcog == "H"] <- "#9CE5AD"
color_pal$color[color_pal$newcog == "I"] <- "#43CBB1"
color_pal$color[color_pal$newcog == "P"] <- "#00AAB6"
color_pal$color[color_pal$newcog == "Q"] <- "#0080B2"
color_pal$color[color_pal$newcog == "X"] <- "#584B9F"

fin_sum <- inner_join(final_summary, color_pal, by = "newcog")
ord <- c("Not_Assigned", "R", "S", "W", "J", "K", "L", "D", "V", "T", "M", "N", "O", "U", "C", "G", "E", "F", "H", "I", "P", "Q", "X")
fin_sum1 <- fin_sum %>% arrange(factor(newcog, levels = ord))

ggplot(fin_sum1, aes(x=bin_name, y=perc, fill =factor(newcog, levels = ord))) + geom_col(width = 0.4) +
  theme_classic() +
  theme(panel.background = element_rect(colour = "black", size=1), text = element_text(size = 5), legend.key.size = unit(5, "mm")) +
  scale_fill_manual(breaks = fin_sum1$newcog, values = fin_sum1$color)


# Fig 3D

## read in data files for functional enrichment

enriched_KEGG_Class <- read_tsv(file = "~/Box Sync/Disha/Projects/carb_utilization/cinnocuum/enriched.txt")
enriched_KEGG_Modules <- read_tsv(file = "~/Box Sync/Disha/Projects/carb_utilization/cinnocuum/functl_enr_fin.txt")

## data modification of kegg_modules

enr_mods <- enriched_KEGG_Modules %>% subset(!(is.na(associated_groups)))
enr_mods1 <- enr_mods[ ,1:11]
enr_mods2 <- pivot_longer(enr_mods1, cols = !c("KEGG_Module", "enrichment_score", "unadjusted_p_value", "adjusted_q_value", "associated_groups", "accession", "gene_clusters_ids"), names_to = "clade", values_to = "fraction")
enr_mods3 <- separate_rows(enr_mods2, KEGG_Module, sep = "!!!")

## ggplot

ggplot(enr_mods3, aes(x=clade, y=KEGG_Module, fill=fraction)) + geom_raster() +
  scale_fill_gradient2(high = "navyblue", low = "white")

## data modification of kegg_class

enr_class <- enriched_KEGG_Class %>% subset(!(is.na(associated_groups)))
enr_class1 <- enr_class[ ,1:11]
enr_class2 <- pivot_longer(enr_class1, cols = !c("KEGG_Class", "enrichment_score", "unadjusted_p_value", "adjusted_q_value", "associated_groups", "accession", "gene_clusters_ids"), names_to = "clade", values_to = "fraction")
enr_class3 <- separate_rows(enr_class2, KEGG_Class, sep = "!!!")

## ggplot

ggplot(enr_class3, aes(x=clade, y=KEGG_Class, fill=fraction)) + geom_raster() +
  scale_fill_gradient2(high = "navyblue", low = "white")

## coloring according to KEGG_Class and COG or kofam if possible

enriched_modules <- as.data.frame(unique(enr_mods3$KEGG_Module))

write_tsv(enriched_modules, "~/Box Sync/Disha/Projects/carb_utilization/cinnocuum/list_enriched_mods.txt")

## I added Kegg_class_cat (where the class category of the KEGG module is listed) and a color for the class_cat in MSExcel; it is easier that way.
list_enr_mod <- read_tsv(file = "~/Box Sync/Disha/Projects/carb_utilization/cinnocuum/list_enriched_mods.txt")


merge2 <- inner_join(enr_mods3, list_enr_mod, by = "KEGG_Module") %>%  arrange(KEGG_Class_Cat)
list_enr_mod1 <- arrange(list_enr_mod, KEGG_Class_Cat)
mycols <- as.vector(list_enr_mod1$color)

## You can add the color according to the class_cat but I removed it so the figure can be fixed in Illustrator
ggplot(merge2, aes(x=clade, y=KEGG_Module, fill=fraction)) + geom_tile(color = "black") + theme_bw() +
  scale_fill_gradient2(high = "navyblue", low = "white") + coord_fixed(ratio = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 2), legend.key.size = unit(2, 'mm')) +
  scale_y_discrete(limits = list_enr_mod1$KEGG_Module)


