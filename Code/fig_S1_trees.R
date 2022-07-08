---
  # Paper: Cinnocuum and its diversity
  # Author: DB
  # Figure S1: tree
  
---
  
# Figure S1: 16S and phylophlan tree
  
# Input files: RAxML_bipartitionsBranchLabels.16S_final, RAxML_bipartitionsBranchLabels.phylo
  
# Output files: plot figures
  
# Libraries
  
library(stringr)
library(stringi)
library(vegan)
library(readxl)
library(cowplot)
library(tidyverse)
library(ggplot2)
library(readr)
library(gplots)
library(ggtree)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(plyr)
library(treeio)
  
# Figure S1A
raxml_16S <- read.raxml("~/Box Sync/Disha/Projects/carb_utilization/cinnocuum/RAxML_bipartitionsBranchLabels.16S_final")

ggtree(raxml_16S, layout = "rectangular") +
  geom_tiplab(size = 1.5) + geom_text2(aes(label = bootstrap), size = 2, color = "red") +
  theme_tree()

# Figure S1B
raxml_cinn_cdiff_faa2 <- read.raxml("~/Box Sync/Disha/Projects/carb_utilization/cinnocuum/raxml_phylo/RAxML_bipartitionsBranchLabels.phylo")

ggtree(raxml_cinn_cdiff_faa2) +
  geom_tiplab(size = 1.5) +
  theme_tree() +
  geom_text2(aes(label = bootstrap), size = 2, color = "red") +
  theme(plot.margin = unit(c(14,8,14,8), "mm")) +
  xlim(-0.1, 38)


  
  
  
  

  
  