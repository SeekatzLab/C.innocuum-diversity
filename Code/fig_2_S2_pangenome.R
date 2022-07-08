---
  
  # Paper: Cinnocuum and its diversity
  # Author: DB
  # Figure 2, Figure S2: pangenome analysis

---
  
# Figure 2A : Number of genes vs total no of genomes for all genomes

# Figure 2B : No of unique genes vs total no of genomes for all genomes

# Figure S2A : No of genes vs total no of genomes for clade 1 & 2

# Figure S2B : No of genes vs total no of genomes for clade 1 & 2

# Input files :  number_of_genes_in_pan_genome.Rtab, gene_presence_absence.Rtab, gene_presence_absence.csv, clade_12_number_of_genes_in_pan_genome.Rtab, clade_12_gene_presence_absence.Rtab, clade_12_gene_presence_absence.csv
  
# Output files: plot figures

# Libraries used in the entire script

library(stringr)
library(stringi)
library(readxl)
library(tidyverse)
library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(plyr)
library(vegan)
library(micropan)


# Figure 2A, Figure 2B all strains are included

# Fig 2C and Fig S2A were produced out of Anvi'o

#setwd("~/Box Sync/Disha/Projects/carb_utilization/cinnocuum/roary_output/")

mydata = read.table("number_of_genes_in_pan_genome.Rtab")

data1 <- mydata
colnames(data1) <- 1:110
data2 <- pivot_longer(data1, cols = everything(), names_to = "no_of_genomes", values_to = "no_of_genes")
data2$no_of_genomes <- as.numeric(data2$no_of_genomes)

### calculating 
model <- lm(log(data2$no_of_genes) ~ log(data2$no_of_genomes))
summary(model)

## Result:
##Call:
##  lm(formula = log(data2$no_of_genes) ~ log(data2$no_of_genomes))

##Residuals:
##  Min       1Q   Median       3Q      Max 
##-0.44060 -0.03448 -0.00211  0.03946  0.31239 

##Coefficients:
##  Estimate Std. Error t value Pr(>|t|)    
##(Intercept)              8.567284   0.010706   800.2   <2e-16 ***
##  log(data2$no_of_genomes) 0.347141   0.002785   124.6   <2e-16 ***
##  ---
##  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##Residual standard error: 0.08573 on 1098 degrees of freedom
##Multiple R-squared:  0.934,	Adjusted R-squared:  0.9339 
##F-statistic: 1.553e+04 on 1 and 1098 DF,  p-value: < 2.2e-16


gene1 <- read.table("./gene_presence_absence.Rtab")

gene2 <- gene1
colnames(gene2) <- gene1[1, ]
rownames(gene2) <- gene1$V1
gene3 <- gene2[-1, ]
gene4 <- gene3[ , -1]
gene5 <- t(gene4)
set.seed(1234)
h.est <- heaps(gene5, n.perm = 500) 

## Result of hest
## Named num [1:2] 3105.488 0.818

### plotting Figure 2A; contains all strains

ggplot(data2, aes(x=factor(no_of_genomes), y=no_of_genes)) + 
  #geom_point() +
  geom_boxplot(outlier.shape = NA, fill= "#E8A419") +
  geom_smooth(method = 'nls', formula = 'y~a*x^b') +
  # geom_text(x = 600, y = 1, label = power_eqn(DD), parse = TRUE) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(colour = "black", size=1), text = element_text(size = 1)) +
  scale_x_discrete(breaks = seq(1, length(data2$no_of_genomes), by =2)) +
  geom_text(x = 15, y= 25000, label= "N = 5256.8328n^0.347141", size = 1) +
  geom_text(x= 30, y = 23000, label = "alpha = 0.818 (Heap's law), p-value < 2e-16", size = 1) +
  ggtitle("Total number of genes vs number of genomes")

## Fig 2B 

gene_presence_absence <- read_csv("~/Box Sync/Disha/Projects/carb_utilization/cinnocuum/roary_output/gene_presence_absence.csv")

colnames(gene_presence_absence)[colnames(gene_presence_absence) == "No. isolates"] <- "no_isolates" 
gene_presence_absence$`Avg sequences per isolate`[gene_presence_absence$`Avg sequences per isolate` == 1.00] <- 1
gene_presence_absence1 <- subset(gene_presence_absence, gene_presence_absence$`Avg sequences per isolate` < 1.02)

gene_count <- gene_presence_absence1 %>% group_by(no_isolates) %>% tally()
ggplot(gene_count, aes(x=no_isolates, y=n)) + geom_col(fill = "#E8A419") + theme_classic() +
  scale_x_continuous(breaks = seq(1, length(gene_count$no_isolates), by =2)) +
  theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(colour = "black", size=1), text = element_text(size = 1)) +
  ggtitle("Gene frequency vs number of genomes") +
  xlab("No. of genomes") +
  ylab("No. of genes") +
  geom_text(x = 100, y= 4000, label = "core", size =1)



## Fig S2A

clade12 = read.table("../clade12/clade_12_number_of_genes_in_pan_genome.Rtab")

clade12_2 <- clade12
colnames(clade12_2) <- 1:89
data_clade12 <- pivot_longer(clade12_2, cols = everything(), names_to = "no_of_genomes", values_to = "no_of_genes")
data_clade12$no_of_genomes <- as.numeric(data_clade12$no_of_genomes)


## calculating equation

model_clade12 <- lm(log(data_clade12$no_of_genes) ~ log(data_clade12$no_of_genomes))
summary(model_clade12)

## Result:
##Call:
##  lm(formula = log(data_clade12$no_of_genes) ~ log(data_clade12$no_of_genomes))

##Residuals:
##  Min        1Q    Median        3Q       Max 
##-0.237805 -0.025868 -0.002583  0.029113  0.169504 

##Coefficients:
##  Estimate Std. Error t value Pr(>|t|)    
##(Intercept)                     8.408014   0.006683    1258   <2e-16 ***
##  log(data_clade12$no_of_genomes) 0.308278   0.001835     168   <2e-16 ***
##  ---
##  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##Residual standard error: 0.05022 on 888 degrees of freedom
##Multiple R-squared:  0.9695,	Adjusted R-squared:  0.9695 
##F-statistic: 2.822e+04 on 1 and 888 DF,  p-value: < 2.2e-16


## calculating hest
gene12 <- read.table("../clade12/clade_12_gene_presence_absence.Rtab")

gene12_2 <- gene12
colnames(gene12_2) <- gene12[1, ]
rownames(gene12_2) <- gene12$V1
gene12_3 <- gene12_2[-1, ]
gene12_4 <- gene12_3[ , -1]
gene12_5 <- t(gene12_4)
set.seed(1234)
h.est <- heaps(gene12_5, n.perm = 500)

## Result:
## Named num [1:2] 1873.184 0.794

## plot FigS2B

ggplot(data_clade12, aes(x=factor(no_of_genomes), y=no_of_genes)) + 
  #geom_point() +
  geom_boxplot(outlier.shape = NA, fill= "#E8A419") +
  geom_smooth(method = 'nls', formula = 'y~a*x^b') +
  # geom_text(x = 600, y = 1, label = power_eqn(DD), parse = TRUE) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(colour = "black", size=1), text = element_text(size = 1)) +
  scale_x_discrete(breaks = seq(1, length(data_clade12$no_of_genomes), by =2)) +
  geom_text(x = 15, y= 17000, label= "N = 4482.8487n^0.308278", size = 1) +
  geom_text(x= 30, y = 16000, label = "alpha = 0.794 (Heap's law), p-value < 2e-16", size =1) +
  ggtitle("Total number of genes vs number of genomes")

# fig S2C

gene_presence_absence_12 <- read_csv("../clade12/clade_12_gene_presence_absence.csv")
colnames(gene_presence_absence_12)[colnames(gene_presence_absence_12) == "No. isolates"] <- "no_isolates" 

gene_count_1 <- subset(gene_presence_absence_12, gene_presence_absence_12$`Avg sequences per isolate` == 1)
gene_count12 <- gene_count_1 %>% group_by(no_isolates) %>% tally()
ggplot(gene_count12, aes(x=no_isolates, y=n)) + geom_col(fill = "#E8A419") + theme_classic() +
  scale_x_continuous(breaks = seq(1, length(gene_count12$no_isolates), by =2)) +
  theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(colour = "black", size=1), text = element_text(size = 1)) +
  ggtitle("Gene frequency vs number of genomes") +
  xlab("No. of genomes") +
  ylab("No. of unique genes") +
  geom_text(x = 80, y= 4000, label = "core", size = 1)




