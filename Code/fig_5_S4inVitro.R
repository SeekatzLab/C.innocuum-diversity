# Figure 5 

library(growthcurver)
library(readxl)
library(dplyr)
library(statmod)
library(ggplot2)

# Read in data - growth_curve_data.xlsx

# Plot growth curves for all sugars and all strains wrapped by sugar type
ggplot(growth_curve_data, aes(x=time, y=mean, color=species)) +
  geom_line(size= 1) + 
  geom_point(size=.1) +
  #geom_smooth() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.001, size = 0.15) +
  facet_wrap(~sugar_type) +
  theme_classic() +
  theme(axis.text.x = element_text(angle =90)) +
  labs(x = "Time (min)", y = "OD600", title = '') +
  scale_color_manual(values=c("#74a365", "#E88219", "#9FC095","#F5191C","#caddc5","#f86264","#505050","#E8A419"))

# Use growthcurver package to get output metrics for all sugars, including AUC, for statistical analysis 

sugar_x <- read_tsv("~/sugar_x.txt")
gcvr_output <- SummarizeGrowthByPlate(sugar)
gcvr_output_metrics <- "~/gcvr_output_metrics.txt"
write.table(gcvr_output, file = gcvr_output_metrics, 
            quote = FALSE, sep = "\t", row.names = FALSE)


# Heatmap of acid production and growth curve data
# Read in data with statistical growth significance (based on AUC ANOVA and Tukey HSD stats), growth_curve_data, sheet = 'heatmap'

ggplot(heatmap, aes(x=strain, y=sugar_type, fill=as.character(growth_acid))) + 
  geom_raster() +
  scale_fill_manual(values = c("snow4","darkorchid3","darkgoldenrod2", "gold"), name = '', guide = guide_legend(reverse = TRUE), labels = c("No growth", "Growth, no acid producion","Growth, slight acid production", "Growth & high acid production")) + 
  theme_classic() +
  labs(title = '', x = 'Bacterial Strain', y = 'Sugar Type') +
  scale_x_discrete(limits = c('CM151C', 'CM152','CM208A', 'CM220', 'CM647', 'CM679','d22_479')) +
  scale_y_discrete(limits = c('salicin','raffinose','lactose','sorbitol','trehalose','maltose','cellobiose','mannitol','sucrose','mannose','fructose','glucose')) +
  theme(legend.key.size = unit(1, 'cm'), legend.key.height = unit(2, 'cm'))


