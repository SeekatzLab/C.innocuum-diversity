#Supplementary Figure for BCP data

library(growthcurver)
library(readxl)
library(dplyr)

# Read in data for BCP heatmap - BCP_supplemental.tsv
# HeatMap using actual pH values 
ggplot(BCP_supplemental, aes(x=species, y=sugar_type, fill=ph_after24)) + 
  geom_raster() + 
  scale_fill_gradient(low = "goldenrod1", high = "darkorchid4", name = 'No Growth') +
  geom_text(aes(species, sugar_type, label=round(ph_after24, digits = 3)),color = "white", check_overlap = TRUE) +
  theme_classic() +
  labs(title = 'Bromocresol Purple pH Values', x = 'Bacterial Strain', y = 'Sugar Type') +
  scale_x_discrete(limits = c('CM151C', 'CM152','CM208A', 'CM220', 'CM647', 'CM679','d22_479', 'Control')) +
  scale_y_discrete(limits = c('BMCA_none','lactose','salicin','raffinose','sorbitol','maltose','cellobiose','mannitol','sucrose','trehalose','mannose','fructose','glucose')) +
  theme(legend.key.size = unit(1, 'cm'), legend.key.height = unit(2, 'cm'))
# Heatmap using calculated pH values 
ggplot(BCP_supplemental, aes(x=species, y=sugar_type, fill=calc_pH)) + 
  geom_raster() + 
  scale_fill_gradient(low = "goldenrod1", high = "darkorchid4", name = 'No Growth') +
  geom_text(aes(species, sugar_type, label=round(calc_pH, digits = 2)),color = "white", check_overlap = TRUE) +
  theme_classic() +
  labs(title = '', x = 'Bacterial Strain', y = 'Sugar Type') +
  scale_x_discrete(limits = c('CM151C', 'CM152','CM208A', 'CM220', 'CM647', 'CM679','d22_479', 'Control')) +
  scale_y_discrete(limits = c('BMCA_none','lactose','salicin','raffinose','sorbitol','maltose','cellobiose','mannitol','sucrose','trehalose','mannose','fructose','glucose')) +
  theme(legend.key.size = unit(1, 'cm'), legend.key.height = unit(2, 'cm'))

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
#read in data for pH regression figure - pH_regression_data.xlsx, sheet = 'controls'

# Calculate regression lines for controls
#BMCA control @ 588
  
fit_control588 <- lm(abs_BMCA_588 ~ pH, data = pH_regression_data)
summary(fit_control588)
ggplot(pH_regression_data,aes(y=abs_BMCA_588,x=pH))+geom_point()+geom_smooth(method="lm")

ggplotRegression_control588 <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("BMCA_BCP_control","Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}
ggplotRegression(fit_control588)

# Use regression equation y=0.18x- 0.73 (where 'x' is pH) to calculate expected pH values

#TCCFB control @ 588

fit_control_TCC_588 <- lm(abs_TCCFB_588 ~ pH, data = pH_regression_data)
summary(fit_control_TCC_588)
ggplot(pH_regression_data,aes(y=abs_TCCFB_588,x=pH))+geom_point()+geom_smooth(method="lm")

ggplotRegression_control_TCC_588 <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("BMCA_BCP_control","Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}
ggplotRegression(fit_control_TCC_588)
  
#Plot regressions with calculated pH values for controls (TCCFB & BMCA)
  ggplot(pH_regression_data, aes(ph_after24, abs_588, shape=exp, colour=exp, fill=exp)) +
  geom_smooth(method="lm") +
  geom_point(size=2) +
  theme_classic() + 
  xlab("pH") +
  ylab("Absorbance OD(588)") + #absorbance is at 588
  ggtitle("") +
  scale_color_manual(breaks = c("f","g","j"), labels = c("TCCFB: y= -0.39 + 0.16x R2 = 0.73", "BMCA: y= -0.73 + 0.18x R2 = 0.94","Calculated pH = (abs+0.73)/0.18"), values = c("blue4", "skyblue","dodgerblue2")) +
  scale_fill_manual(breaks = c("f","g","j"), labels = c("TCCFB: y= -0.39 + 0.16x R2 = 0.73", "BMCA: y= -0.73 + 0.18x R2 = 0.94","Calculated pH = (abs+0.73)/0.18"), values = c("blue4", "skyblue","dodgerblue2")) +
  scale_shape_manual(breaks = c("f","g","j"), labels = c("TCCFB: y= -0.39 + 0.16x R2 = 0.73","BMCA: y= -0.73 + 0.18x R2 = 0.94","Calculated pH = (abs+0.73)/0.18"), values = c(16, 17, 18))





