library(readxl)
library(tidyverse)
library(fitdistrplus)
library(logspline)
library(hongR)
#------part 3 strain phenotype analysis based on the growth rate
#temperature analysis
pheno <- read_excel("data/phenoMatrix_35ConditionsNormalizedByYPD.xlsx")
conditionList <- colnames(pheno)
conditionList <- conditionList[-1]

for (i in conditionList){
  print(i)
  ggplot(pheno, aes_string(i)) +
  geom_density(fill="lightblue") +
  xlab(i) + 
  ylab("Relative growth rate")
  ggsave(out <- paste('result/',i,'.png', sep = ""), width=4, height=4, dpi=300)

}

# distribution under glycerol+YP
j = 'YPGLYCEROL'
ggplot(pheno, aes_string(j)) +
  geom_density(fill="lightblue") +
  xlab("Relative growth rate") + 
  ylab("Density") +
  theme_bw() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=22))

#classification of strain based on substrate usage
#acetate
strain_acetate <- pheno[,c('strain_name','YPACETATE')]
mean(strain_acetate$YPACETATE)
strain_acetate0 <- filter(strain_acetate, YPACETATE>=0.60)
strain_acetate1 <- filter(strain_acetate, YPACETATE <= 0.17)
strain_acetate2 <- filter(strain_acetate, 0.40<=YPACETATE & YPACETATE <= 0.412)

strain_acetate0$type <- 'acetate_high'
strain_acetate1$type <- 'acetate_low'
strain_acetate2$type <- 'acetate_medium'

strain_acetate_classification <- rbind.data.frame(strain_acetate0, strain_acetate1,strain_acetate2)


#xylose
strain_xylose <- pheno[,c('strain_name','YPXYLOSE')]
mean(strain_xylose$YPXYLOSE)
strain_xylose0 <- filter(strain_xylose, YPXYLOSE>=0.60)
strain_xylose1 <- filter(strain_xylose, YPXYLOSE <= 0.17)
strain_xylose2 <- filter(strain_xylose, 0.394 <=YPXYLOSE & YPXYLOSE <= 0.41)

strain_xylose0$type <- 'xylose_high'
strain_xylose1$type <- 'xylose_low'
strain_xylose2$type <- 'xylose_medium'
strain_xylose_classification <- rbind.data.frame(strain_xylose0, strain_xylose1, strain_xylose2)


#sorbitol
strain_sorbitol <- pheno[,c('strain_name','YPSORBITOL')]
mean(strain_sorbitol$YPSORBITOL)
strain_sorbitol0 <- filter(strain_sorbitol, YPSORBITOL>=0.68)
strain_sorbitol1 <- filter(strain_sorbitol, YPSORBITOL <= 0.21)
strain_sorbitol2 <- filter(strain_sorbitol, 0.42<= YPSORBITOL & YPSORBITOL <= 0.435)

strain_sorbitol0$type <- 'sorbitol_high'
strain_sorbitol1$type <- 'sorbitol_low'
strain_sorbitol2$type <- 'sorbitol_medium'
strain_sorbitol_classification <- rbind.data.frame(strain_sorbitol0, strain_sorbitol1,strain_sorbitol2)

#glycerol
strain_glycerol <- pheno[,c('strain_name','YPGLYCEROL')]
mean(strain_glycerol$YPGLYCEROL)

strain_glycerol0 <- filter(strain_glycerol, YPGLYCEROL>=0.80)
strain_glycerol1 <- filter(strain_glycerol, YPGLYCEROL <= 0.25)
strain_glycerol2 <- filter(strain_glycerol, 0.515<= YPGLYCEROL & YPGLYCEROL <= 0.53)

strain_glycerol0$type <- 'high'
strain_glycerol1$type <- 'low'
strain_glycerol2$type <- 'medium'

strain_glycerol_classification <- rbind.data.frame(strain_glycerol0 ,strain_glycerol2, strain_glycerol1)
strain_glycerol_classification$type <- as.factor(strain_glycerol_classification$type)
strain_glycerol_classification %>%
  ggplot(aes(x=type,y=YPGLYCEROL, fill=type))+
  geom_boxplot(outlier.size=0) +  scale_x_discrete(limits = c('high','medium','low')) +
  geom_jitter(position=position_jitter(h=.1), aes(color = type), alpha=0.3) +
  xlab('Strain group') + ylab( 'Relative growth rate') +
  theme_bw() + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=22)) +
  theme(legend.position = "none")




#ethanol
strain_ethanol <- pheno[,c('strain_name','YPETHANOL')]
mean(strain_ethanol$YPETHANOL)

strain_ethanol0 <- filter(strain_ethanol, YPETHANOL>=0.9)
strain_ethanol1 <- filter(strain_ethanol, YPETHANOL <= 0.3)
strain_ethanol2 <- filter(strain_ethanol, 0.515<= YPETHANOL & YPETHANOL <= 0.53)

strain_ethanol0$type <- 'ethanol_high'
strain_ethanol1$type <- 'ethanol_low'
strain_ethanol2$type <- 'ethanol_medium'
strain_ethanol_classification <- rbind.data.frame(strain_ethanol0, strain_ethanol1,strain_ethanol2)

#ribose
strain_ribose <- pheno[,c('strain_name','YPRIBOSE')]
mean(strain_ribose$YPRIBOSE)
summary(strain_ribose$YPRIBOSE)
strain_ribose0 <- filter(strain_ribose, YPRIBOSE >=0.62)
strain_ribose1 <- filter(strain_ribose, YPRIBOSE <= 0.17)
strain_ribose2 <- filter(strain_ribose, 0.40<= YPRIBOSE & YPRIBOSE <= 0.42)

strain_ribose0$type <- 'ribose_high'
strain_ribose1$type <- 'ribose_low'
strain_ribose2$type <- 'ribose_medium'
strain_ribose_classification <- rbind.data.frame(strain_ribose0, strain_ribose1,strain_ribose2)



#galactose
strain_galactose <- pheno[,c('strain_name','YPGALACTOSE')]
mean(strain_galactose$YPGALACTOSE)
summary(strain_galactose$YPGALACTOSE)
plot(density(strain_galactose$YPGALACTOSE))
strain_galactose0 <- filter(strain_galactose, YPGALACTOSE >=1.37)
strain_galactose1 <- filter(strain_galactose, YPGALACTOSE <= 0.4)
strain_galactose2 <- filter(strain_galactose, 0.89 <= YPGALACTOSE & YPGALACTOSE <= 0.93)

strain_galactose0$type <- 'galactose_high'
strain_galactose1$type <- 'galactose_low'
strain_galactose2$type <- 'galactose_medium'
strain_galactose_classification <- rbind.data.frame(strain_galactose0, strain_galactose1,strain_galactose2)


#glucose
strain_glucose <- pheno[,c('strain_name','YPD40')] #here we choose the relative growth at 40 degree
plot(density(strain_glucose$YPD40))
mean(strain_glucose$YPD40)
summary(strain_glucose$YPD40)
strain_glucose0 <- filter(strain_glucose, YPD40 >=1)
strain_glucose1 <- filter(strain_glucose, YPD40 <= 0.28)
strain_glucose2 <- filter(strain_glucose, 0.61 <= YPD40 & YPD40 <= 0.633)

strain_glucose0$type <- 'glucose_high'
strain_glucose1$type <- 'glucose_low'
strain_glucose2$type <- 'glucose_medium'
strain_glucose_classification <- rbind.data.frame(strain_glucose0, strain_glucose1,strain_glucose2)


# save the result
write.table(strain_xylose_classification, "result/strain_xylose_classification.txt", row.names = FALSE, sep = "\t")
write.table(strain_glycerol_classification, "result/strain_glycerol_classification.txt", row.names = FALSE, sep = "\t")
write.table(strain_sorbitol_classification, "result/strain_sorbitol_classification.txt", row.names = FALSE, sep = "\t")
write.table(strain_acetate_classification, "result/strain_acetate_classification.txt", row.names = FALSE, sep = "\t")
write.table(strain_ethanol_classification, "result/strain_ethanol_classification.txt", row.names = FALSE, sep = "\t")
write.table(strain_ribose_classification, "result/strain_ribose_classification.txt", row.names = FALSE, sep = "\t")
write.table(strain_galactose_classification, "result/strain_galactose_classification.txt", row.names = FALSE, sep = "\t")
write.table(strain_glucose_classification, "result/strain_glucose_classification.txt", row.names = FALSE, sep = "\t")



  