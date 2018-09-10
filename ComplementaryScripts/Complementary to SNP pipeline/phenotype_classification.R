library(readxl)
library(tidyverse)
library(fitdistrplus)
library(logspline)
library(hongR)
#------part 3 strain phenotype analysis based on the growth rate
#temperature analysis
pheno <- read_excel("data/phenoMatrix_35ConditionsNormalizedByYPD.xlsx")
df40 <- pheno$YPD40
df42 <- pheno$YPD42


descdist(df40, discrete =  FALSE)


#a normal distribution could describe the growth rate
fit.norm1 <- fitdist(df40,'norm')
plot(fit.norm1)
mean0 = fit.norm1$estimate[1]
sd0 = fit.norm1$estimate[2]

#a normal distribution could describe the growth rate
fit.norm1 <- fitdist(df42,'norm')
plot(fit.norm1)
mean1 = fit.norm1$estimate[1]
sd1 = fit.norm1$estimate[2]

GetStrainFromPvalue <- function(q_filter) {
  # obtain the Quantile based on p value
  # q_filter = 0.25
  q42 <- qnorm(q_filter, mean = mean1, sd = sd1, lower.tail = FALSE, log.p = FALSE)
  print(q42)
  # obtain the strain which have a higher growth rate on 42 and 40 degree
  sample2 <- filter(pheno, YPD42 >= q42)
  StrainAdaptHighTemp0 <- sample2$strain_name
  StrainAdaptHighTemp <- pheno[pheno$strain_name %in% StrainAdaptHighTemp0, ]
  StrainAdaptHighTemp <- StrainAdaptHighTemp[, c("strain_name")]
  return(StrainAdaptHighTemp)
}

StrainAdaptHighTemp_25 <- GetStrainFromPvalue(0.25)
StrainAdaptHighTemp_50 <- GetStrainFromPvalue(0.5)
StrainAdaptHighTemp_75 <- GetStrainFromPvalue(0.75)
StrainAdaptHighTemp_100 <- GetStrainFromPvalue(1)
StrainAdaptHighTemp_25_50 <- setdiff(StrainAdaptHighTemp_50, StrainAdaptHighTemp_25)
StrainAdaptHighTemp_50_75 <- setdiff(StrainAdaptHighTemp_75, StrainAdaptHighTemp_50)
StrainAdaptHighTemp_75_100 <- setdiff(StrainAdaptHighTemp_100, StrainAdaptHighTemp_75)
StrainAdaptHighTemp_25$type <- "g1"
StrainAdaptHighTemp_25_50$type <- "g2"
StrainAdaptHighTemp_50_75$type <- "g3"
StrainAdaptHighTemp_75_100$type <- "g4"
StrainAdaptHighTemp_combine <- rbind.data.frame(StrainAdaptHighTemp_25,StrainAdaptHighTemp_25_50,
                                                StrainAdaptHighTemp_50_75, StrainAdaptHighTemp_75_100)
StrainAdaptHighTemp_combine$growth <- getSingleReactionFormula(pheno$YPD42,pheno$strain_name,StrainAdaptHighTemp_combine$strain_name)
StrainAdaptHighTemp_combine$growth <- as.numeric(StrainAdaptHighTemp_combine$growth)
StrainAdaptHighTemp_combine$type <- as.factor(StrainAdaptHighTemp_combine$type)
#calculate the average value based on 4 groups
p <- ggplot(StrainAdaptHighTemp_combine, aes(x=type, y=growth)) + 
  geom_boxplot() +
  labs(x="group",y="relative growth") + 
  theme(axis.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=16)) 
p
write.table(StrainAdaptHighTemp_combine, "result/strain_classification based on relative growth under high temperature.txt", row.names = FALSE, sep = "\t")










#high ethonal analysis
df_ethonal <- pheno$YPETHANOL
descdist(df_ethonal, discrete =  FALSE)
fit.norm1 <- fitdist(df_ethonal,'norm')
plot(fit.norm1)
mean0 = fit.norm1$estimate[1]
sd0 = fit.norm1$estimate[2]
q40 <- qnorm(q_filter, mean = mean0, sd = sd0, lower.tail = FALSE, log.p = FALSE)
StrainAdaptHighEthonal <- filter(pheno, YPETHANOL >=q40)
StrainAdaptHighEthonal$type <- 'High_ethonal'
StrainAdaptHighEthonal <- StrainAdaptHighEthonal[,c('strain_name', 'type')]


#high xylose
df_ethonal <- pheno$YPXYLOSE
descdist(df_ethonal, discrete =  FALSE)
fit.norm1 <- fitdist(df_ethonal,'norm')
plot(fit.norm1)
mean0 = fit.norm1$estimate[1]
sd0 = fit.norm1$estimate[2]
q40 <- qnorm(q_filter, mean = mean0, sd = sd0, lower.tail = FALSE, log.p = FALSE)
StrainAdaptHighXylose <- filter(pheno, YPXYLOSE >=q40)
StrainAdaptHighXylose$type <- 'High_xylose'
StrainAdaptHighXylose <- StrainAdaptHighXylose[,c('strain_name', 'type')]

#merge the above three source
Strain_interest_top25 <- rbind.data.frame(StrainAdaptHighTemp,StrainAdaptHighEthonal,StrainAdaptHighXylose)

#save the data
write.table(Strain_interest_top25, "result/Strain_interest_top25.txt", row.names = FALSE, sep = "\t")



#part 4 strain classification
strain_classification <- read_excel("data/strain_classification.xls")
strain_classification <- strain_classification[, c('Standardized_name', 'Ecological_origins', 'Clades')]

clade_assignment <- read_excel("data/strain_classification.xls", sheet = "Table S19")
intersect(clade_assignment$Clade, unique(strain_classification$Clades))
strain_classification$origin <- getSingleReactionFormula(clade_assignment$`Clade assignment`,clade_assignment$Clade,strain_classification$Clades)
write.table(strain_classification, "result/strain_classification with origin.txt", row.names = FALSE, sep = "\t")















  
  