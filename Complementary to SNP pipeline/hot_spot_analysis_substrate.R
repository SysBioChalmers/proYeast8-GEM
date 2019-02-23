# this code is used for handling with the hot spot analysis result for strain classification based on substrate

library(readxl)
library(tidyverse)
library(fitdistrplus)
library(logspline)
library(hongR)
library(superheat)

# input the data
hotspot_ex_high <- read.table("data/hotspot from pdb_ex for glycerol_high", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
hotspot_ex_high <- filter(hotspot_ex_high, str_detect(hotspot_ex_high$cluster,"cluster")==FALSE ) %>%
  filter(., pvalue <1)

hotspot_ex_medium <- read.table("data/hotspot from pdb_ex for glycerol_medium", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
hotspot_ex_medium <- filter(hotspot_ex_medium, str_detect(hotspot_ex_medium$cluster,"cluster")==FALSE ) %>%
  filter(., pvalue <1)

hotspot_ex_low <- read.table("data/hotspot from pdb_ex for glycerol_low", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
hotspot_ex_low <- filter(hotspot_ex_low, str_detect(hotspot_ex_low$cluster,"cluster")==FALSE ) %>%
  filter(., pvalue <1)



# homo
hotspot_homo_high <- read.table("data/hotspot from pdb_homo for glycerol_high", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
hotspot_homo_high <- filter(hotspot_homo_high, str_detect(hotspot_homo_high$cluster,"cluster")==FALSE ) %>%
  filter(., pvalue <1)

hotspot_homo_medium <- read.table("data/hotspot from pdb_homo for glycerol_medium", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
hotspot_homo_medium <- filter(hotspot_homo_medium, str_detect(hotspot_homo_medium$cluster,"cluster")==FALSE ) %>%
  filter(., pvalue <1)

hotspot_homo_low <- read.table("data/hotspot from pdb_homo for glycerol_low", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
hotspot_homo_low <- filter(hotspot_homo_low, str_detect(hotspot_homo_low$cluster,"cluster")==FALSE ) %>%
  filter(., pvalue <1)


# merging all the hotspot
# comparing in closeness distribution
hotspot <- rbind.data.frame(hotspot_ex_high,hotspot_homo_high,hotspot_ex_medium,hotspot_homo_medium,
                            hotspot_ex_low,hotspot_homo_low)
hotspot$closeness <- as.numeric(hotspot$closeness)
hotspot$pvalue <- as.numeric(hotspot$pvalue)


# input the gene annotation data
gene_annotation <- read.delim2('/Users/luho/Google Drive/R application and code/protein 3D structure QC and QA/SNP analysis pipeline/data/all_gene_yeast with annotation from different database.txt', header = TRUE, stringsAsFactors = FALSE)
kinetics_analysis <- read.table("/Users/luho/Documents/GitHub/ecModels/ecYeastGEM/figures/yeast8/results/KcatSensitivities_YEP.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(kinetics_analysis) <- c("Protein", "glucose", "acetate", "ethanol", "glycerol", "sorbitol", "galactose", "ribose", "xylose")

plot(density(kinetics_analysis$glycerol))


hotspot$kinectics <- getSingleReactionFormula(kinetics_analysis$glycerol,kinetics_analysis$Protein, hotspot$gene)
hotspot$kinectics <- as.numeric(hotspot$kinectics)

hotspot <- hotspot %>% filter(.,kinectics >= 0.000000001)



























# POST-analysis
# choose based on group
hotspot_high <- hotspot[hotspot$stain_type =="glycerol_high",]
hotspot_medium <- hotspot[hotspot$stain_type =="glycerol_medium",]
hotspot_low <- hotspot[hotspot$stain_type =="glycerol_low",]


high <- density(hotspot$closeness[hotspot$stain_type =="glycerol_high"])
medium <- density(hotspot$closeness[hotspot$stain_type =="glycerol_medium"])
low <- density(hotspot$closeness[hotspot$stain_type =="glycerol_low"])


plot(high,col = "steelblue",
     main = "Density of closeness in different strain stype",
     xlim = c(1,15),
     xlab="closeness",
     lty=1,
     lwd=3)
lines(medium, lty=1,lwd=3)
lines(low, lty=1,lwd=3)




# mapping based on gene pdb id
high0 <- hotspot_high %>% unite(.,id_mapping,c('gene','structure'), sep = "#")
medium0 <- hotspot_medium %>% unite(.,id_mapping,c('gene','structure'), sep = "#")
low0 <- hotspot_low %>% unite(.,id_mapping,c('gene','structure'), sep = "#")
g40 <- hotspot_g4 %>% unite(.,id_mapping,c('gene','structure'), sep = "#")
high0$cluster_medium <- getMultipleReactionFormula(medium0$cluster,medium0$id_mapping,high0$id_mapping)
high0$cluster_low <- getMultipleReactionFormula(low0$cluster,low0$id_mapping,high0$id_mapping)
high0$cluster_g4 <- getMultipleReactionFormula(g40$cluster,g40$id_mapping,high0$id_mapping)



# calculate the ratio of residues in hotspot of high0 which occured in other groups
e2 <- vector()
e3 <- vector()
e4 <- vector()
for (i in seq_along(high0$cluster)) {
  s1 <- high0$cluster[i]
  s2 <- high0$cluster_medium[i]
  s3 <- high0$cluster_low[i]
  s4 <- high0$cluster_g4[i]

  s10 <- unlist(str_split(s1, ";"))
  s20 <- unlist(str_split(s2, ";"))
  s30 <- unlist(str_split(s3, ";"))
  s40 <- unlist(str_split(s4, ";"))

  if (length(intersect(s10, s20)) >= 1) {
    e20 <- length(intersect(s10, s20))/length(s10)
  } else {
    e20 <- 0
  }

  if (length(intersect(s10, s30)) >= 1) {
    e30 <-length(intersect(s10, s30))/length(s10)
  } else {
    e30 <- 0
  }

  if (length(intersect(s10, s40)) >= 1) {
    e40 <- length(intersect(s10, s40))/length(s10)
  } else {
    e40 <- 0
  }

  e2[i] <- e20
  e3[i] <- e30
  e4[i] <- e40
}

high0$in_high <- 1
high0$in_medium <- e2
high0$in_low <- e3
high0$in_g4 <- e4


# question: how to filter the redundent hoptspot
high1 <- filter(high0,closeness >= 6)


# draw the heatmap
rownames(high1) <- paste("hotspot",1:length(high1$in_high),sep = "")
high0_heat <- high1[,c('in_high', 'in_medium','in_low','in_g4')]
colnames(high0_heat) <- c('high','medium','low','g4')

superheat(high0_heat,
          # scale the matrix columns
          scale = FALSE,
          # add row dendrogram
          row.dendrogram = TRUE,
          left.label.size = 0.4,
          bottom.label.size = 0.05,
          left.label.text.size = 4,
          bottom.label.text.size = 5,
          grid.hline.col = "white",
          grid.vline.col = "white")


write.table(high0, "result/high0 cluster mapping.txt", row.names = FALSE, sep = "\t")
write.table(high1, "result/high0 cluster mapping filtering.txt", row.names = TRUE, sep = "\t")

