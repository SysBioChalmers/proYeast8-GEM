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
hotspot0 <- hotspot %>% filter(.,kinectics >= 0.000000001)

# initial relation between the closenes and kinetics
plot(hotspot$closeness, hotspot$kinectics)



# POST-analysis
# choose based on group
hotspot_high <- hotspot0[hotspot0$stain_type =="glycerol_high",]
hotspot_medium <- hotspot0[hotspot0$stain_type =="glycerol_medium",]
hotspot_low <- hotspot0[hotspot0$stain_type =="glycerol_low",]


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
high0 <- hotspot_high %>% unite(.,id_mapping,c('gene','structure'), sep = ";")
medium0 <- hotspot_medium %>% unite(.,id_mapping,c('gene','structure'), sep = ";")
low0 <- hotspot_low %>% unite(.,id_mapping,c('gene','structure'), sep = ";")

cluster_sum <- rbind.data.frame(high0, medium0, low0)
cluster_sum$cluster_high <- getMultipleReactionFormula(high0$cluster,high0$id_mapping,cluster_sum$id_mapping)
cluster_sum$cluster_medium <- getMultipleReactionFormula(medium0$cluster,medium0$id_mapping,cluster_sum$id_mapping)
cluster_sum$cluster_low <- getMultipleReactionFormula(low0$cluster,low0$id_mapping,cluster_sum$id_mapping)




# calculate the ratio of residues in hotspot of high0 which occured in other groups
e1 <- vector()
e2 <- vector()
e3 <- vector()

for (i in seq_along(cluster_sum$cluster)) {
  s1 <- cluster_sum$cluster[i]
  s2 <- cluster_sum$cluster_high[i]
  s3 <- cluster_sum$cluster_medium[i]
  s4 <- cluster_sum$cluster_low[i]
  

  s10 <- unique(unlist(str_split(s1, ";")))
  s20 <- unique(unlist(str_split(s2, ";")))
  s30 <- unique(unlist(str_split(s3, ";")))
  s40 <- unique(unlist(str_split(s4, ";")))

  if (length(intersect(s10, s20)) >= 1) {
    e10 <- length(intersect(s10, s20))/length(s10)
  } else {
    e10 <- 0
  }

  if (length(intersect(s10, s30)) >= 1) {
    e20 <-length(intersect(s10, s30))/length(s10)
  } else {
    e20 <- 0
  }

  if (length(intersect(s10, s40)) >= 1) {
    e30 <-length(intersect(s10, s40))/length(s10)
  } else {
    e30 <- 0
  }
  
  e1[i] <- e10
  e2[i] <- e20
  e3[i] <- e30

}

cluster_sum$ratio_high <- e1
cluster_sum$ratio_medium <- e2
cluster_sum$ratio_low <- e3

#remove the duplicated based on the gene and the cluster
cluster_sum <- cluster_sum %>% separate(.,id_mapping, into = c('gene','structure'), sep = ";")
cluster_sum$unique_sign <- paste(cluster_sum$gene, cluster_sum$cluster, sep = "_")
cluster_sum1 <- cluster_sum[!duplicated(cluster_sum$unique_sign), ]


# question: how to filter the redundent hoptspot
cluster_sum1 <- filter(cluster_sum1,kinectics >= 0.001)


# draw the heatmap
cluster_sum1$unique_sign <- str_replace_all(cluster_sum1$unique_sign, "@@", "@") %>%
  str_replace_all(.,";","+")
rownames(cluster_sum1) <- cluster_sum1$unique_sign
cluster_sum2 <- cluster_sum1[,c('ratio_high', 'ratio_medium','ratio_low')]
colnames(cluster_sum2) <- c('high','medium','low')


superheat(cluster_sum2,
          # scale the matrix columns
          scale = FALSE,
          # add row dendrogram
          row.dendrogram = TRUE,
          left.label.size = 0.8,
          bottom.label.size = 0.05,
          left.label.text.size = 2,
          bottom.label.text.size = 3,
          grid.hline.col = "white",
          grid.vline.col = "white")




write.table(cluster_sum1, "result/hotspot analysis for glycerol.txt", row.names = FALSE, sep = "\t")


#plot the heat map for the fcc of the hotspot
hotspot_fcc <- read.table('data/hotspot_fcc.txt', sep = "\t", stringsAsFactors = FALSE)
colnames(hotspot_fcc) <- "hotspot"

hotspot_fcc1  <- hotspot_fcc  %>% separate(., hotspot, into = c('gene','mutation'), sep="_")
hotspot_fcc1$hotspot <- hotspot_fcc$hotspot
hotspot_fcc1$fcc <- getSingleReactionFormula(kinetics_analysis$glycerol,kinetics_analysis$Protein,hotspot_fcc1$gene)
write.table(hotspot_fcc1, "result/hotspot_fcc1.txt", row.names = FALSE, sep = "\t")
























