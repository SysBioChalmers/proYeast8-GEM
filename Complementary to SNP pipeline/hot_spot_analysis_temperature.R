# this code is used for handling with the hot spot analysis result for strain classification based on temperature

library(readxl)
library(tidyverse)
library(fitdistrplus)
library(logspline)
library(hongR)
library(superheat)

# input the data
hotspot_ex_g1 <- read.table("data/hotspot_ex_g1.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
hotspot_ex_g1 <- filter(hotspot_ex_g1, str_detect(hotspot_ex_g1$cluster,"cluster")==FALSE ) %>%
  filter(., pvalue <1)

hotspot_ex_g2 <- read.table("data/hotspot_ex_g2.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
hotspot_ex_g2 <- filter(hotspot_ex_g2, str_detect(hotspot_ex_g2$cluster,"cluster")==FALSE ) %>%
  filter(., pvalue <1)

hotspot_ex_g3 <- read.table("data/hotspot_ex_g3.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
hotspot_ex_g3 <- filter(hotspot_ex_g3, str_detect(hotspot_ex_g3$cluster,"cluster")==FALSE ) %>%
  filter(., pvalue <1)

hotspot_ex_g4 <- read.table("data/hotspot_ex_g4.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
hotspot_ex_g4 <- filter(hotspot_ex_g4, str_detect(hotspot_ex_g4$cluster,"cluster")==FALSE ) %>%
  filter(., pvalue <1)


# homo
hotspot_homo_g1 <- read.table("data/hotspot_homo_g1.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
hotspot_homo_g1 <- filter(hotspot_homo_g1, str_detect(hotspot_homo_g1$cluster,"cluster")==FALSE ) %>%
  filter(., pvalue <1)

hotspot_homo_g2 <- read.table("data/hotspot_homo_g2.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
hotspot_homo_g2 <- filter(hotspot_homo_g2, str_detect(hotspot_homo_g2$cluster,"cluster")==FALSE ) %>%
  filter(., pvalue <1)

hotspot_homo_g3 <- read.table("data/hotspot_homo_g3.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
hotspot_homo_g3 <- filter(hotspot_homo_g3, str_detect(hotspot_homo_g3$cluster,"cluster")==FALSE ) %>%
  filter(., pvalue <1)

hotspot_homo_g4 <- read.table("data/hotspot_homo_g4.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
hotspot_homo_g4 <- filter(hotspot_homo_g4, str_detect(hotspot_homo_g4$cluster,"cluster")==FALSE ) %>%
  filter(., pvalue <1)


# merging all the hotspot
# comparing in closeness distribution
hotspot <- rbind.data.frame(hotspot_ex_g1,hotspot_homo_g1,hotspot_ex_g2,hotspot_homo_g2,
                            hotspot_ex_g3,hotspot_homo_g3,hotspot_ex_g4,hotspot_homo_g4)
hotspot$closeness <- as.numeric(hotspot$closeness)
hotspot$pvalue <- as.numeric(hotspot$pvalue)

hotspot[hotspot$gene == "YLR056W", ]

# choose based on group
hotspot_g1 <- hotspot[hotspot$stain_type =="g1",]
hotspot_g2 <- hotspot[hotspot$stain_type =="g2",]
hotspot_g3 <- hotspot[hotspot$stain_type =="g3",]
hotspot_g4 <- hotspot[hotspot$stain_type =="g4",]


g1 <- density(hotspot$closeness[hotspot$stain_type =="g1"])
g2 <- density(hotspot$closeness[hotspot$stain_type =="g2"])
g3 <- density(hotspot$closeness[hotspot$stain_type =="g3"])
g4 <- density(hotspot$closeness[hotspot$stain_type =="g4"])

plot(g1,col = "steelblue",
     main = "Density of closeness in different strain stype",
     xlim = c(1,15),
     xlab="closeness",
     lty=1,
     lwd=3)
lines(g2, lty=1,lwd=3)
lines(g3, lty=1,lwd=3)
lines(g4, lty=1,lwd=3)

# mapping based on gene pdb id
g10 <- hotspot_g1 %>% unite(.,id_mapping,c('gene','structure'), sep = "#")
g20 <- hotspot_g2 %>% unite(.,id_mapping,c('gene','structure'), sep = "#")
g30 <- hotspot_g3 %>% unite(.,id_mapping,c('gene','structure'), sep = "#")
g40 <- hotspot_g4 %>% unite(.,id_mapping,c('gene','structure'), sep = "#")
g10$cluster_g2 <- getMultipleReactionFormula(g20$cluster,g20$id_mapping,g10$id_mapping)
g10$cluster_g3 <- getMultipleReactionFormula(g30$cluster,g30$id_mapping,g10$id_mapping)
g10$cluster_g4 <- getMultipleReactionFormula(g40$cluster,g40$id_mapping,g10$id_mapping)



# calculate the ratio of residues in hotspot of g10 which occured in other groups
e2 <- vector()
e3 <- vector()
e4 <- vector()
for (i in seq_along(g10$cluster)) {
  s1 <- g10$cluster[i]
  s2 <- g10$cluster_g2[i]
  s3 <- g10$cluster_g3[i]
  s4 <- g10$cluster_g4[i]

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

g10$in_g1 <- 1
g10$in_g2 <- e2
g10$in_g3 <- e3
g10$in_g4 <- e4


# question: how to filter the redundent hoptspot
g11 <- filter(g10,closeness >= 6)


# draw the heatmap
rownames(g11) <- paste("hotspot",1:length(g11$in_g1),sep = "")
g10_heat <- g11[,c('in_g1', 'in_g2','in_g3','in_g4')]
colnames(g10_heat) <- c('g1','g2','g3','g4')

superheat(g10_heat,
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


write.table(g10, "result/g10 cluster mapping.txt", row.names = FALSE, sep = "\t")
write.table(g11, "result/g10 cluster mapping filtering.txt", row.names = TRUE, sep = "\t")

