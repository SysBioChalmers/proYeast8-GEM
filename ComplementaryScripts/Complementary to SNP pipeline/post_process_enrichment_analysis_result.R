library(readxl)
library(tidyverse)
library(hongR)
#------part 1 enrichment analysis of gene from genome with least and largest SNP number
gene_geneome_leastSNP <- read_excel("data/gene enrichment analysis in geneome level.xlsx",  sheet = "gene with least SNP")
gene_geneome_largestSNP <- read_excel("data/gene enrichment analysis in geneome level.xlsx",  sheet = "gene with largest SNP")

#filter one p_value <= 0.05
gene_geneome_leastSNP$PValue <- as.numeric(gene_geneome_leastSNP$PValue)
gene_geneome_largestSNP$PValue <- as.numeric(gene_geneome_largestSNP$PValue)

gene_geneome_leastSNP0 <- filter(gene_geneome_leastSNP, PValue <= 0.05) %>%
  filter(., Category == 'KEGG_PATHWAY'| Category == 'GOTERM_BP_FAT')#filter two: choose the function related to bioprocess

gene_geneome_largestSNP0 <- filter(gene_geneome_largestSNP, PValue <= 0.05) %>%
  filter(., Category == 'KEGG_PATHWAY'| Category == 'GOTERM_BP_FAT')#filter two: choose the function related to bioprocess

#enrichment analysis of gene from metabolic gene with least and largest SNP number
gene_gem_leastSNP <- read_excel("data/gene enrichment analysis in GEM level.xlsx", sheet = "GENE_lessSNP")
gene_gem_largestSNP <- read_excel("data/gene enrichment analysis in GEM level.xlsx", sheet = "GENE_moreSNP")

#filter one p_value <= 0.05
gene_gem_leastSNP$PValue <- as.numeric(gene_gem_leastSNP$PValue)
gene_gem_largestSNP$PValue <- as.numeric(gene_gem_largestSNP$PValue)

gene_gem_leastSNP0 <- filter(gene_gem_leastSNP, PValue <= 0.05) %>%
  filter(., Category == 'KEGG_PATHWAY'| Category == 'GOTERM_BP_DIRECT')#filter two: choose the function related to bioprocess

gene_gem_largestSNP0 <- filter(gene_gem_largestSNP, PValue <= 0.05) %>%
  filter(., Category == 'KEGG_PATHWAY'| Category == 'GOTERM_BP_DIRECT')#filter two: choose the function related to bioprocess








#------ part 2 enrichment analysis of gene in wine, ethonal and wild
gene_wine <- read_excel("data/enrichment analysis of wine clumps.xlsx", sheet = "gene_wine")
gene_bioethonal <- read_excel("data/enrichment analysis of wine clumps.xlsx", 
           sheet = "gene_bioethonal")
gene_wild <- read_excel("data/enrichment analysis of wine clumps.xlsx", 
           sheet = "gene_wild")

intersect(gene_wild$gene, gene_wine$gene)
intersect(gene_wild$gene, gene_bioethonal$gene)
intersect(gene_wine$gene, gene_bioethonal$gene)


# plot the result for the analysis from wine, bioethanol and wild type sce strains
df= data.frame(strain_type=c('wild','wine','bioethonal'), protein_num=c(65,41,22))
ggplot(df, aes(x=strain_type, y=protein_num, fill=strain_type)) +
  geom_bar(stat="identity") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  labs(x ="Strain type", y = "Protein number") +
  geom_text(aes(label=protein_num), vjust=-0.5) +
  theme(panel.background = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  coord_cartesian(ylim = c(0, 80))


df2= data.frame(strain_type=c('wild','wine','bioethonal'), strain_num=c(54,390,35))
ggplot(df2, aes(x=strain_type, y=strain_num, fill=strain_type)) +
  geom_bar(stat="identity") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  labs(x ="Strain type", y = "Strain number") +
  geom_text(aes(label=strain_num), vjust=-0.5) +
  theme(panel.background = element_blank()) +
  theme(axis.line = element_line(colour = "black"))





go_wine <- read_excel("data/enrichment analysis of wine clumps.xlsx", sheet = "Wine")
go_bioethonal <- read_excel("data/enrichment analysis of wine clumps.xlsx", sheet = "bioethonal")
go_wild <- read_excel("data/enrichment analysis of wine clumps.xlsx", sheet = "Wild")

signficantGO <- function(go_analysis){
  go_analysis$PValue <- as.numeric(go_analysis$PValue)
  go_analysis0 <- filter(go_analysis, PValue <= 0.05) %>%
  filter(., Category == 'KEGG_PATHWAY'| Category == 'GOTERM_BP_FAT')#filter two: choose the function related to bioprocess
  go_term <- unique(go_analysis0$Term)
  return(go_term)
}

go_term_wine <- signficantGO(go_wine)
go_term_bioethonal <- signficantGO(go_bioethonal)
go_term_wild <- signficantGO(go_wild)
intersect(go_term_wine, go_term_wild)
wine_vip <-setdiff(go_term_wine,go_term_wild)

#obtain genes for the different go term from wine
go_wine_different <- go_wine[go_wine$Term %in% wine_vip,]








#------part 3 enrichment analysis of gene from genome with least and largest relative SNP number

####compare the result from metabolic and the whole genome
genome_lessSNP <- read_excel("data/gene enrichment analysis in geneome level (relative SNP number).xlsx", 
                             sheet = "less_SNP")
####compare the result from metabolic and the whole genome
genome_moreSNP <- read_excel("data/gene enrichment analysis in geneome level (relative SNP number).xlsx", 
                             sheet = "more_SNP")
######compare the result from metabolic and the whole genome
genome_less_snSNP <- read_excel("data/gene enrichment analysis in geneome level (relative SNP number).xlsx", 
                             sheet = "less_snSNP")
######compare the result from metabolic and the whole genome
genome_more_snSNP <- read_excel("data/gene enrichment analysis in geneome level (relative SNP number).xlsx", 
                                sheet = "more_snSNP")



filterGeneral <- function(enrichment0){
  enrichment0$PValue <- as.numeric(enrichment0$PValue)
  s1 <- filter(enrichment0, PValue <= 0.05) %>%
    filter(., Category == 'GOTERM_BP_FAT')#filter two: choose the function related to bioprocess
  return(s1)
}

genome_lessSNP0 <- filterGeneral(genome_lessSNP)
genome_moreSNP0 <- filterGeneral(genome_moreSNP)

genome_less_snSNP0 <- filterGeneral(genome_less_snSNP)
genome_more_snSNP0 <- filterGeneral(genome_more_snSNP)

paste0(genome_lessSNP0$Term,collapse = "; ")
















