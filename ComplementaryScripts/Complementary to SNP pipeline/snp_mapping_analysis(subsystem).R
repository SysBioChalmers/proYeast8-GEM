library(readxl)
library(tidyverse)
library(fitdistrplus)
library(logspline)
library(hongR)
library(hexbin)
#input the genome annotation of gene from SGD database
s288_sgd_annotation <- read_tsv('data/s288_genome.tsv')


# function to standardize the snp number for each gene
calc_relative_snp <- function(geneName, snp_number, annotation = s288_sgd_annotation) {
  #gene_list <- geneGEM$geneNames
  #snp_list <- geneGEM$SNP_NUM
  gene_list <- geneName
  snp_list <- snp_number
  snp_list <- as.numeric(snp_list)
  protein_length <- getSingleReactionFormula(annotation$protein_length, annotation$systematic_name, gene_list)
  protein_length <- as.numeric(protein_length)
  snp_relative <- snp_list / protein_length
  return(snp_relative)
}


# nsSNP analysis
geneGEM <- read_excel("data/geneGEM with SNP number.xlsx")
# standard the SNP number based on the protein lenght
geneGEM$SNP_NUM <- calc_relative_snp(geneGEM$locus_tag,geneGEM$SNP_NUM)
geneGEM$nsSNP <- calc_relative_snp(geneGEM$locus_tag,geneGEM$nsSNP)
geneGEM0 <- geneGEM[geneGEM$SNP_NUM > 0 & !is.na(geneGEM$SNP_NUM),]

x <- geneGEM0$SNP_NUM
dens1 <- density(x)

# give the subsystem information
sce_kegg <- read_excel("data/sce_kegg.xlsx", 
                       sheet = "gene_pathway")
pathwayList <- read_excel("data/sce_kegg.xlsx", 
                       sheet = "pathwayList")

geneGEM0$pathwayID <- getMultipleReactionFormula(sce_kegg$pathwayID,sce_kegg$geneID,geneGEM0$locus_tag)
geneGEM0 <- geneGEM0 %>% filter(.,!is.na(pathwayID))
gene_pathway <- splitAndCombine(geneGEM0$pathwayID, geneGEM0$locus_tag, sep0 = ";")
colnames(gene_pathway) <- c('pathwayID','locus_tag')


gene_pathway$pathwayName <- getSingleReactionFormula(pathwayList$pathwayName,pathwayList$pathwayID,gene_pathway$pathwayID)
gene_pathway$nsSNP <- getSingleReactionFormula(geneGEM0$nsSNP,geneGEM0$locus_tag,gene_pathway$locus_tag)
gene_pathway$nsSNP <- as.numeric(gene_pathway$nsSNP)

#calculate the average nsSNP value based on subsystem group
gene_pathway$pathwayName <- as.factor(gene_pathway$pathwayName)


nsSNP_subsytem <- gene_pathway %>%
  group_by(pathwayName) %>%
  dplyr::summarize(Mean = mean(nsSNP, na.rm=TRUE))

gene_pathway0 <- gene_pathway %>% filter(.,pathwayID=='sce01110')
ggplot(gene_pathway0, aes(x=nsSNP, y=..density..)) +
  geom_histogram(fill="red", colour="red", size=.2) +
  geom_density() + 
  ggtitle("biosynthesis of secondary metabolites")

for (i in seq_along(nsSNP_subsytem$pathwayName)){
  print(i)
  s <- nsSNP_subsytem$pathwayName[i]
  print(s)
  gene_pathway0 <- gene_pathway %>% filter(.,pathwayName==s) 
   ggplot(gene_pathway0, aes(x=nsSNP, y=..density..)) +
    geom_histogram(fill="red", colour="red", size=.2) +
    geom_density() + 
    ggtitle(s)
  s <- str_replace_all(s, "\\ / ", "_")
   ggsave(out <- paste('result/',s,'.png', sep = ""), width=4, height=4, dpi=300)
  
}



#correlation analysis between kinectics analysis result and nsSNP number
uniprotGeneID_mapping <- read_excel("data/uniprotGeneID_mapping.xlsx")
# input the result for the kinetics sensitivity analysis
kinetics_analysis <- read.table("/Users/luho/Documents/GitHub/ecModels/ecYeastGEM/figures/yeast8/results/KcatSensitivities_YEP.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(kinetics_analysis) <- c("Protein", "glucose", "acetate", "ethanol", "glycerol", "sorbitol", "galactose", "ribose", "xylose")
# id-mapping
kinetics_analysis$locus <- getSingleReactionFormula(uniprotGeneID_mapping$GeneName,uniprotGeneID_mapping$Entry,kinetics_analysis$Protein)


geneGEM0$kinetics <- getSingleReactionFormula(kinetics_analysis$xylose,kinetics_analysis$locus,geneGEM0$locus_tag)
geneGEM0$kinetics <- as.numeric(geneGEM0$kinetics)
geneGEM1 <- geneGEM0 %>% filter(.,kinetics >=0.000001050)
plot(geneGEM1$kinetics, geneGEM1$nsSNP)


sp <- ggplot(geneGEM1, aes(x=kinetics, y=nsSNP))
sp + geom_point()
sp + stat_binhex() +
  scale_fill_gradient(low="black", high="red",
                      breaks=c(0,2,4,6,8,10),
                      limits=c(0, 10))


#correlation analysis between protein abundence and nsSNP number
protein_abundence <- read.table('data/protein_abundence_paxdb.txt', header = TRUE, stringsAsFactors = FALSE)
protein_abundence$string_external_id <- str_replace_all(protein_abundence$string_external_id,'4932.','')

geneGEM$abundence <- getSingleReactionFormula(protein_abundence$abundance,protein_abundence$string_external_id,geneGEM$locus_tag)
geneGEM$abundence <- as.numeric(geneGEM$abundence)
sp <- ggplot(geneGEM, aes(x=abundence, y=nsSNP))
sp + geom_point()
sp + stat_binhex() +
  scale_fill_gradient(low="black", high="red",
                      breaks=c(0,2,4,6,8,10),
                      limits=c(0, 10))



#correlation analysis between rxn number connect with gene and nsSNP number
rxn_gene_mapping <- read_excel("data/rxn_gene_mapping.xlsx")
rxn_gene_mapping$rxn_gene <- paste(rxn_gene_mapping$rxnID, rxn_gene_mapping$gene, sep = ";")
rxn_gene_mapping0 <- rxn_gene_mapping[!duplicated(rxn_gene_mapping$rxn_gene),]
gene_rxn_number <- as.data.frame(table(rxn_gene_mapping0$gene))


geneGEM$rxn_number <- getSingleReactionFormula(gene_rxn_number$Freq,gene_rxn_number$Var1,geneGEM$locus_tag)
geneGEM$rxn_number <- as.numeric(geneGEM$rxn_number) 
geneGEM2<-  geneGEM %>% filter(.,!is.na(rxn_number))
sp <- ggplot(geneGEM2, aes(x=rxn_number, y=nsSNP))
sp + geom_point()
sp + stat_binhex() +
  scale_fill_gradient(low="black", high="red",
                      breaks=c(0,2,4,6,8,10),
                      limits=c(0, 10))

















