library(readxl)
library(readr)
library(tidyverse)
library(stringr)
library(hongR)

gene_all <- read_excel("data/gene_list_yeastGEM.xlsx",  sheet = "Sheet1")
gene_all$geneNames <- str_trim(gene_all$geneNames, side = "both")

# pdb data from gangLI
yeast_3D_information <- read_excel("data/yeast_3D_swiss.xls",  sheet = "Protein")
yeast_3D_information$locus <- str_trim(yeast_3D_information$locus, side = "both")

# pdb data from uniprot
uniprot_yeast_protein_3D <- read_excel("data/yeast_3D_uniprot.xlsx", sheet = "Sheet0")
uniprot_yeast_protein_3D$`Gene names  (ordered locus )` <- str_trim(uniprot_yeast_protein_3D$`Gene names  (ordered locus )`, side = "both")
uniprot_yeast_protein_3D <- filter(uniprot_yeast_protein_3D, is.na(`Cross-reference (PDBsum)`)==FALSE)
uniprot_protein_3D <- splitAndCombine(uniprot_yeast_protein_3D$`Cross-reference (PDBsum)`,uniprot_yeast_protein_3D$`Gene names  (ordered locus )`, sep0 = ";")
colnames(uniprot_protein_3D) <- c("pdb_id", "gene")
uniprot_protein_3D <- filter(uniprot_protein_3D, pdb_id != "")

# pdb file statistics analysis
gene_all$pdb_id <- getMultipleReactionFormula(yeast_3D_information$pdb_id,yeast_3D_information$locus,gene_all$geneNames)
gene_all$pdb_number <- str_count(gene_all$pdb_id, ";") + 1

gene_all$pdb_id_uniprot <- getMultipleReactionFormula(uniprot_protein_3D$pdb_id,uniprot_protein_3D$gene,gene_all$geneNames)
gene_all$pdb_number_uniprot <- str_count(gene_all$pdb_id_uniprot, ";") +1

# find the gene without 3D in both swiss and sgd
protein_sequence <- read_tsv(file = 'data/yeast_protein_sequence_SGD.tsv')
gene_3D_swiss <- filter(gene_all, pdb_number >=1)
gene_no3D_swiss <- filter(gene_all, is.na(pdb_number) ==TRUE)
gene_3D_pdb <- filter(gene_all, pdb_number_uniprot >=1)
gene_no3D_pdb <- filter(gene_all, is.na(pdb_number_uniprot) ==TRUE)


gene_swiss <- gene_3D_swiss$geneNames
gene_pdb <- gene_3D_pdb$geneNames

gene_without_3D_swiss <- gene_no3D_swiss$geneNames
gene_without_3D_pdb <- gene_no3D_pdb$geneNames

gene_with_3D <- unique(c(gene_swiss,gene_pdb))
gene_without_3D <- intersect(gene_without_3D_swiss,gene_without_3D_pdb)

#obtain the protein sequence for these without 3D
index0 <- which(protein_sequence$`Gene Systematic name` %in% gene_without_3D ==TRUE)
protein_gene_without3D <- data.frame(geneName=gene_without_3D,stringsAsFactors = FALSE)
protein_gene_without3D$sequence <- getSingleReactionFormula(protein_sequence$`Sequence residue`,protein_sequence$`Gene Systematic name`,protein_gene_without3D$geneName)
protein_gene_without3D$length <- getSingleReactionFormula(protein_sequence$`Sequence length`,protein_sequence$`Gene Systematic name`,protein_gene_without3D$geneName)


write.table(protein_gene_without3D, "result/protein_gene_without3D.txt", sep = "\t", row.names = TRUE)

