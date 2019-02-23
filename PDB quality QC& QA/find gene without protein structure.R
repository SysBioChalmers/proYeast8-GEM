# this code was mainly used to find sequences without PDB files based on PDB files
# information from uniprot, swiss and pdb files
# sometimes it can be found maybe the information in the uniprot is not updated on time

library(readxl)
library(readr)
library(tidyverse)
library(stringr)
library(hongR)

gene_all <- read_excel("data/gene_list_yeastGEM.xlsx",  sheet = "Sheet1")
gene_all$geneNames <- str_trim(gene_all$geneNames, side = "both")

# pdb data from gangLI
yeast_3D_information <- read_excel("data/yeast_3D_swiss_July.xls")
yeast_3D_information$locus <- str_trim(yeast_3D_information$locus, side = "both")

# pdb data from uniprot
uniprot_yeast_protein_3D <- read_excel("data/yeast_3D_uniprot.xlsx", sheet = "Sheet0")
uniprot_yeast_protein_3D$`Gene names  (ordered locus )` <- str_trim(uniprot_yeast_protein_3D$`Gene names  (ordered locus )`, side = "both")
uniprot_yeast_protein_3D <- filter(uniprot_yeast_protein_3D, is.na(pdb)==FALSE)
uniprot_protein_3D <- splitAndCombine(uniprot_yeast_protein_3D$pdb,uniprot_yeast_protein_3D$`Gene names  (ordered locus )`, sep0 = ";")
colnames(uniprot_protein_3D) <- c("pdb_id", "gene")
uniprot_protein_3D <- filter(uniprot_protein_3D, pdb_id != "")

# pdb file statistics analysis
gene_all$pdb_swiss <- getMultipleReactionFormula(yeast_3D_information$pdb_id,yeast_3D_information$locus,gene_all$geneNames)
gene_all$pdb_ex <- getMultipleReactionFormula(uniprot_protein_3D$pdb_id,uniprot_protein_3D$gene,gene_all$geneNames)


# find the gene without 3D in both swiss and PDB database
gene_with_3D <- filter(gene_all, !is.na(pdb_swiss) | !is.na(pdb_ex))
gene_without_3D <- filter(gene_all, is.na(pdb_swiss) & is.na(pdb_ex))


# obtain the protein sequence for these without 3D
protein_sequence <- read_tsv(file = 'data/yeast_protein_sequence_SGD.tsv')
gene_without_3D$sequence <- getSingleReactionFormula(protein_sequence$`Sequence residue`,protein_sequence$`Gene Systematic name`,gene_without_3D$geneNames)
gene_without_3D$length <- getSingleReactionFormula(protein_sequence$`Sequence length`,protein_sequence$`Gene Systematic name`,gene_without_3D$geneNames)

write.table(gene_without_3D, "result/protein_gene_without3D.txt", sep = "\t", row.names = TRUE)

# summary
# it can be found that 181 gene without PDB from swiss model and pdb database, thus we will do the model simulation for these
# structure manually




