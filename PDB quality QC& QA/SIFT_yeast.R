# this code is used to check the chainID annotation for each pair of PDBid and proteinid
# hongzhong lu
# 12th, November


# load main package
library(readr)
library(readxl)
library(hongR)

pdb_chain_uniprot <- read_delim("data/pdb_chain_uniprot.csv",  ";", escape_double = FALSE, trim_ws = TRUE)
pdb_chain_uniprot$id_mapping <- paste(pdb_chain_uniprot$SP_PRIMARY,pdb_chain_uniprot$PDB, sep = "@")
uniprotGeneID_mapping <- read_excel("data/uniprotGeneID_mapping.xlsx")
pdb_chain_sce <- pdb_chain_uniprot[pdb_chain_uniprot$SP_PRIMARY %in% uniprotGeneID_mapping$Entry, ]

#compare the PDB id from two sources
#this analysis result showed that the sift result is more comprehensive
#thus the sift database can be a very good resource
#the website of sift database:
#http://www.ebi.ac.uk/pdbe/docs/sifts/
yeast_3D_uniprot <- read_excel("data/yeast_3D_uniprot.xlsx")
yeast_3D_uniprot$pdb_sift <- getMultipleReactionFormula(pdb_chain_sce$PDB,pdb_chain_sce$SP_PRIMARY,yeast_3D_uniprot$Entry)

pdb_sce <- unique(pdb_chain_sce$PDB)
length(pdb_sce)

#produce the pdbID summary to download the pdb sequence
pdbID <- paste0(pdb_sce,collapse = ",")
