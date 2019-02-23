# this code is used to classify the pdb file into two types
# the first type contains xx.pdb
# the second type not contains xx.pdb
# 12th, November, 2018
# hongzhong lu

source('Preprocess data from swiss database.R')
source('main_function_for_quality_analysis.R')

# function used in this part
# calculate the number between a range (number1, number2)
# input the gene information in yeastGEM
gene_all <- read_excel("data/gene_list_yeastGEM.xlsx",  sheet = "Sheet1")
gene_all$geneNames <- str_trim(gene_all$geneNames, side = "both")


index_experiment <- which(model_EXP$locus %in% gene_all$geneNames ==TRUE)
pdb_EX <- model_EXP[index_experiment,]
gene_EX <- unique(pdb_EX$locus)

# the followed code is to divide the pdb files into two files:
# file 1 with right pdb format
# file 2 with wrong pdb format
pdb_EX_template <- unique(pdb_EX$template)
pdb_EX_file <- list.files("PDB experimental pdb files")
pdb_EX_template0 <- paste(pdb_EX_template, ".pdb", sep = "")
pdb_EX_template01 <- intersect(pdb_EX_file,pdb_EX_template0)

#establish a file to store the pdb_EX with right format
dir.create("pdb_EX_right_format")
for (i in seq_along(pdb_EX_template01)) {
  s <- pdb_EX_template01[i]
  file0 <- paste("PDB experimental pdb files/", s, sep = "")
  file.copy(file0, "pdb_EX_right_format")
}

#establish a file to store the pdb_EX with wrong format
pdb_EX_template1 <- setdiff(pdb_EX_template0, pdb_EX_file) 
pdb_EX_template11 <- str_replace_all(pdb_EX_template1,".pdb", "-pdb-bundle.tar")

intersect(pdb_EX_template11, pdb_EX_file)

dir.create("pdb_EX_other_format")
for (i in seq_along(pdb_EX_template11)) {
  s <- pdb_EX_template11[i]
  file0 <- paste("PDB experimental pdb files/", s, sep = "")
  file.copy(file0, "pdb_EX_other_format")
}





