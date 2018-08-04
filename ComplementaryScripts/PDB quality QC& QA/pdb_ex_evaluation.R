library(stringr)
library(readxl)
library(hongR)
library(tidyverse)
getwd()
PDB_ex <- read.table('result/pdb_EX for PDB structure.txt', header = TRUE, stringsAsFactors = FALSE)
PDB_ex_template <- unique(PDB_ex$template)
PDB_ex$template_id <- paste(PDB_ex$template, PDB_ex$chain, sep = '.')
PDB_ex_template_id <- unique(PDB_ex$template_id)


#this function is used to filter fasta file of genes which belong to metabolic
list.files('result')
PDB_ex_file <- list.files("PDB experimental pdb files")
PDB_ex_template0 <- paste(PDB_ex_template, ".pdb", sep = "")
PDB_ex_template01 <- intersect(PDB_ex_file,PDB_ex_template0)

#establish a file to store the pdb_ex with right format
dir.create("pdb_ex_right_format")
for (i in seq_along(PDB_ex_template01)) {
  s <- PDB_ex_template01[i]
  file0 <- paste("PDB experimental pdb files/", s, sep = "")
  file.copy(file0, "pdb_ex_right_format")
}


#establish a file to store the pdb_ex with wrong format
PDB_ex_template1 <- setdiff(PDB_ex_template0, PDB_ex_file) 
PDB_ex_template11 <- str_replace_all(PDB_ex_template1,".pdb", "-pdb-bundle.tar")

intersect(PDB_ex_template11, PDB_ex_file)

dir.create("pdb_ex_other_format")
for (i in seq_along(PDB_ex_template11)) {
  s <- PDB_ex_template11[i]
  file0 <- paste("PDB experimental pdb files/", s, sep = "")
  file.copy(file0, "pdb_ex_other_format")
}

#produce file to download the pdb file
pdb_ex_download <- str_replace_all(PDB_ex_template1,".pdb", "")
pdb_ex_download <- paste0(pdb_ex_download, collapse = ",")

#after we download all the pdb structure, then we can re-run the above codes


#Evaluate the structure quality based on the seq extracted from pdb files

PDB_ex$id0 <- paste(PDB_ex$locus,PDB_ex$template_id,sep = ".")
pdb_seq_blast_sgd <- read_excel("data/pdb_seq_blast_sgd.xlsx")
pdb_seq_blast_sgd$id0 <- paste(pdb_seq_blast_sgd$orf,pdb_seq_blast_sgd$PDBid, sep = '.')

PDB_ex$pident2 <-  getSingleReactionFormula(pdb_seq_blast_sgd$pident,pdb_seq_blast_sgd$id0,PDB_ex$id0)
PDB_ex$mismatch2 <-  getSingleReactionFormula(pdb_seq_blast_sgd$mismatch,pdb_seq_blast_sgd$id0,PDB_ex$id0)

PDB_ex_100_new <- filter(PDB_ex, pident2 ==100)
s1 <- unique(PDB_ex_100_new$locus)


PDB_ex_100_old <- filter(PDB_ex, pident ==100)
s2 <- unique(PDB_ex_100_old$locus)

s3 <- intersect(s1,s2)






