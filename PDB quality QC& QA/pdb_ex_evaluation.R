library(stringr)
library(readxl)
library(hongR)
library(tidyverse)

# function to compare whether two lists contain the same elements
SameElements <- function(a, b) return(identical(sort(a), sort(b)))
a=c('a',1,'b')
b=c('b',1,'a')
SameElements(a,b)


PDB_ex <- read.table('result/pdb_EX for PDB structure.txt', header = TRUE, stringsAsFactors = FALSE)
PDB_ex_template <- unique(PDB_ex$template)
PDB_ex$template_id <- paste(PDB_ex$template, PDB_ex$chain, sep = '.')
PDB_ex_template_id <- unique(PDB_ex$template_id)
uniprotID <- unique(PDB_ex$UniProtKB_ac)

#compare the pdb_ex with uniprot database annotation
#it can be found the experimental pdb id from swiss is not the same to uniprot
#which may be due to that the pdb id from uniprot is not updated on-time
#so we mainly based on the swiss model database
yeast_3D_uniprot <- read_excel("data/yeast_3D_uniprot.xlsx")
yeast_3D_uniprot0 <- select(yeast_3D_uniprot, Entry, pdb) %>%
  filter(.,!is.na(pdb))

yeast_3D_uniprot1 <- yeast_3D_uniprot0[which(yeast_3D_uniprot0$Entry %in% uniprotID ==TRUE),]
yeast_3D_uniprot1$pdb_ex_swiss <- getMultipleReactionFormula(PDB_ex$template,PDB_ex$UniProtKB_ac,yeast_3D_uniprot1$Entry)
yeast_3D_uniprot1$pdb_ex_swiss <- str_to_upper(yeast_3D_uniprot1$pdb_ex_swiss)

pdb_uniprot <- str_split(yeast_3D_uniprot1$pdb, ";")
for (i in seq_along(pdb_uniprot)){
  pdb_uniprot[[i]] <- pdb_uniprot[[i]][-length(pdb_uniprot[[i]])]
}

pdb_swiss <- str_split(yeast_3D_uniprot1$pdb_ex_swiss, ";")
for(i in seq_along(yeast_3D_uniprot1$Entry)){
yeast_3D_uniprot1$equal_sign[i] <- SameElements(pdb_uniprot[[i]], pdb_swiss[[i]])

}




#this function is used to filter fasta file of genes which belong to metabolic
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
#it should be careful about the residue sequence in the pdb structures
PDB_ex$id0 <- paste(PDB_ex$locus,PDB_ex$template_id,sep = ".")
pdb_seq_blast_sgd <- read_excel("data/pdb_seq_blast_sgd.xlsx")
pdb_seq_blast_sgd$id0 <- paste(pdb_seq_blast_sgd$orf,pdb_seq_blast_sgd$PDBid, sep = '.')

PDB_ex$pident2 <-  getSingleReactionFormula(pdb_seq_blast_sgd$pident,pdb_seq_blast_sgd$id0,PDB_ex$id0)
PDB_ex$mismatch2 <-  getSingleReactionFormula(pdb_seq_blast_sgd$mismatch,pdb_seq_blast_sgd$id0,PDB_ex$id0)

PDB_ex_100_new <- filter(PDB_ex, pident2 ==100)
s1 <- unique(PDB_ex_100_new$locus)


#compared with the previous version
#in the previous version, the residue sequence is directly from pdb database, but in fact
#this sequence maybe different from the the real residue sequence contained in the pdb structures
#from the followed comparison we could find that the number of structure with pidenty =100 is different
#in these two methods calculation
PDB_ex_100_old <- filter(PDB_ex, pident ==100)
s2 <- unique(PDB_ex_100_old$locus)
s3 <- intersect(s1,s2)






