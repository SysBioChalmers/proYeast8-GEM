library(readxl)
library(readr)
library(tidyverse)
library(stringr)
library(hongR)

# First download the latest meta information from swiss model database
# The PDB files can also be found in the link
# https://swissmodel.expasy.org/repository/species/559292


model <- read.csv2("data/metaData_swiss_April.txt", sep = "\t", stringsAsFactors = FALSE)
geneIDmapping <- read_excel("data/uniprotGeneID_mapping.xlsx", sheet = "Sheet0")

model$id_mapping <- paste(model$UniProtKB_ac, model$template, sep = "@")
model$locus <- getSingleReactionFormula(geneIDmapping$GeneName,geneIDmapping$Entry,model$UniProtKB_ac)



homology_model <- read_excel("data/Homology_model_info_gangLi_April_2018.xls",
                      sheet = "model_info")
homology_model$id_mapping <- paste(homology_model$UniProtKB_ac, homology_model$SMTLE,sep = "@")
homology_model$Built_with <- paste(homology_model$ENGIN, homology_model$VERSN,sep = " ")



# obtain the model simulated by swiss
model_EXP <- filter(model,provider =="PDB")
model_swiss <- filter(model,provider =="SWISSMODEL")

# add missing information for model_EXP from old file provided by Gang
yeast_3D_swiss <- read_excel("data/yeast_3D_swiss.xls")
yeast_3D_EXP <- filter(yeast_3D_swiss,struct_is_experimental==TRUE)
yeast_3D_EXP$id_mapping <- paste(yeast_3D_EXP$seq_uniprot,yeast_3D_EXP$pdb_id,sep = "@")
model_EXP$Resolution <- getSingleReactionFormula(yeast_3D_EXP$Resolution,yeast_3D_EXP$id_mapping,model_EXP$id_mapping)
model_EXP$Ligands <- getSingleReactionFormula(yeast_3D_EXP$Ligands,yeast_3D_EXP$id_mapping,model_EXP$id_mapping)
model_EXP$locus <- getSingleReactionFormula(geneIDmapping$GeneName,geneIDmapping$Entry,model_EXP$UniProtKB_ac)
model_EXP$Seq_Identity <- "NA"
model_EXP$Seq_similarity <- "NA"

mutation <- read_excel("data/map_exppdb_prot_id_April.xlsx") #input the mutation information for the Experimental files
mutation0 <- mutation %>% separate(qseqid, into = c("PDBid","chain"), sep = "_" )
mutation0$UniProtKB_ac <- getSingleReactionFormula(geneIDmapping$Entry,geneIDmapping$GeneName,mutation0$sseqid)
mutation0$id_mapping <- paste(mutation0$UniProtKB_ac, mutation0$PDBid,sep = "@")

model_EXP$pident <- getSingleReactionFormula(mutation0$pident,mutation0$id_mapping,model_EXP$id_mapping)
model_EXP$length <- getSingleReactionFormula(mutation0$length,mutation0$id_mapping,model_EXP$id_mapping)
model_EXP$mismatch <- getSingleReactionFormula(mutation0$mismatch,mutation0$id_mapping,model_EXP$id_mapping)
model_EXP$gapopen <- getSingleReactionFormula(mutation0$gapopen,mutation0$id_mapping,model_EXP$id_mapping)
model_EXP$qstart <- getSingleReactionFormula(mutation0$qstart,mutation0$id_mapping,model_EXP$id_mapping)
model_EXP$qend <- getSingleReactionFormula(mutation0$qend,mutation0$id_mapping,model_EXP$id_mapping)
model_EXP$sstart <- getSingleReactionFormula(mutation0$sstart,mutation0$id_mapping,model_EXP$id_mapping)
model_EXP$send <- getSingleReactionFormula(mutation0$send,mutation0$id_mapping,model_EXP$id_mapping)
model_EXP$chain <- getSingleReactionFormula(mutation0$chain,mutation0$id_mapping,model_EXP$id_mapping)





# add missing information for model_swiss from homology_model
yeast_3D_swiss <- read_excel("data/yeast_3D_swiss.xls")
yeast_3D_homo <- filter(yeast_3D_swiss,struct_is_experimental==FALSE)
yeast_3D_homo$id_mapping <- paste(yeast_3D_homo$seq_uniprot,yeast_3D_homo$pdb_id,sep = "@")
model_swiss$Resolution <- getSingleReactionFormula(yeast_3D_homo$Resolution,yeast_3D_homo$id_mapping,model_swiss$id_mapping)
model_swiss$Ligands <- getSingleReactionFormula(yeast_3D_homo$Ligands,yeast_3D_homo$id_mapping,model_swiss$id_mapping)
model_swiss$locus <- getSingleReactionFormula(geneIDmapping$GeneName,geneIDmapping$Entry,model_swiss$UniProtKB_ac)
model_swiss$Seq_Identity <- getSingleReactionFormula(homology_model$SID,homology_model$id_mapping,model_swiss$id_mapping)
model_swiss$Seq_similarity <- getSingleReactionFormula(homology_model$SIM,homology_model$id_mapping,model_swiss$id_mapping)
model_swiss$Oligo_State <- getSingleReactionFormula(homology_model$OSTAT,homology_model$id_mapping,model_swiss$id_mapping)
model_swiss$GEMQE <- getSingleReactionFormula(homology_model$GMQE,homology_model$id_mapping,model_swiss$id_mapping)
model_swiss$Built_with <- getSingleReactionFormula(homology_model$Built_with,homology_model$id_mapping,model_swiss$id_mapping)
model_swiss$Found_by <- getSingleReactionFormula(homology_model$FOUND,homology_model$id_mapping,model_swiss$id_mapping)
model_swiss$Method <- getSingleReactionFormula(homology_model$MTHD,homology_model$id_mapping,model_swiss$id_mapping)



# merge the personal simulation information for proteins which can't be found from present swiss database

     # input the gene information in yeastGEM
gene_all <- read_excel("data/gene_list_yeastGEM.xlsx",  sheet = "Sheet1")
gene_all$geneNames <- str_trim(gene_all$geneNames, side = "both")

    #classify the genes from yeastGEM
    #type1: genes whose pdb files can be found from swiss database directly
    #type1: genes whose pdb files can not be found from swiss database directly
index0 <- which(model$locus %in% gene_all$geneNames ==TRUE)
gene_t1 <- unique(model$locus[index0])
gene_t2 <- setdiff(gene_all$geneNames, gene_t1) #total 207 metabolic gene need manual PDB files simulation


    #sumamrize the PDB files and sequence for gene in type2
    #refine the PDB files for 216 metabolic genes which can't be obtained from swiss database
    #These PDB files are obtained one by one using swiss model web service
model2 <- read_excel("data/proteinStructure_manual check.xlsx", sheet = "Sheet1")
seqence_gene_without3D <- read_excel("data/protein_gene_without3D.xlsx", sheet = "Sheet1")
model2$geneID0 <- str_replace_all(model2$geneID,"_2018-03-25","") #geneID0 For mapping of gene with '-'
seqence_gene_without3D$geneID0 <- str_replace_all(seqence_gene_without3D$geneName,"-","")

model2$wild_sequence <- getSingleReactionFormula(seqence_gene_without3D$sequence,seqence_gene_without3D$geneID0,model2$geneID0)
model2$uniprot_seq_length <- getSingleReactionFormula(seqence_gene_without3D$length,seqence_gene_without3D$geneID0,model2$geneID0)
model2$UniProtKB_ac <- getSingleReactionFormula(geneIDmapping$Entry,geneIDmapping$GeneName,model2$locus)
model2$locus <- getSingleReactionFormula(seqence_gene_without3D$geneName,seqence_gene_without3D$geneID0,model2$geneID0)
model2 <- model2 %>%
  separate(Range, c("from", "to"), " - ")
model2$provider <- "SWISSMODEL"
model2$pdb_sequence <- NA
    # filter based on gene_t2
index1 <- which(model2$locus %in% gene_t2 ==TRUE)
model2 <- model2[index1,]


#merge model_swiss and model2 to form the model_homo
model_swiss0 <- select(model_swiss,
                      UniProtKB_ac, 
                      locus,
                      uniprot_seq_length,
                      provider,
                      from,
                      to,
                      coverage,
                      template,
                      qmean,
                      Seq_similarity,
                      Seq_Identity,
                      Resolution,
                      Ligands,
                      Oligo_State,
                      GEMQE,
                      Built_with,
                      Method
                      )

model20 <- select(model2,
                  UniProtKB_ac,
                  locus,
                  uniprot_seq_length,
                  provider,
                  from,
                  to,
                  coverage,
                  template,
                  qmean,
                  Seq_similarity,
                  Seq_Identity,
                  Resolution,
                  Ligands,
                  Oligo_State,
                  GEMQE,
                  Built_with,
                  Method
)


model_homo <- rbind.data.frame(model_swiss0,model20)

##for the homo model, the structure ID given by swiss model should be provided







