library(tidyverse)
library(stringr)
library(readxl)

#which gene have mutation fromm 1011 project
geneName <- list.files('1011_project')
geneName <- str_replace_all(geneName,'\\.fasta','')
geneMetabolic <- intersect(geneName, gene_feature_GEM$locus_tag)

#YAL012W
ss = 'YAL012W'
ss %in% gene_feature_GEM$locus_tag
gene_snp <- getGeneCoordinate(gene_name = ss, genesum = gene_feature_GEM)
length(gene_snp$gene)


#count the mutant number in each protein position
#in the first step. read fasta file for each gene and change it into dataframe or list format for next mapping analysis
#first run a single gene-YAL012W.fasta from 1011 project
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
library("Biostrings")
fastaFile <- readDNAStringSet("data/YAL012W.fasta")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence, stringsAsFactors = FALSE)
df['geneName'] <- 'YAL012W'

df_list <- list()
list_name <- df$geneName[1]
for (i in seq(length(df$sequence))){
  df_list[i] <- str_split(df$sequence[i],"")
}


findPPosition0 <- function(alted_seq, geneName){
  #this function is used to find the postion of mutated amino acids based on genomics mutation
  #alted_seq <- gene1_s0
  #geneName <- 'YAL012W'
  gene_snp <- getGeneCoordinate(gene_name = geneName, genesum = gene_feature_GEM)
  gene_snp[['gene']] <- alted_seq
  
  #translation
  library(seqinr)
  realcds <- str_to_lower(paste(gene_snp[['gene']],collapse = ""))
  toycds <- s2c(realcds)
  gene_snp[['protein_mutated']] <- translate(seq = toycds)
  
  #find the relative postion of mutated amino acids
  aa_position <- which(gene_snp[['protein']] != gene_snp[['protein_mutated']] )
  
  #calculate the mutation number in the mutated postion (for specific strain -x)
  gene_snp[['mutation_position']] <- rep(0,length(gene_snp[['protein']])) #initialize the start value for each positions
  gene_snp[['mutation_position']][aa_position] <- 1
  result <- unlist(gene_snp[['mutation_position']])
  return(result)
}


countMutationProtein0 <- function (gene_name, mutated_gene_seq){

  #mutated_gene_seq <- df_list
  #geneName = 'YAL012W'
  
  gene_snp <- getGeneCoordinate(gene_name = geneName, genesum = gene_feature_GEM)
  tt <- rep(0,length(gene_snp[['protein']]))
  for (i in seq(length(mutated_gene_seq))){
    
    tt <- tt + findPPosition0(df_list[[i]],geneName)
  }
  
  return(tt)
}


ss = 'YAL012W'
gene_snp <- getGeneCoordinate(gene_name = ss, genesum = gene_feature_GEM)
gene_snp[['pro_mutation_count']] <- countMutationProtein0(gene_name = ss, mutated_gene_seq = df_list)
pos_mutation <- which(gene_snp[['pro_mutation_count']] != 0)







































