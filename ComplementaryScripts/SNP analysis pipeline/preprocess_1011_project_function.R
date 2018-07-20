library(tidyverse)
library(stringr)
library(readxl)
library("Biostrings")
#count the mutant number in each protein position
#in the first step. read fasta file for each gene and change it into dataframe or list format for next mapping analysis
#first run a single gene-YAL012W.fasta from 1011 project
#source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")

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

#fuction to calculate the standard samples contained the mutation in the specific postion
sampleStand <- function (sample_num){
  ss <- sample_num^3/(2^3+sample_num^3)
  return(ss)
}


#function to calculate the initial WAP value
getTotalWAP <- function (mutated_pos, sample0, distance){
  #input
  #mutated_pos   a vector contained all mutated psotion
  #sample0       a vector contained the strandard num of sample contains the above mutation
  #distance      a matrix contrains the pair distance
  
  #output
  #wap1          a num
  
  all_pair <- combn(mutated_pos, 2)
  wap1 <- 0
  t=6
  for (i in 1:ncol(all_pair)){
    s1 <- all_pair[,i]
    n1 <- sample0[s1[1]]
    n2 <- sample0[s1[2]]
    d0 <- distance[s1[1],s1[2]]
    #calculate wap for each pair
    wap1 <- wap1 + n1*n2*exp(-d0[1]^2/2/t^2)
  }
  
  return(wap1)
  
}

#function to calculate the sample WAP
getSampleWAP <- function(mutated_pos, sample0, distance, seq=seq0, n=10000){
  m <- length(pos_mutation_t)
  fixed_sample <- sample0[mutated_pos]
  
  wap_sample <- vector()
  for (i in 1:n){
    sample_position <- sample(seq, m, replace = FALSE, prob = NULL)
    sample_standard_zero <- rep(0,length(seq))
    sample_standard_zero[sample_position] <- fixed_sample
    pos_mutation_t0 <- which(sample_standard_zero != 0)
    wap_sample[i] <- getTotalWAP(pos_mutation_t0,sample_standard_zero,distance)
    print(sample_position)
  }
  
  return(wap_sample)
}

#function to plot the desity graph
plotNullDistribution <- function(wap_sample) {
  plot(density(wap_sample),
       frame = FALSE, col = "steelblue",
       main = "Density",
       xlab = "coverage",
       ylab = "Density"
  )
  
  
  plot(ecdf(wap_sample),
       main = "Cumulative density",
       xlab = "Coverage",
       ylab = "Cumulative"
  )
}

#function to calculate the right-tailed p value
getPvalue <- function(wap_initial, wap_sampling) {
  
  #input:
  #wap_initial  a num
  #wap_sampling a vector of wap obtained sampling method
  
  #output:
  #right tailed p value
  
  index1 <- which(wap_sampling >= wap_initial)
  tail.prob <- (length(index1) + 1) / length(wap_sampling)
  print(tail.prob)
  return(tail.prob)
}


#this function is used in DEMO_CLUMPS.r
calculateMutNear <- function(col0) {
  
  #col1 <- unlist(c(ResidueDistance_1n8p[, 1]))
  col1 <- unlist(col0)
  index1 <- which(col1 <= 10)
  count0 <- pos_mutation_3D[index1]
  mutant_total <- sum(count0)
  return(mutant_total)
}























