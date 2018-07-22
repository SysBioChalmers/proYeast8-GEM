library(tidyverse)
library(stringr)
library(readxl)
library("Biostrings")
library(filesstrings) #move the files
#source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")


#this function is used to filter the whole sample
filterMetabolicGene <- function() {
  geneName0 <- list.files("1011_project")
  geneMetabolic <- paste(str_trim(gene_feature_GEM$locus_tag, side = "both"), ".fasta", sep = "")
  geneMetabolic0 <- intersect(geneName0, geneMetabolic)
  dir.create("target_gene")
  
  for (i in seq_along(geneMetabolic0)) {
    geneX <- geneMetabolic0[i]
    file0 <- paste("1011_project/", geneX, sep = "")
    file.copy(file0, "target_gene")
  }
}

#this function is used to just preprocess the fasta file without filteration
processFasta <- function(gene_test) {
  # read each fasta file and change it into a dataframe
  # exampel:gene_test <- "YAL012W.fasta"
  gene_name_test <- str_replace_all(gene_test, ".fasta", "")
  fastaFile <- readDNAStringSet(paste("target_gene/", gene_test, sep = ""))
  # obtain the strain name information and sequence information
  seq_name <- names(fastaFile)
  sequence <- paste(fastaFile)
  
  # establish a dataframe contains the strain name and sequnece information
  df <- data.frame(seq_name, sequence, stringsAsFactors = FALSE)
  df_list <- list()
  for (i in seq(length(df$sequence))) {
    df_list[i] <- str_split(df$sequence[i], "")
  }
  
  return(df_list)
}

#this function is used to choose the gene based on strain phenotype information
filterMutationStrainType <- function(gene_test, strain_select) {
  # read each fasta file and change it into a dataframe
  # then the filter can be used for each dataframe to obtain the strains we need
  # gene name
  # exampel:gene_test <- "YAL012W.fasta"
  gene_name_test <- str_replace_all(gene_test, ".fasta", "")
  fastaFile <- readDNAStringSet(paste("target_gene/", gene_test, sep = ""))
  # obtain the strain name information and sequence information
  seq_name <- names(fastaFile)
  sequence <- paste(fastaFile)
  
  # establish a dataframe contains the strain name and sequnece information
  df <- data.frame(seq_name, sequence, stringsAsFactors = FALSE)
  strain_select['index_strain'] <- paste(strain_select$Standardized_name, "_", gene_name_test, "_", sep = "")
  
  
  for (j in seq_along(df$seq_name)) {
    exist_sign <- vector()
    for (i in seq_along(strain_select$index_strain)) {
      exist_sign[i] <- str_detect(df$seq_name[j], strain_select$index_strain[i])
    }
    
    if (any(exist_sign) == TRUE) {
      df$choosed[j] <- "YES"
    } else {
      df$choosed[j] <- "NO"
    }
  }
  
  df_refine <- filter(df, choosed == "YES")
  
  
  df_list <- list()
  for (i in seq(length(df_refine$sequence))) {
    df_list[i] <- str_split(df_refine$sequence[i], "")
  }
  
  dir.create("target_gene_processed")
  filename0 <- paste("target_gene_processed/",gene_name_test, ".RData", sep = "")
  save(df_list, file = filename0)
  return(df_list)
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


countMutationProtein0 <- function (geneName, mutated_gene_seq){

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
  m <- length(mutated_pos)
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
calculateMutNear <- function(col0, pos_mutation_3D) {
  

}























