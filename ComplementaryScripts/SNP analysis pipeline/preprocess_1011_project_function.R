library(tidyverse)
library(stringr)
library(readxl)
library(Biostrings)
library(filesstrings) #move the files
#source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
library(centiserve) #this package is used to calculate the closeness centrality
library(igraph)#form the unique clust based on Floyd-Warshall shortest-paths algorithm

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


#example to calculate closeness centrality
#install.packages("CINNA")
#library(CINNA)
#g <- graph(c(1,2,2,3,3,4,4,2))
#plot(g, edge.arrow.size=.4)
#closeness.residual(g)

#this function is used to conduct the hotspot analysis for each PDB structure
hotSpotAnalysis <- function(pos_mutation_t,
                            sample_standard,
                            distance,
                            pos_initial) {
  #input:
  #pos_mutation_t  a vector contains the mutation position of PDB structure
  #sample_standard  a vector contains the standard mutated sample in each postion
  #distance    a matrix contains the dstance between each paired amino acid
  #pos_initial  a vector contains the postion of each residue (coordinate) in the PDB structure
  
  #output
  #a dataframe contains the targeted clustered mutated position, the closess of each cluster and p-value of each cluster
  

  # calculate the residue distance between all pair (how many amino acids separate the pair)
  all_pair <- combn(pos_mutation_t, 2)
  all_pair_distance <- vector()
  all_pair_list <- vector()
  aa_distance <- vector()
  for (i in 1:ncol(all_pair)) {
    s1 <- all_pair[, i]
    d0 <- distance[s1[1], s1[2]]
    all_pair_distance[i] <- d0
    all_pair_list[i] <- paste(s1[1], s1[2], sep = "@")
    aa_distance[i] <- abs(s1[1] - s1[2])
  }
  
  # calculate the distance between all pair (how many amino acids separate the pair)
  all_distance <- vector()
  all_pair_ini <- combn(pos_initial, 2)
  for (i in 1:ncol(all_pair_ini)) {
    s1 <- all_pair_ini[, i]
    d0 <- distance[s1[1], s1[2]]
    all_distance[i] <- d0
  }
  
  
  # calculate the P value for each pair of amino acids based on the distance
  pvalue_pair <- vector()
  for (i in seq_along(all_pair_distance)) {
    distance0 <- all_pair_distance[i]
    pvalue_pair[i] <- length(which(all_distance <= distance0)) / length(all_distance)
  }
  
  # choose the cluste based on p value, Distance between two pair and sperated residues(>20)
  target_pair <- vector()
  index00 <- which(all_pair_distance <= 10 & pvalue_pair <= 0.05 & aa_distance >= 20)
  target_pair <- all_pair_list[index00]
  target_pair_distance <- all_pair_distance[index00]
  
  # change the pair into the link relation
  target_pair0 <- str_split(target_pair, "@")
  links <- data_frame(from = vector(length = length(target_pair0)), to = vector(length = length(target_pair0)), weight = target_pair_distance)
  for (i in seq_along(target_pair0)) {
    links$from[i] <- target_pair0[[i]][1]
    links$to[i] <- target_pair0[[i]][2]
    links$weight[i] <- target_pair_distance[i]
  }
  
  g <- graph_from_data_frame(d = links, directed = FALSE)
  plot(g)
  
  # split the graph into subgraph and get the unique cluster
  # calculate the closeness centrality for each clust
  dg <- decompose.graph(g)
  detail_mutant_position0 <- list()
  position_combine <- vector()
  for (i in seq_along(dg)) {
    clust2 <- dg[[i]][1]
    detail_mutant_position0[[i]] <- as.integer(names(clust2))
    position_combine[i] <- paste0(as.integer(names(clust2)), collapse = ";")
  }
  
  closeness0 <- list()
  cluster_closeness <- vector()
  dg <- decompose.graph(g)
  for (i in seq_along(dg)) {
    closeness0[[i]] <- closeness.residual(dg[[i]])
    cluster_closeness[i] <- sum(closeness.residual(dg[[i]]))
  }
  
  
  # calculate the pvalue for each clust
  pvalue <- vector()
  for (i in seq_along(detail_mutant_position0)) {
    if (length(detail_mutant_position0[[i]]) >= 2) {
      pos_mutation_t0 <- detail_mutant_position0[[i]]
      wap_original <- getTotalWAP(pos_mutation_t0, sample_standard, distance)
      wap_sample0 <- getSampleWAP(pos_mutation_t0, sample_standard, distance, seq = pos_initial, n = 10000)
      pvalue[i] <- getPvalue(wap_original, wap_sample0)
    } else {
      pvalue[i] <- 1
    }
  }
  
  # summarize the result
  spotSummary <- data_frame(cluster = position_combine, closeness = cluster_closeness, p_value = pvalue)
  
  return(spotSummary)
}





















