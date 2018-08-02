library(tidyverse)
library(stringr)
library(readxl)
library(Biostrings)
library(filesstrings) #move the files
library(hongR)
#source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
library(centiserve) #this package is used to calculate the closeness centrality
library(igraph)#form the unique clust based on Floyd-Warshall shortest-paths algorithm
library(seqinr)

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


#---------------version 1
#the followed two function  were used to estimate the mutation information based on input DNA fasta file
findPPosition0 <- function(alted_seq, geneName){
  #this function is used to find the postion of mutated amino acids based on genomics mutation
  #alted_seq <- df_list[[1]]
  #geneName <- 'YAL012W'
  gene_snp <- getGeneCoordinate(gene_name = geneName, genesum = gene_feature_GEM)
  gene_snp[['gene']] <- alted_seq
  
  #translation
  #using package seqinr
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
#function to obtain the mutation of amino acids residue in the 3d structure based on sequence blast anlysis
#but it can be very dangeous using this method, as the insertion or deletion could lead to many mutation in a seq
countMutationProtein0 <- function (geneName, mutated_gene_seq){

  #mutated_gene_seq <- df_list
  #geneName = 'YAL012W'
  df_list <- mutated_gene_seq
  gene_snp <- getGeneCoordinate(gene_name = geneName, genesum = gene_feature_GEM)
  tt <- rep(0,length(gene_snp[['protein']]))
  for (i in seq(length(mutated_gene_seq))){
    if(length(gene_snp[['gene']]) != length(df_list[[i]])) {
      tt <- tt + rep(0,length(gene_snp[['protein']]))
       }
    ##to avoide the insertion or deletion in the seq
    else{
      tt <- tt + findPPosition0(df_list[[i]],geneName)
      }
  }
  
  return(tt)
}

#----------------vesion 2 
#the followed two function  were used to estimate the mutation information based on single SNP information
#using function to obtain the each gene's mutation information based on the processed mutation data
#These function is used to obtain the count of mutation on portein level
countMutationProtein <- function (gene_name, mutation_annotation=mutated_gene1){
  #this function could produce the all the results about mutated amino acids information
  #gene_name <- "YDL205C"
  mutated_gene0 <- filter(mutation_annotation, Gene2==gene_name)
  tt <- rep(0,length(gene_snp[['protein']]))
  for (i in seq(length(mutated_gene0$Gene2))){
    tt <- tt + findPPosition(mutated_gene0$Pos[i],mutated_gene0$Alt[i],gene_name)
  }
  
  return(tt)
}
findPPosition <- function(mutatedPosition, alted, geneName){
  #this function is used to find the postion of mutated amino acids based on genomics mutation
  #mutatedPosition = 93192
  #alted ='A'
  mutation_position <- which(gene_snp[['gene_coordinate']]==mutatedPosition)
  gene_snp <- getGeneCoordinate(gene_name = geneName, genesum = gene_feature_GEM)
  gene_snp[['gene']][mutation_position] <- alted
  
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


#if the cds from the minus strand, then the functin changeATCG should be used to firstly 
#get the mutation information on the minus strand based on that from the positive strand
changeATCG <- function (ss){
  # this function was used to get the mutation information from the minus strand based on the mutation information
  # on the positive strand
  if (ss =="A"){
    ss <- "T"
  } else if(ss=="C"){
    ss <- "G"
  } else if(ss=="T"){
    ss <-"A"
  } else if(ss=="G"){
    ss <-"C"
  }
  return(ss)
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


#part 2 function related to hot spot analysis
#this function return the mutated informaiton based on genomics fasta information
#this function is used to establish the relation between the mutated amino acid and related position
PositionResidue <- function(alted_seq, geneName) {
  #alted_seq <- df_list[[6]]
  #geneName <- "YAL012W"
  gene_snp <- getGeneCoordinate(gene_name = geneName, genesum = gene_feature_GEM)
  gene_snp[["gene"]] <- alted_seq
  
  # translation
  library(seqinr)
  realcds <- str_to_lower(paste(gene_snp[["gene"]], collapse = ""))
  toycds <- s2c(realcds)
  gene_snp[["protein_mutated"]] <- translate(seq = toycds)
  
  # find the relative postion of mutated amino acids
  aa_position <- which(gene_snp[["protein"]] != gene_snp[["protein_mutated"]])
  aa_type <- gene_snp[["protein_mutated"]][aa_position]
  
  # built the relation between aa_position and aa_type
  # aa_postion and aa_type should contain one element
  mutatedAA <- paste(aa_type, aa_position, sep = "@@") # this estabolish the relation between the postion and mutated amino acids
  return(mutatedAA)
}


#this function return the mutated informaiton based on simple SNP information: mutated position and changed amino acids
PositionResidue0 <- function(mutatedPosition, alted, geneName) {
  #mutatedPosition = 130975
  #alted ='A'
  #geneName = "YAL012W"
  gene_snp <- getGeneCoordinate(gene_name = geneName, genesum = gene_feature_GEM)
  mutation_position <- which(gene_snp[['gene_coordinate']]==mutatedPosition)
  
  gene_snp[['gene']][mutation_position] <- alted
  
  # translation
  #library(seqinr)
  realcds <- str_to_lower(paste(gene_snp[["gene"]], collapse = ""))
  toycds <- s2c(realcds)
  gene_snp[["protein_mutated"]] <- translate(seq = toycds)
  
  # find the relative postion of mutated amino acids
  aa_position <- which(gene_snp[["protein"]] != gene_snp[["protein_mutated"]])
  aa_type <- gene_snp[["protein_mutated"]][aa_position]
  
  # built the relation between aa_position and aa_type
  # aa_postion and aa_type should contain one element
  mutatedAA <- paste(aa_type, aa_position, sep = "@@") # this estabolish the relation between the postion and mutated amino acids
  return(mutatedAA)
}


#this function is used to put the residue from the same postion together
ResidueSum <- function(pos_residue) {
  pos_residue <- unlist(pos_residue)
  pos_residue_df <- data.frame(pos_aa = pos_residue, pos_aa1 = pos_residue, stringsAsFactors = FALSE)
  pos_residue_df <- pos_residue_df %>% separate(., pos_aa1, into = c("residue", "position"), sep = "@@")
  
  pos_residue_df0 <- data.frame(pos = unique(pos_residue_df$position), stringsAsFactors = FALSE)
  pos_residue_df0$residue <- getMultipleReactionFormula(pos_residue_df$pos_aa, pos_residue_df$position, pos_residue_df0$pos)
  pos_residue_df0$pos <- as.numeric(pos_residue_df0$pos)
  
  return(pos_residue_df0)
}


#this function is from hongR
splitAndCombine <- function(gene, rxn,sep0) { ##one rxn has several genes, this function was used to splite the genes
  
  gene <- str_split(gene, sep0)
  tt<- length(gene)
  gene0 <- list()
  for (i in 1:tt){
    gene0[[i]] <- paste(rxn[i], gene[[i]], sep = "@@@")
    
  }
  
  gene1 <- unique(unlist(gene0)) # the duplicated element is not deleted
  gene2 <- str_split(gene1, "@@@" )
  rxnGene <- data.frame(v1=character(length(gene2)),stringsAsFactors = FALSE)
  tt1 <- length(gene2)
  for (j in 1:tt1){
    rxnGene$v1[j] <- gene2[[j]][2]
    rxnGene$v2[j] <- gene2[[j]][1]
  }
  
  return(rxnGene)
}

#this function is used to get the vertices
getHotVertice <- function(aa_3d, residue0, aa_pro, distance0 = ResidueDistance_1n8p) {
  #input
  #aa_3d  a vector for the coordinate of PDB structure
  #residue  a vector for the muated residue information of the stucture
  #aa_pro a vector for the original coordinate of protein aa sequence
  #ditance  a matrix for the distance of the paired residue of pdb structure
  
  #output
  #dataframe contains the inforation of the mutated residues
  
  # establish the structure coordinate and all the residues
  pos_residue_3d <- data.frame(pos = aa_3d, residue = residue0, stringsAsFactors = FALSE)
  pos_residue_3d <- splitAndCombine(pos_residue_3d$residue, pos_residue_3d$pos, sep0 = "\\;")
  colnames(pos_residue_3d) <- c("residue", "pos_3d")
  pos_residue_3d <- pos_residue_3d[which(pos_residue_3d$residue != "NA"), ]
  # replace the protein aa coordinate into the structure coordinate
  pos_residue_3d <- pos_residue_3d %>% separate(., residue, into = c("residue", "pos_pro"))
  pos_residue_3d$residue <- paste(pos_residue_3d$residue, pos_residue_3d$pos_3d, sep = "@@")
  all_pair <- combn(pos_residue_3d$residue, 2)
  all_pair0 <- as.data.frame(t(all_pair))
  
  # choose the cluste based on p value, Distance between two pair and sperated residues(>20)
  all_pair1 <- all_pair0 %>%
    separate(., V1, into = c("residue", "c1"), sep = "@@") %>%
    separate(., V2, into = c("residue2", "c2"), sep = "@@")
  all_pair1[1:10, ]
  all_pair1$c1 <- as.numeric(all_pair1$c1)
  all_pair1$c2 <- as.numeric(all_pair1$c2)
  
  all_pair_distance <- vector()
  for (i in seq_along(all_pair1$residue)) {
    all_pair_distance[i] <- distance0[all_pair1$c1[i], all_pair1$c2[i]]
  }
  
  # calculate the distance between all pair (how many amino acids separate the pair)
  pos_initial <- 1:length(aa_pro)
  all_distance <- vector()
  all_pair_ini <- combn(pos_initial, 2)
  for (i in 1:ncol(all_pair_ini)) {
    s1 <- all_pair_ini[, i]
    d0 <- distance0[s1[1], s1[2]]
    all_distance[i] <- d0
  }
  
  # first filter based on aa_distance and space distance
  all_pair0$aa_distance <- abs(all_pair1$c1 - all_pair1$c2)
  all_pair0$distance <- all_pair_distance
  all_pair2 <- all_pair0[which((all_pair0$aa_distance > 20 & all_pair0$distance <= 10) | all_pair0$aa_distance == 0), ]
  
  # calculate the P value for each pair of amino acids based on the distance
  for (i in seq_along(all_pair2$distance)) {
    distance0 <- all_pair2$distance[i]
    all_pair2$pvalue_pair[i] <- length(which(all_distance <= distance0)) / length(all_distance)
  }
  
  all_pair3 <- all_pair2[which(all_pair2$pvalue_pair < 0.05), ]
  return(all_pair3)
}


#this function is used to obtain the cluster analysis results
clusterAnalysis <- function(residueInf) {
  #input
  #residueInf  dataframe contains the detailed information about the residue
  
  #output
  # a dataframe contains the cluster information
  
  # obtain the clusters
  links <- data_frame(from = residueInf$V1, to = residueInf$V2, weight = residueInf$distance)
  g <- graph_from_data_frame(d = links, directed = FALSE)
  plot(g)
  # library(networkD3)
  # simpleNetwork(links, fontSize = 18, opacity = 1, zoom = TRUE,fontFamily = "sans-serif")
  # split the graph into subgraph and get the unique cluster
  # calculate the closeness centrality for each clust
  dg <- decompose.graph(g)
  detail_mutant_position0 <- list()
  position_combine <- vector()
  for (i in seq_along(dg)) {
    clust2 <- dg[[i]][1]
    detail_mutant_position0[[i]] <- names(clust2)
    position_combine[i] <- paste0(names(clust2), collapse = ";")
  }
  
  
  closeness0 <- list()
  cluster_closeness <- vector()
  dg <- decompose.graph(g)
  for (i in seq_along(dg)) {
    ##weights represent the distance between the node, has been contained in dg, E(dg).weight could be
    ##used to check the weiht information in each subgraph
    closeness0[[i]] <- closeness.residual(dg[[i]]) 
    cluster_closeness[i] <- sum(closeness.residual(dg[[i]]))
  }
  spotSummary <- data_frame(cluster = position_combine, closeness = cluster_closeness)
  return(spotSummary)
}


#this function is used to calculate p value for each clust
getHotPvalue <- function(cluster0, sample_standard, distance, seq) {
  #cluster0 <- important_hot$cluster
  #input
  # cluster0    a vector contains the detailed mutation postion information for each cluster
  # sample_standard   a vector contains the standard number of mutated residue in each position
  # ditance a matrix contains the ditance for each paired residue
  # seq   a vector contains the coordinate information of each amino acid
  #output
  # a vector contains the calculated pvalue for each cluster
  
  
  # obtain the detaile postion for each cluster
  cluster1 <- str_split(cluster0, ";")
  str_replace_all("X@@75", "[:alpha:]@@", "")
  for (i in seq_along(cluster1)) {
    cluster1[[i]] <- str_replace_all(cluster1[[i]], "[:alpha:]@@", "")
  }
  
  for (i in seq_along(cluster1)) {
    cluster1[[i]] <- unique(as.numeric(cluster1[[i]]))
  }
  
  
  # calculate the p_value
  pvalue <- vector()
  for (i in seq_along(cluster1)) {
    if (length(cluster1[[i]]) >= 2) {
      pos_mutation_t0 <- cluster1[[i]]
      wap_original <- getTotalWAP(pos_mutation_t0, sample_standard, distance)
      wap_sample0 <- getSampleWAP(pos_mutation_t0, sample_standard, distance, seq, n = 10000)
      pvalue[i] <- getPvalue(wap_original, wap_sample0)
    } else {
      pvalue[i] <- 1
    }
  }
  
  return(pvalue)
}











