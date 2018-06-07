library(tidyverse)
library(stringr)


#get the gene name
#try to calculate the mutation on the amino acids based on the coordination on the chromosome
mutated_test <- read.csv('mutated_proteins.txt', header = TRUE, sep = "\t", stringsAsFactors = FALSE)
mutated_test$Chromosome <- str_trim(mutated_test$Chromosome, side = "both")
mutated_test$Position <- as.numeric(mutated_test$Position)
#function get the gene name based on the mutation position
getGeneName <- function(chr,mutated_positions,gene_annotation = gene_feature0){
  ss <- filter(gene_feature0, 
               chromosome == chr & 
                 start <= mutated_positions & 
                 end >= mutated_positions)
  if(length(ss$locus_tag)){
    ss0 <- ss$locus_tag
  } else{
    ss0 <- "INTERGENIC"
  }
  return(ss0)
}
for (i in seq(length(mutated_test$Chromosome))){
  mutated_test$Gene2[i] <- getGeneName(mutated_test$Chromosome[i],mutated_test$Position[i])
}

mutated_test0 <- filter(mutated_test, Gene2 != "INTERGENIC") ##filter the mutated test





#choose the metabolic gene
index_m <- which(mutated_proteins0$Gene2 %in% gene_feature_GEM$locus_tag ==TRUE)
mutated_gene <- mutated_proteins0[index_m,]
mutated_gene$Ref <- str_trim(mutated_gene$Ref, side = "both")
mutated_gene$Alt <- str_trim(mutated_gene$Alt, side = "both")




#first run the program for each gene from different conditions or strains
#get the mutated gene based on the mutation postion 
gene_annotation <- filter(gene_feature_GEM, locus_tag=="YDL205C")

#pre-process the gene annotation data before mutation mapping
gene_feature1 <- gene_annotation
gene_snp <- list()
gene_snp[['gene']] <- unlist(strsplit(gene_feature1$cds_seq[1], split = ""))
# cds_location
gene_snp[['gene_coordinate']] <- seq(93745,92762,-1) #should using a function
gene_snp[['protein']] <- unlist(strsplit(gene_feature1$aa_seq[1], split = ""))
gene_snp[['protein_coordinate']] <- seq(as.numeric(gene_feature1$aa_length[1]))


#if mutation_position existed, get the mutated gene
#input the mutated information of gene from different conditons or strains
mutated_gene0 <- filter(mutated_gene, Gene2=="YDL205C")
     #if the cds from the minus strand

changeATCG <- function (ss){
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
 

for (i in seq(length(mutated_gene0$Chromosome))){
  mutated_gene0$Ref[i] <- changeATCG(mutated_gene0$Ref[i])
  mutated_gene0$Alt[i] <- changeATCG(mutated_gene0$Alt[i])
  }

mutation_position <- which(gene_snp[['gene_coordinate']]==93192)
gene_snp[['gene']][mutation_position] <- "A"








#translation
library(seqinr)
realcds <- str_to_lower(paste(gene_snp[['gene']],collapse = ""))
toycds <- s2c(realcds)
gene_snp[['protein_mutated']] <- translate(seq = toycds)










#find the relative postion of mutated amino acids
aa_position <- which(gene_snp[['protein']] != gene_snp[['protein_mutated']] )









#calculate the mutation number in the mutated postion (for specific strain -x)
gene_snp[['mutation_num_x']] <- rep(0,length(gene_snp[['protein']])) #initialize the start value for each positions
gene_snp[['mutation_num_x']][aa_position] <- 1





#save and load list data
save(gene_snp, file="YFL011W.RData")
load("YFL011W.RData")
