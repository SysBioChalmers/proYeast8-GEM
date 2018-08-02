# this script is showing how to analyze the data for the SNP data in the excel format
# it is also base for hotspot and clumps analysis method
# the SNP format is as followed:
# Chr       pos      Ref       Alt
# chr1      19971    C         T    
# chr1      19972    T         C


#library
library(tidyverse)
library(stringr)
library(readxl)

#common function
getSingleReactionFormula <- function(description, reaction_ko, ko) {###description can be any charater of metabolite
  index <- vector()
  result <- vector()
  tt <- vector()
  for (i in 1:length(ko)){
    if(length(match(ko[i],reaction_ko))){
      index <- match(ko[i],reaction_ko)
      tt <- description[index]
      result[i] <- paste0(tt, collapse = ";")
    } else{
      
      result[i] <- NA
    }
  }
  return(result)
}

#get the gene name
#try to calculate the mutation on the amino acids based on the coordination on the chromosome
mutated_test <- read_excel("data/snp_adaption_to_high_ethanol.XLS")
mutated_test$Chr <- str_trim(mutated_test$Chr, side = "both")
mutated_test$Pos <- as.numeric(mutated_test$Pos)

#function get the gene name based on the mutation position
getGeneName <- function(chr,mutated_positions,gene_annotation = gene_feature0){
  #input:
  #1. chr: chromsome
  #2. mutated_positiion
  #3. gene_featured0: contains the gene sequence information from chromsome of sec-s288c ,like the start and end
  #output:
  # the gene name contained this mutation
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


for (i in seq(length(mutated_test$Chr))){
  mutated_test$Gene2[i] <- getGeneName(mutated_test$Chr[i],mutated_test$Pos[i])
}

mutated_test0 <- filter(mutated_test, Gene2 != "INTERGENIC") ##filter the mutated test

#choose the metabolic gene
#if the gene is type of "complement", then the complement_sign is "TRUE"
#else the complement_sign is "FALSE"
gene_feature_GEM$complement_sign <- str_detect(gene_feature_GEM$cds_location,"complement")
index_m <- which(mutated_test0$Gene2 %in% gene_feature_GEM$locus_tag ==TRUE)
mutated_gene <- mutated_test0[index_m,]


mutated_gene$Ref <- str_trim(mutated_gene$Ref, side = "both")
mutated_gene$Alt <- str_trim(mutated_gene$Alt, side = "both")



#The followed part is based on one gene under different condition
#mutated gene information preprocess
#if mutation_position existed, get the mutated gene
#input the mutated information of gene from different conditons or strains

mutated_gene$complement_sign <- getSingleReactionFormula(gene_feature_GEM$complement_sign,gene_feature_GEM$locus_tag,mutated_gene$Gene2)
mutated_gene1 <- mutated_gene

for (i in seq(length(mutated_gene1$Chr))){
  if(mutated_gene1$complement_sign[i]){
  mutated_gene1$Ref[i] <- changeATCG(mutated_gene1$Ref[i])
  mutated_gene1$Alt[i] <- changeATCG(mutated_gene1$Alt[i])
  
  } else{
    mutated_gene1$Ref[i] <- mutated_gene1$Ref[i]
    mutated_gene1$Alt[i] <- mutated_gene1$Alt[i]
  }
}






#first run the program for each gene from different conditions or strains
#pre-process the gene annotation data before mutation mapping
#update the mutation information in the protein level
gene_list  <- unique(mutated_gene1$Gene2)
tt <- vector()
for (i in seq_along(gene_list)){
ss = gene_list[i]
gene_snp <- getGeneCoordinate(gene_name = ss, genesum = gene_feature_GEM)
gene_snp[['pro_mutation_count']] <- countMutationProtein(gene_name = ss, mutation_annotation=mutated_gene1)
tt[i] <- sum(gene_snp[['pro_mutation_count']])
print(gene_snp)
}


#results analysis
num_gene_with_nsSNP <- tt[tt > 0]
num_nsSNP <- sum(num_gene_with_nsSNP)
protein_mutation <- data.frame(orf=gene_list,nsSNP=tt)
#import the annotation of these protein
gene_annotation <- read.delim2("data/all_gene_yeast with annotation from different database.txt")
protein_mutation$annotation_sgd <- getSingleReactionFormula(gene_annotation$annotation_SGD,gene_annotation$gene,protein_mutation$orf)
protein_mutation$annotation_uni <- getSingleReactionFormula(gene_annotation$annotation_uniprot,gene_annotation$gene,protein_mutation$orf)
protein_mutation0 <- protein_mutation[protein_mutation$nsSNP >= 4,]
write.table(protein_mutation0,"result/protein_mutation0.txt", row.names = FALSE, sep = "\t" )



#analyze the mutation information with the structure
ss = 'YBR115C'
seq_from_3D <- 2:919 #this is the coordinated of original protein sequence and should changed into 3D structure coordinates
dirForDistanceMatrix <- paste("data/ResidueDistance_",ss,".xlsx", sep = "")


gene_snp <- getGeneCoordinate(gene_name = ss, genesum = gene_feature_GEM)
gene_snp[['pro_mutation_count']] <- countMutationProtein(gene_name = ss, mutation_annotation=mutated_gene1)
pos_mutation <- which(gene_snp[['pro_mutation_count']] != 0)


#input the distance of all the pired residues
ResidueDistance <- read_excel(dirForDistanceMatrix,col_names = FALSE) #in the followed calculation, the matrix dosen't have the col and row names
ResidueDistance <- as.matrix(ResidueDistance)


#obtain the mutation information for the structure
residueIn3D <- gene_snp[['protein']][seq_from_3D]
pos_mutation_3D <- gene_snp[['pro_mutation_count']][seq_from_3D]


#mutation position on structure and #mutation number on structure
pos_mutation_c <- which(pos_mutation_3D != 0)
seq0 <- 1:length(pos_mutation_3D) #seq0 is the coordinate of PDB structure
pos_count_num <- pos_mutation_3D[pos_mutation_c]


#calculate p_values using UPMS method
sample_standard1 <- sampleStand(pos_mutation_3D)
wap_original <- getTotalWAP(pos_mutation_c,sample_standard1,ResidueDistance)
wap_sample0 <- getSampleWAP(pos_mutation_c,sample_standard1,ResidueDistance, seq=seq0,n=10000)
plotNullDistribution(wap_sample0)
Strain_3D <- getPvalue(wap_original,wap_sample0)







