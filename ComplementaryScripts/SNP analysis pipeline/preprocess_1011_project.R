#which gene have mutation fromm 1011 project
geneName0 <- list.files('1011_project')
geneName0 <- str_replace_all(geneName,'\\.fasta','')
geneMetabolic <- intersect(geneName0, gene_feature_GEM$locus_tag)

#YAL012W
ss = 'YAL012W'
ss %in% gene_feature_GEM$locus_tag
gene_snp <- getGeneCoordinate(gene_name = ss, genesum = gene_feature_GEM)


#count the mutant number in each protein position
#in the first step. read fasta file for each gene and change it into dataframe or list format for next mapping analysis
#first run a single gene-YAL012W.fasta from 1011 project

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


ss = 'YAL012W'
gene_snp <- getGeneCoordinate(gene_name = ss, genesum = gene_feature_GEM)
gene_snp[['pro_mutation_count']] <- countMutationProtein0(gene_name = ss, mutated_gene_seq = df_list)

ResidueDistance_1n8p <- read_excel("data/ResidueDistance_1n8p.xlsx",col_names = FALSE)
ResidueDistance_1n8p <- as.matrix(ResidueDistance_1n8p)

# mapping the mutation postion onto the protein structure
# for 'YAL012W', the PDB structure id is 'https://www.rcsb.org/structure/1n8p'
# the amino acid sequence in structure is from 2:394 while  the original sequence is from 1:394


#obtain the mutation information for the structure
amino_acid_3D <- gene_snp[['protein']][2:394]
pos_mutation_3D <- gene_snp[['pro_mutation_count']][2:394]

#mutation position on structure
pos_mutation_t <- which(pos_mutation_3D != 0)
seq0 <- 1:length(pos_mutation_3D)
#mutation number on structure
pos_count_t <- pos_mutation_3D[pos_mutation_t]


#wap calculation for each pair mutated residue
#step1 calculate the standardard sample number
sample_standard <- sampleStand(pos_mutation_3D)

#step2 calculate the wap for each pair of mutated residues based on mutation postion
wap_original <- getTotalWAP(pos_mutation_t,sample_standard,ResidueDistance_1n8p)


#step3 change the postion of mutation while keep the mutation number in each postion
#only change the postion but not change the mutated number???
wap_sample0 <- getSampleWAP(pos_mutation_t,sample_standard,ResidueDistance_1n8p, seq=seq0,n=10000)

#step4 analyze the result
plotNullDistribution(wap_sample0)
getPvalue(wap_original,wap_sample0)







#calculate the mutant number near a mutated residue
# filter the mutantion num smaller than 3
# just consider the structure affect
# pos_mutation_3D[which(pos_mutation_3D < 3)] <- 0
# pos_mutation_3D[which(pos_mutation_3D >=3)] <- 1

residue_forAnalysis <- pos_mutation_3D[pos_mutation_3D >=3]


#import the residue distance obtained from python code

calculateMutNear <- function(col0) {
  
  #col1 <- unlist(c(ResidueDistance_1n8p[, 1]))
  col1 <- unlist(col0)
  index1 <- which(col1 <= 10)
  count0 <- pos_mutation_3D[index1]
  mutant_total <- sum(count0)
  return(mutant_total)
}


total_mutant <- vector()
total_mutant <- lapply(ResidueDistance_1n8p , calculateMutNear)
total_mutant <- unlist(total_mutant)
mutant_1n8p <- data.frame(amino_acid = amino_acid_3D, position = 2:394, mutant_num = total_mutant )
