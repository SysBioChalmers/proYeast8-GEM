#step 1 fasta file preparation
#calculation of the mutated residue in each position

#choose the metabolic gene file
filterMetabolicGene()

#in this project, a fasta file could contained different strains, so we need a function to refine each fasta file
#so that this file just contained the strain or mutation we are interested
strain_classification <- read_excel("data/strain_classification.xls")
unique(strain_classification$Clades)

#run the analysis pipeline using an example
#choose the strain in the clades of bioethanol
strain_select1 <- filter(strain_classification, Clades =='3. Brazilian bioethanol') %>%
  select(.,Standardized_name)

#read each fasta file and change it into a dataframe
#then the filter can be used for each dataframe to obtain the strains we need
gene_test1 <- "YAL012W.fasta"
df_list_filter <- filterMutationStrainType(gene_test = gene_test1, strain_select = strain_select1)

#if no filter, we just preprocess the fasta file with all strain
df_list <- processFasta(gene_test = gene_test1)

#count the mutant number in each protein position
ss <- str_replace_all(gene_test1, ".fasta","")
gene_snp <- getGeneCoordinate(gene_name = ss, genesum = gene_feature_GEM)
gene_snp[['pro_mutation_count']] <- countMutationProtein0(geneName  = ss, mutated_gene_seq = df_list)



#step 2 input the structure information
# mapping the mutation postion onto the protein structure
# for 'YAL012W', the PDB structure id is 'https://www.rcsb.org/structure/1n8p'
#input the distance of all the pired residues
ResidueDistance_1n8p <- read_excel("data/ResidueDistance_1n8p.xlsx",col_names = FALSE)
ResidueDistance_1n8p <- as.matrix(ResidueDistance_1n8p)

#the amino acid sequence in structure is from 2:394 while  the original sequence is from 1:394
#obtain the mutation information for the structure
amino_acid_3D <- gene_snp[['protein']][2:394]
pos_mutation_3D <- gene_snp[['pro_mutation_count']][2:394]

#mutation position on structure
pos_mutation_t <- which(pos_mutation_3D != 0)
seq0 <- 1:length(pos_mutation_3D)
#mutation number on structure
pos_count_t <- pos_mutation_3D[pos_mutation_t]



#step 3
#wap calculation for each pair mutated residue
#calculate the standardard sample number
sample_standard <- sampleStand(pos_mutation_3D)

#calculate the wap for each pair of mutated residues based on mutation postion
wap_original <- getTotalWAP(pos_mutation_t,sample_standard,ResidueDistance_1n8p)


#change the postion of mutation while keep the mutation number in each postion
#only change the postion but not change the mutated number???
wap_sample0 <- getSampleWAP(pos_mutation_t,sample_standard,ResidueDistance_1n8p, seq=seq0,n=10000)

#analyze the result
plotNullDistribution(wap_sample0)
getPvalue(wap_original,wap_sample0)



#step 4 evaluate the each residue
#calculate the mutant number near a mutated residue
total_mutant_num <- vector()
total_mutant_position <- vector()
detail_mutant_position <- list()

for (i in 1:ncol(ResidueDistance_1n8p)) {
  col1 <- unlist(c(ResidueDistance_1n8p[, i]))
  index1 <- which(col1 <= 10 & col1 > 0)
  count0 <- pos_mutation_3D[index1]
  
  if(length(which(count0 !=0))){
    total_mutant_position[i] <- length(which(count0 !=0))
    detail_mutant_position[[i]] <- index1[count0 !=0] ##obtain the detailed position which has mutated residue
  } else{
    total_mutant_position[i] <- 0
    detail_mutant_position[[i]] <- 0
  }
  
  total_mutant_num[i] <- sum(count0)
}

mutant_1n8p <- data.frame(amino_acid = amino_acid_3D, position = 1:393, mutant_num = total_mutant,
                          mutant_position = total_mutant_position, stringsAsFactors = FALSE)


# calculate p value for each postion

for (i in seq_along(detail_mutant_position)) {
  if (length(detail_mutant_position[[i]]) >= 3) {
    pos_mutation_t0 <- detail_mutant_position[[i]]
    wap_original <- getTotalWAP(pos_mutation_t0, sample_standard, ResidueDistance_1n8p)
    wap_sample0 <- getSampleWAP(pos_mutation_t0, sample_standard, ResidueDistance_1n8p, seq = seq0, n = 10000)
    mutant_1n8p$p_value[i] <- getPvalue(wap_original, wap_sample0)
  } else {
    mutant_1n8p$p_value[i] <- 1
  }
}


#still need a fuction to filter the residue cluster that prepared for staticallly analysis






#step5 hotspot analysis
#5.1 estimate the mutation position
#mutation position on structure
pos_mutation_t <- which(pos_mutation_3D != 0)
pos_count_t <- pos_mutation_3D[pos_mutation_t]

distance <- ResidueDistance_1n8p
  all_pair <- combn(pos_mutation_t, 2)
  all_pair_distance <- vector()
  all_pair_list <- vector()
  aa_distance <- vector()
  for (i in 1:ncol(all_pair)){
    s1 <- all_pair[,i]
    d0 <- distance[s1[1],s1[2]]
    all_pair_distance[i] <- d0
    all_pair_list[i] <- paste(s1[1], s1[2], sep = "@")
    aa_distance[i] <- abs(s1[1]-s1[2])
    
    }
 

## 
pos_initial <- 1:393  
all_distance <- vector()
all_pair_ini <- combn(pos_initial, 2)
for (i in 1:ncol(all_pair_ini)){
  s1 <- all_pair_ini[,i]
  d0 <- distance[s1[1],s1[2]]
  all_distance[i] <- d0
}

##
pvalue_pair <- vector()
for (i in seq_along(all_pair_distance)){
   distance0 <- all_pair_distance[i]
   pvalue_pair[i] <- length(which(all_distance <= distance0))/length(all_distance)
}


#choose the cluste based on p value, Distance between two pair and sperated residues(>20)
target_pair <- vector()
index00 <- which(all_pair_distance <= 10 & pvalue_pair <= 0.05 & aa_distance >=20)
target_pair <- all_pair_list[index00]
target_pair_distance <- all_pair_distance[index00]
#if N > 20
target_pair0 <- str_split(target_pair, "@")


#form the unique clust based on Floyd-Warshall shortest-paths algorithm
library(igraph)
#change the pair into the link relation
links <- data_frame(from= vector(length = length(target_pair0)), to = vector(length = length(target_pair0)), weight=target_pair_distance)
for (i in seq_along(target_pair0)){
  links$from[i] <- target_pair0[[i]][1]
  links$to[i] <- target_pair0[[i]][2]
  links$weight[i] <- target_pair_distance[i]
}

g <- graph_from_data_frame(d=links,directed=FALSE) 
plot(g, edge.arrow.size=.4)


detail_mutant_position0 <- list()
dg <- decompose.graph(g)
for (i in seq_along(dg)){
  clust2 <- dg[[i]][1]
  detail_mutant_position0[[i]] <- as.integer(names(clust2))
}



pvalue <- vector()
for (i in seq_along(detail_mutant_position0)) {
  if (length(detail_mutant_position0[[i]]) >= 2) {
    pos_mutation_t0 <- detail_mutant_position0[[i]]
    wap_original <- getTotalWAP(pos_mutation_t0, sample_standard, ResidueDistance_1n8p)
    wap_sample0 <- getSampleWAP(pos_mutation_t0, sample_standard, ResidueDistance_1n8p, seq = seq0, n = 10000)
    pvalue[i] <- getPvalue(wap_original, wap_sample0)
  } else {
    pvalue[i] <- 1
  }
}

#in the analysis, we can futher calculate the closeness



