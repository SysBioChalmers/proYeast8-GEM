#step 1 fasta file preparation
#calculation of the mutated residue in each position

#choose the metabolic gene file
filterMetabolicGene() # a target gene file will occured after this file

#in this project, a fasta file could contained different strains, so we need a function to refine each fasta file
#so that this file just contained the strain or mutation we are interested
strain_classification <- read_excel("data/strain_classification.xls")
unique(strain_classification$Clades)

#run the analysis pipeline using an example
#read each fasta file and change it into a dataframe

#if no filter, we just preprocess the fasta file with all strain
#obtain the sequence for the gene from each strain
gene_test1 <- "YAL012W.fasta"
df_list <- processFasta(gene_test = gene_test1) 


#then the filter can be used for each dataframe to obtain the strains we need
#choose the strain in the clades of bioethanol
strain_select1 <- filter(strain_classification, str_detect(strain_classification$Clades,"Wine")) %>%
  select(.,Standardized_name)

df_list_filter <- filterMutationStrainType(gene_test = gene_test1, strain_select = strain_select1)


#count the mutant number in each protein position
ss <- str_replace_all(gene_test1, ".fasta","")
gene_snp <- getGeneCoordinate(gene_name = ss, genesum = gene_feature_GEM)
gene_snp[['pro_mutation_count']] <- countMutationProtein0(geneName  = ss, mutated_gene_seq = df_list)
gene_snp[['pro_coordinate']] <- 1:length(gene_snp[['protein']])


#step 2 input the structure information
# mapping the mutation postion onto the protein structure
# for 'YAL012W', the PDB structure id is 'https://www.rcsb.org/structure/1n8p'

#input the distance of all the pired residues
ResidueDistance_1n8p <- read_excel("data/ResidueDistance_YAL012W.xlsx",col_names = FALSE) #in the followed calculation, the matrix dosen't have the col and row names
ResidueDistance_1n8p <- as.matrix(ResidueDistance_1n8p)

#the amino acid sequence in structure is from 2:394 while  the original sequence is from 1:394
#obtain the mutation information for the structure
seq_from_3D <- 2:394#seq_from_3D <- 2:394 #"YAL012W.fasta"#this is the coordinated of original protein sequence and should changed into 3D structure coordinates
amino_acid_3D <- gene_snp[['protein']][seq_from_3D]
count_mutation_3D <- gene_snp[['pro_mutation_count']][seq_from_3D]

#mutation position on structure and #mutation number on structure
pos_mutation_c <- which(count_mutation_3D != 0)
seq0 <- 1:length(count_mutation_3D) #seq0 is the coordinate of PDB structure
pos_count_num <- count_mutation_3D[pos_mutation_c]
sample_standard1 <- sampleStand(count_mutation_3D)


#step 3 hot spot analysis
#this main function will be divided into different parts for easy understanding

# with the protein coordinate-mutated amino acid
pos_residue1 <- list()
for (i in seq_along(df_list)){
  pos_residue1[[i]] <- PositionResidue(df_list[[i]], ss)
}

pos_residue_df <- ResidueSum(pos_residue1)  

#mapping the mutate residue onto the original protein sequence
gene_snp[['residue']] <- getMultipleReactionFormula(pos_residue_df$residue, pos_residue_df$pos,gene_snp[['pro_coordinate']]) 
residue_3D <- gene_snp[['residue']][seq_from_3D]

#obtain the paired residue 
residue_pair <- getHotVertice(aa_3d = seq0, residue0 = residue_3D, aa_pro = seq_from_3D, distance0 = ResidueDistance_1n8p)

#calculate closeness of each cluster
important_hot <- clusterAnalysis(residue_pair)


#further calculate the p value of each choosed cluster
# calculate the pvalue for each clust
important_hot$pvalue <- getHotPvalue(cluster0 = important_hot$cluster, 
                                     sample_standard = sample_standard1, 
                                     distance=ResidueDistance_1n8p, 
                                     seq = seq0)




