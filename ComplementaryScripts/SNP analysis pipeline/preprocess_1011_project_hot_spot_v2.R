#----------------note
#this main script is used to handle with the gene mutation only from SNP information
#in this process, the gene with SNP will be translated into protein, based on which
#the SNP could be classified into nsSNP and sSNP
#Only nsSNP is used to mapping onto protein 3D structure

source('genomics annotation summary.R')
source('getGeneCoordinate.R')
source('preprocess_1011_project_function.R')


#step0 choose samples that need to be analyzed
strain_classification <- read_excel("data/strain_classification.xls")
unique(strain_classification$Clades)

strain_select1 <- filter(strain_classification, str_detect(strain_classification$Clades,"Wine")) %>%
  select(.,Standardized_name)


#------------new version----------------------------------------------------------
#------------this version is used to preprocess data from 1011 project
# step 1 
#preprocess the SNP information
ss <- "YAL012W"
mutated_gene1 <- preprocessSNP(ss)


gene_snp <- getGeneCoordinate(gene_name = ss, genesum = gene_feature_GEM)
gene_snp[['pro_mutation_count']] <- countMutationProtein(gene_name = ss, mutation_annotation=mutated_gene1)
sum(gene_snp[['pro_mutation_count']])
gene_snp[['pro_coordinate']] <- 1:length(gene_snp[['protein']])
pos_mutation <- which(gene_snp[['pro_mutation_count']] != 0)


#step 2 input the structure information
#input the distance of all the pired residues
dirForDistanceMatrix <- paste("data/ResidueDistance_",ss,".txt", sep = "")
ResidueDistance <- read.table(dirForDistanceMatrix,sep = ",") #in the followed calculation, the matrix dosen't have the col and row names
ResidueDistance <- as.matrix(ResidueDistance)

#the amino acid sequence in structure is from 2:394 while  the original sequence is from 1:394
#obtain the mutation information for the structure
seq_3D_origin <- 2:394#seq_from_3D <- 2:394 #"YAL012W.fasta"#this is the coordinated of original protein sequence and should changed into 3D structure coordinates
amino_acid_3D <- gene_snp[['protein']][seq_3D_origin]
count_mutation_3D <- gene_snp[['pro_mutation_count']][seq_3D_origin]

#mutation position on structure and #mutation number on structure
pos_mutation_3D <- which(count_mutation_3D != 0)
seq_3D <- 1:length(count_mutation_3D) #seq0 is the coordinate of PDB structure
mutation_count_3D <- count_mutation_3D[pos_mutation_3D]

#wap calculation for each pair mutated residue
#calculate the standardard sample number
sample_standard1 <- sampleStand(count_mutation_3D)


#step 3 hot spot analysis
#this main function will be divided into different parts for easy understanding
pos_residue1 <- list()
for (i in seq_along(mutated_gene1$strain)){
  pos_residue1[[i]] <- PositionResidueSNP(mutated_gene1$Pos[i],mutated_gene1$Alt[i], ss)
}


pos_residue_df <- ResidueSum(pos_residue1)


#mapping the mutate residue onto the original protein sequence
gene_snp[['residue']] <- getMultipleReactionFormula(pos_residue_df$residue, pos_residue_df$pos,gene_snp[['pro_coordinate']]) 
residue_3D <- gene_snp[['residue']][seq_3D_origin]

#obtain the paired residue 
residue_pair <- getHotVertice(aa_3d = seq_3D, residue0 = residue_3D, aa_pro = seq_3D_origin, distance0 = ResidueDistance)

#calculate closeness of each cluster
important_hot <- clusterAnalysis(residue_pair)


#further calculate the p value of each choosed cluster
#calculate the pvalue for each clust
important_hot$pvalue <- getHotPvalue(cluster0 = important_hot$cluster, 
                                     sample_standard = sample_standard1, 
                                     distance=ResidueDistance, 
                                     seq = seq_3D)




