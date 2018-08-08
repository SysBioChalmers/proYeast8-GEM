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
pos_mutation <- which(gene_snp[['pro_mutation_count']] != 0)


#step 2 input the structure information
# mapping the mutation postion onto the protein structure
#input the distance of all the pired residues
#ResidueDistance_1n8p <- read_excel("data/ResidueDistance_YAL012W.xlsx",col_names = FALSE) #in the followed calculation, the matrix dosen't have the col and row names
dirForDistanceMatrix <- paste("data/ResidueDistance_",ss,".txt", sep = "")
ResidueDistance <- read.table(dirForDistanceMatrix,sep = ",") #in the followed calculation, the matrix dosen't have the col and row names
ResidueDistance <- as.matrix(ResidueDistance)

#the amino acid sequence in structure is from 2:394 while  the original sequence is from 1:394
#obtain the mutation information for the structure
seq_3D_origin <- 2:394 #"YAL012W.fasta"#this is the coordinated of original protein sequence and should changed into 3D structure coordinates
amino_acid_3D <- gene_snp[['protein']][seq_3D_origin]
count_mutation_3D <- gene_snp[['pro_mutation_count']][seq_3D_origin]

#mutation position on structure and #mutation number on structure
pos_mutation_3D <- which(count_mutation_3D != 0)
seq_3D <- 1:length(count_mutation_3D) #seq_3D is the coordinate of PDB structure
mutation_count_3D <- count_mutation_3D[pos_mutation_c] # this parameter is not used



#step 3
#wap calculation for each pair mutated residue
#calculate the standardard sample number
sample_standard1 <- sampleStand(count_mutation_3D)

#calculate the wap for each pair of mutated residues based on mutation postion
wap_original <- getTotalWAP(pos_mutation_3D,sample_standard1,ResidueDistance)

#change the postion of mutation while keep the mutation number in each postion
#only change the postion but not change the mutated number???
wap_sample0 <- getSampleWAP(pos_mutation_3D,sample_standard1,ResidueDistance, seq=seq_3D, n=10000)

#analyze the result
plotNullDistribution(wap_sample0)
Strain_3D <- getPvalue(wap_original,wap_sample0)






##batch process for the above whole process
#input we need
# gene id
# snp
# gene_feature_GEM
# coordinate of pdb on original sequence
# coordinate of pdb itself
# residue distance matrix


#output
#p_value with id

























