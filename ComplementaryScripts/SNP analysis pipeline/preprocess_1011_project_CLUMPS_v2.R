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
# for 'YAL012W', the PDB structure id is 'https://www.rcsb.org/structure/1n8p'

#input the distance of all the pired residues
ResidueDistance_1n8p <- read_excel("data/ResidueDistance_YAL012W.xlsx",col_names = FALSE) #in the followed calculation, the matrix dosen't have the col and row names
ResidueDistance_1n8p <- as.matrix(ResidueDistance_1n8p)

#the amino acid sequence in structure is from 2:394 while  the original sequence is from 1:394
#obtain the mutation information for the structure
seq_from_3D <- 2:394 #"YAL012W.fasta"#this is the coordinated of original protein sequence and should changed into 3D structure coordinates
amino_acid_3D <- gene_snp[['protein']][seq_from_3D]
count_mutation_3D <- gene_snp[['pro_mutation_count']][seq_from_3D]

#mutation position on structure and #mutation number on structure
pos_mutation_c <- which(count_mutation_3D != 0)
seq0 <- 1:length(count_mutation_3D) #seq0 is the coordinate of PDB structure
pos_count_num <- count_mutation_3D[pos_mutation_c]



#step 3
#wap calculation for each pair mutated residue
#calculate the standardard sample number
sample_standard1 <- sampleStand(count_mutation_3D)

#calculate the wap for each pair of mutated residues based on mutation postion
wap_original <- getTotalWAP(pos_mutation_c,sample_standard1,ResidueDistance_1n8p)

#change the postion of mutation while keep the mutation number in each postion
#only change the postion but not change the mutated number???
wap_sample0 <- getSampleWAP(pos_mutation_c,sample_standard1,ResidueDistance_1n8p, seq=seq0,n=10000)

#analyze the result
plotNullDistribution(wap_sample0)
Strain_3D <- getPvalue(wap_original,wap_sample0)


























