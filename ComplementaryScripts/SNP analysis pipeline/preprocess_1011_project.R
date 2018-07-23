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
gene_snp[['pro_mutation_count']] <- countMutationProtein0(geneName  = ss, mutated_gene_seq = df_list_filter)



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
pos_mutation_c <- which(pos_mutation_3D != 0)
seq0 <- 1:length(pos_mutation_3D)
#mutation number on structure
pos_count_num <- pos_mutation_3D[pos_mutation_c]



#step 3
#wap calculation for each pair mutated residue
#calculate the standardard sample number
sample_standard1 <- sampleStand(pos_mutation_3D)

#calculate the wap for each pair of mutated residues based on mutation postion
wap_original <- getTotalWAP(pos_mutation_c,sample_standard1,ResidueDistance_1n8p)

#change the postion of mutation while keep the mutation number in each postion
#only change the postion but not change the mutated number???
wap_sample0 <- getSampleWAP(pos_mutation_c,sample_standard1,ResidueDistance_1n8p, seq=seq0,n=10000)

#analyze the result
plotNullDistribution(wap_sample0)
getPvalue(wap_original,wap_sample0)


#step4 hotspot analysis
hotSpotAnalysis(pos_mutation_t = pos_mutation_c,
                sample_standard = sample_standard1,
                distance = ResidueDistance_1n8p,
                pos_initial = seq0
                )








