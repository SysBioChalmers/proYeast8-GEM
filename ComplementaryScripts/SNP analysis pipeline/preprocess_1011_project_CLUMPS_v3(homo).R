#----------------note
#this main script is used to handle with the gene mutation only from SNP information
#in this process, the gene with SNP will be translated into protein, based on which
#the SNP could be classified into nsSNP and sSNP
#Only nsSNP is used to mapping onto protein 3D structure
source('genomics annotation summary.R')
source('preprocess_1011_project_function.R')


#step0 choose samples that need to be analyzed
strain_classification <- read_excel("data/strain_classification.xls")
unique(strain_classification$Clades)

strain_type <- "bioethanol"
strain_select1 <- filter(strain_classification, str_detect(strain_classification$Clades, strain_type)) %>%
  select(., Standardized_name)

#------------new version----------------------------------------------------------
#------------this version is used to preprocess data from 1011 project
# step 1 
#preprocess the SNP information
ss <- "YPR110C"
mutated_gene0 <- preprocessSNP(ss)
mutated_gene1 <- mutated_gene0[which(mutated_gene0$strain %in% strain_select1$Standardized_name), ]

gene_snp <- getGeneCoordinate(gene_name = ss, genesum = gene_feature_GEM)
gene_snp[["pro_mutation_count"]] <- countMutationProtein(gene_name = ss, mutation_annotation = mutated_gene1)
gene_snp[['pro_mutation_count']] <- countMutationProtein(gene_name = ss, mutation_annotation=mutated_gene1)
pos_mutation <- which(gene_snp[['pro_mutation_count']] != 0)


#step 2 input the structure information
# mapping the mutation postion onto the protein structure
#input the distance of all the pired residues
#ResidueDistance_1n8p <- read_excel("data/ResidueDistance_YAL012W.xlsx",col_names = FALSE) #in the followed calculation, the matrix dosen't have the col and row names
pdbID <- "4c2m@C"
r1 <- 1
r2 <- 305 # input the corrected residue sequnence
r3 <- paste(r1, r2, sep = "-")
dirForDistanceMatrix <- paste("residue_distance/", pdbID, ".txt", sep = "")
ResidueDistance0 <- read.table(dirForDistanceMatrix, sep = ",") # in the followed calculation, the matrix dosen't have the col and row names
ResidueDistance0 <- as.matrix(ResidueDistance0)
ResidueDistance <- ResidueDistance0[r1:r2, r1:r2]


#the amino acid sequence in structure is from 2:394 while  the original sequence is from 1:394
#obtain the mutation information for the structure
p1 <- 31
p2 <- 335
p3 <- paste(p1, p2, sep = "-")
seq_3D_origin <- p1:p2 #"YAL012W.fasta"#this is the coordinated of original protein sequence and should changed into 3D structure coordinates
amino_acid_3D <- gene_snp[['protein']][seq_3D_origin]
count_mutation_3D <- gene_snp[['pro_mutation_count']][seq_3D_origin]

#mutation position on structure and #mutation number on structure
pos_mutation_3D <- which(count_mutation_3D != 0)
seq_3D <- 1:length(count_mutation_3D) #seq_3D is the coordinate of PDB structure
mutation_count_3D <- count_mutation_3D[pos_mutation_3D] # this parameter is not used



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
#------------batch process--------------------------------------------------------
#------------new version----------------------------------------------------------
#------------this version is used to preprocess data from 1011 project

# step 0
# input the gene information
pdb_info  <- read_excel("data/pdb_homo summary for manual check.xlsx")
pdb_info$pdbid <- paste(pdb_info$sstart2, pdb_info$send2, pdb_info$template, pdb_info$coordinate_id,sep = "_")
pdb_info <- filter(pdb_info, is.na(pdb_info$with_distance))
pdb_info <- select(pdb_info, locus, pdbid, qstart2, qend2, sstart2, send2)
geneWithSNP <- getGeneNameWithSNP()
pdb_info <- pdb_info[which(pdb_info$locus %in% geneWithSNP ==TRUE),]
pdb_info$pdbid[506:513] <- paste('s', 5:12, sep = "") #chainid for some homo-pdb

#add two more clumns
pdb_info$strain_type <- strain_type
pdb_info$p_value <- NA


#creat new file to store the results
outfile0 <- paste('result/CLUMPS from pdb_homo for ', strain_type, sep = "")
dir.create(outfile0)
print(outfile0)


for (i in 1:723) {
  # step 1
  # preprocess the SNP information
  #i <- 210
  print(i)
  ss <- pdb_info$locus[i]
  mutated_gene0 <- preprocessSNP(ss)
  mutated_gene1 <- mutated_gene0[which(mutated_gene0$strain %in% strain_select1$Standardized_name), ]

  gene_snp <- getGeneCoordinate(gene_name = ss, genesum = gene_feature_GEM)
  gene_snp[["pro_mutation_count"]] <- countMutationProtein(gene_name = ss, mutation_annotation = mutated_gene1)
  gene_snp[["pro_coordinate"]] <- 1:length(gene_snp[["protein"]])
  pos_mutation <- which(gene_snp[["pro_mutation_count"]] != 0)


  # step 2 input the structure information
  # input the distance of all the pired residues
  pdbID <- pdb_info$pdbid[i]
  r1 <- pdb_info$qstart2[i]
  r2 <- pdb_info$qend2[i] # input the corrected residue sequnence
  r3 <- paste(r1, r2, sep = "-")
  dirForDistanceMatrix <- paste("residue_distance/pdb_homo/", pdbID, ".pdb.txt", sep = "")
  ResidueDistance0 <- read.table(dirForDistanceMatrix, sep = ",") # in the followed calculation, the matrix dosen't have the col and row names
  ResidueDistance0 <- as.matrix(ResidueDistance0)
  ResidueDistance <- ResidueDistance0 # [r1:r2,r1:r2]
  dim(ResidueDistance)

  # the amino acid sequence in structure is from 2:394 while  the original sequence is from 1:394
  # obtain the mutation information for the structure
  p1 <- pdb_info$sstart2[i]
  p2 <- pdb_info$send2[i]
  p3 <- paste(p1, p2, sep = "-")
  seq_3D_origin <- p1:p2 # seq_from_3D <- 2:394 #"YAL012W.fasta"#this is the coordinated of original protein sequence and should changed into 3D structure coordinates
  amino_acid_3D <- gene_snp[["protein"]][seq_3D_origin]
  count_mutation_3D <- gene_snp[["pro_mutation_count"]][seq_3D_origin]

  # mutation position on structure and #mutation number on structure
  pos_mutation_3D <- which(count_mutation_3D != 0)
  seq_3D <- 1:length(count_mutation_3D) # seq0 is the coordinate of PDB structure
  mutation_count_3D <- count_mutation_3D[pos_mutation_3D]

  # there should be two postions which have mutations
  if (length(pos_mutation_3D) >= 2) {
    # wap calculation for each pair mutated residue
    # calculate the standardard sample number
    sample_standard1 <- sampleStand(count_mutation_3D)

    # step 3
    # calculate the standardard sample number
    sample_standard1 <- sampleStand(count_mutation_3D)

    # calculate the wap for each pair of mutated residues based on mutation postion
    wap_original <- getTotalWAP(pos_mutation_3D, sample_standard1, ResidueDistance)

    # change the postion of mutation while keep the mutation number in each postion
    # only change the postion but not change the mutated number???
    wap_sample0 <- getSampleWAP(pos_mutation_3D, sample_standard1, ResidueDistance, seq = seq_3D, n = 10000)

    # analyze the result
    plotNullDistribution(wap_sample0)
    pdb_info$p_value[i] <- getPvalue(wap_original, wap_sample0)
    print(paste('-------p_value=', pdb_info$p_value[i], sep = ""))
  } else {
    pdb_info$p_value[i] <- 'NA'
    print("Not enough mutation")
    next
  }

}


# save the result
write.table(pdb_info, paste(outfile0,'/','pdb_info.txt', sep = ""), row.names = FALSE, sep = "\t")


