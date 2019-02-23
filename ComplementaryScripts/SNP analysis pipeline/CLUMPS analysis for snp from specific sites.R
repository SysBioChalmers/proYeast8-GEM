#----------------note
#this main script is used to conduct CLUMPS analysis for mutations from specific sites.
#this step can be very interesting to enlarge the hotspot zones obtained from hotspot analysis pipeline.
source('genomics annotation summary.R')
source('preprocess_1011_project_function.R')


#step0 choose samples that need to be analyzed
strain_classification <- read.table("data/strain_glycerol_classification.txt", header = TRUE)
strain_classification <- strain_classification[, c('strain_name','type')]
strain_type <-"glycerol_high"
strain_select1 <- chooseStrain(type = "glycerol_high")


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

which(pdb_info$locus=='YJL052W')




# preprocess the SNP information
i <- 172 # for YJL052W
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
# seq_3D_origin <- c(31, 73)
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
} else {
  print("Not enough mutation")
}



#method to calculate the p value for the choosed cluster
all_cluster <- vector()
all_cluster[1] <- "V@@31;A@@73"
all_cluster[2] <- "V@@31;A@@73;K@@24"
all_cluster[3] <- "V@@31;A@@73;V@@70"
all_cluster[4] <- "V@@31;A@@73;K@@24;V@@70"
all_cluster[5] <- "V@@31;A@@73;K@@24;V@@70;S@@125"
all_cluster[6] <- "V@@31;A@@73;K@@24;V@@70;S@@125;E@@248"
all_cluster[7] <- "K@@24;V@@70"
all_cluster[8] <- "S@@125;E@@248"

pvalue_sum <- vector()

for(i in 1:length(all_cluster)){
  print(i)
  pvalue_sum[i] <- getHotPvalue(
  cluster0 = all_cluster[i],
  sample_standard = sample_standard1,
  distance = ResidueDistance,
  seq = seq_3D_origin
)

}

cluster_pvalue <- data.frame(cluster=all_cluster, p_value=pvalue_sum, stringsAsFactors = FALSE)





