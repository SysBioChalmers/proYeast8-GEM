#----------------note
# this main script is used to handle with the gene mutation only from SNP information
# in this process, the gene with SNP will be translated into protein, based on which
# the SNP could be classified into nsSNP and sSNP
# Only nsSNP is used to mapping onto protein 3D structure
# the nsSNP will be printed used for the mutation function prediction based on mutFunc
source("preprocess_1011_project_function.R")
source("genomics annotation summary.R")


# function to obtain the interesting protein list from clumps analysis result
# here we only choose the proteinID with p_value <= 0.05
# part 0
substrate0 <- 'glycerol'
strain_type0 <- 'glycerol_high'
infile0 = paste("data/strain_", substrate0,"_classification.txt", sep = "")
strain_classification <- read.table(infile0, header = TRUE, stringsAsFactors = FALSE)
strain_classification <- strain_classification[, c("strain_name", "type")]



# input the gene annotation data
gene_annotation <- read.delim2('data/all_gene_yeast with annotation from different database.txt', header = TRUE, stringsAsFactors = FALSE)
uniprotGeneID_mapping <- read_excel("data/uniprotGeneID_mapping.xlsx")
# input the result for the kinetics sensitivity analysis
kinetics_analysis <- read.table("/Users/luho/Documents/GitHub/ecModels/ecYeastGEM/figures/yeast8/results/KcatSensitivities_YEP.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(kinetics_analysis) <- c("Protein", "glucose", "acetate", "ethanol", "glycerol", "sorbitol", "galactose", "ribose", "xylose")
# id-mapping
kinetics_analysis$locus <- getSingleReactionFormula(uniprotGeneID_mapping$GeneName,uniprotGeneID_mapping$Entry,kinetics_analysis$Protein)
# density graph
plot(density(kinetics_analysis$glycerol))



# part 1
hotspot_glycerol  <- read.table('data/hotspot analysis for glycerol.txt', header = TRUE, stringsAsFactors = FALSE)




# part 2
# print the mutated residue information for the strain classification
geneset1 <- 'YJL052W'
geneset1_index <- which(hotspot_glycerol$gene %in% geneset1)
geneset1_detail <- hotspot_glycerol[geneset1_index,]
geneset1_detail <- geneset1_detail %>% separate(.,seq_3D_origin, into = c('sstart2', 'send2'), sep = "-")

# function to choose the strain list with specific phenotype
strainList <- function(substrate, strain_type) {
  #input
  #substrate <-'glycerol'
  #strain_type: glycerol_high, glycerol_low, glycerol_medium. which represent three different types of 
  #strains with high, medium and low growth rate on the medium with glycerol as the substrate
  infile0 = paste("data/strain_", substrate,"_classification.txt", sep = "")
  strain_classification <- read.table(infile0, header = TRUE, stringsAsFactors = FALSE)
  strain_classification <- strain_classification[, c("strain_name", "type")]
  # strain_type <-"glycerol_high"
  strain_select0 <- chooseStrain(type = strain_type, strain0=strain_classification)
  return( strain_select0)
}

# prepare the strain
strain_select11 <- strainList(substrate=substrate0, strain_type=strain_type0)

for (i in 1:length(geneset1_detail$locus)) {
  print(i)
  # step 1
  # preprocess the SNP information
  i=1
  ss <- geneset1_detail$gene[i]
  mutated_gene0 <- preprocessSNP(ss)
  mutated_gene1 <- mutated_gene0[which(mutated_gene0$strain %in% strain_select11$Standardized_name), ]
  # input the protein coordinate ID
  p1 <- as.numeric(geneset1_detail$sstart2[i])
  p2 <- as.numeric(geneset1_detail$send2[i])
  p3 <- paste(p1, p2, sep = "-")
  seq_3D_origin <- p1:p2 # this is the coordinated of original protein sequence and should changed into 3D structure coordinates

  # produce the mutaed residue and its position
  pos_residue1 <- list()
  for (j in seq_along(mutated_gene1$strain)) {
    pos_residue1[[j]] <- PositionResidueSNP(mutated_gene1$Pos[j], mutated_gene1$Alt[j], ss)
  }

  pos_residue_df <- ResidueSum(pos_residue1)

  # mapping the mutate residue onto the original protein sequence
  gene_snp <- getGeneCoordinate(gene_name = ss, genesum = gene_feature_GEM)
  gene_snp[["pro_coordinate"]] <- 1:length(gene_snp[["protein"]])
  gene_snp[["residue"]] <- getMultipleReactionFormula(pos_residue_df$residue, pos_residue_df$pos, gene_snp[["pro_coordinate"]])
  residue_3D <- gene_snp[["residue"]][seq_3D_origin]

  # analyse the mutation frequence of each residue
  residue_3D0 <- residue_3D[!is.na(residue_3D)]
  residue_3D0 <- str_split(residue_3D0, ";")
  residue_3D0 <- unlist(residue_3D0)
  unique(residue_3D0)
  tmp <- table(residue_3D0)
  result <- as.data.frame(tmp)
  result <- result %>% separate(., residue_3D0, into = c("alt", "position"), sep = "@@")
  result$position <- as.numeric(result$position)
  result <- result %>% arrange(., position)
  # obtain the original ref residue
  gene_snp <- getGeneCoordinate(gene_name = ss, genesum = gene_feature_GEM)
  result$ref <- getSingleReactionFormula(gene_snp[["protein"]], gene_snp[["protein_coordinate"]], result$position)
  result$orf <- ss
  result$pdbID <- geneset1_detail$structure[i]
  result <- select(result, orf, ref, position, alt, Freq, pdbID)
  print(result)
  write.table(result, "result/gene_mutation_related_to_growth.txt", row.names = FALSE, sep = " ")
  
}










# save the results in format which can predict the protein function
# result0 <- select(result, orf, ref, position, alt)
# result0 <- result0 %>% unite(mutation,  c('ref','position','alt'), sep="")

# write.table(result0, "result/print the snp for mutation effects prediction.txt", row.names = FALSE, sep = " ")

#note
# YMR207C ref exist error which is due to the protein seq is not consitant between uniprot and sgd
        