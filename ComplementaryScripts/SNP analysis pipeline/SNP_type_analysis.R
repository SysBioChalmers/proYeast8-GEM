#----------------note
# this main script is used to calculate nsSNP number occured for each gene based on SNP number
source("genomics annotation summary.R")
source("preprocess_1011_project_function.R")


# step0 choose samples that need to be analyzed
strain_classification <- read_excel("data/strain_classification.xls")
strain_type <- "all_strain"

geneGEM_with_SNP_number <- read_csv("data/geneGEM with SNP number.csv")


for (i in 1:nrow(geneGEM_with_SNP_number)) {
  print(i)
  ss <- geneGEM_with_SNP_number$geneNames[i]
  strain_select1 <- chooseStrain(type = strain_type)
  
  if(geneGEM_with_SNP_number$SNP_NUM[i] > 0){
  mutated_gene0 <- preprocessSNP(ss)
  mutated_gene1 <- mutated_gene0[which(mutated_gene0$strain %in% strain_select1$Standardized_name), ]

  # produce the mutaed residue and its position
  nsSNP_num <- 0
  for (j in seq_along(mutated_gene1$strain)) {
    s <- PositionResidueSNP(mutated_gene1$Pos[j], mutated_gene1$Alt[j], ss)

    if (length(s) != 0) {
      nsSNP_num <- nsSNP_num + 1
    }
  }
  geneGEM_with_SNP_number$nsSNP[i] <- nsSNP_num
  } else{
    geneGEM_with_SNP_number$nsSNP[i]
  }
}






