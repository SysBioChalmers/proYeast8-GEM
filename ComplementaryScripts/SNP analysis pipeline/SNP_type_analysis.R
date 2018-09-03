#----------------note
# this main script is used to handle with the gene mutation only from SNP information
# in this process, the gene with SNP will be translated into protein, based on which
# the SNP could be classified into nsSNP and sSNP
# Only nsSNP is used to mapping onto protein 3D structure
source("genomics annotation summary.R")
source("preprocess_1011_project_function.R")

# newly added function, will integrate the main function
chooseStrain <- function(type,strain0=strain_classification){
  if(type=="all_strain"){
    return(strain0)
  } else{
    strain_select <- filter(strain_classification, str_detect(strain_classification$Clades, strain_type)) %>%
      select(., Standardized_name)
    return(strain_select)
  }
  
}

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
#note
# YMR207C ref exist error which is due to the protein seq is not consitant between uniprot and sgd
        