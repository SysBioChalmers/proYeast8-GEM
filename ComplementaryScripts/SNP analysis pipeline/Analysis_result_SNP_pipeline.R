
#----------------note
# this main script is used to handle with the gene mutation only from SNP information
# in this process, the gene with SNP will be translated into protein, based on which
# the SNP could be classified into nsSNP and sSNP
# Only nsSNP is used to mapping onto protein 3D structure

source("genomics annotation summary.R")
source("preprocess_1011_project_function.R")
source("parse_active_site_for_hotspot_analysis.R")
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
strain_select1 <- chooseStrain(type = strain_type)





# input the data of clumps method
clumps_ex <- read.table('result0/CLUMPS from pdb_ex for all_strain/pdb_EX.txt', header= TRUE, sep = "\t", stringsAsFactors = FALSE)
clumps_ex$pdb_source <- "Exp"
clumps_homo <- read.table('result0/CLUMPS from pdb_homo for all_strain/pdb_info.txt', header= TRUE, sep = "\t", stringsAsFactors = FALSE)
clumps_homo$pdb_source <- "Homo"

clumps_all <- rbind.data.frame(clumps_ex, clumps_homo)
clumps_all_fiter <- filter(clumps_all, p_value <= 0.05)

# remove pdb file with the same gene and same residue coordinate
clumps_all_fiter$combineINF <- paste(clumps_all_fiter$locus, clumps_all_fiter$sstart2, clumps_all_fiter$send2, sep = "@@")
clumps_all_fiter0 <- clumps_all_fiter[!duplicated(clumps_all_fiter$combineINF),]


for (i in 37:38) {
  print(i)
  # step 1
  # preprocess the SNP information
  ss <- clumps_all_fiter0$locus[i]
  mutated_gene0 <- preprocessSNP(ss)
  mutated_gene1 <- mutated_gene0[which(mutated_gene0$strain %in% strain_select1$Standardized_name), ]

  # input the protein coordinate ID
  p1 <- clumps_all_fiter0$sstart2[i]
  p2 <- clumps_all_fiter0$send2[i]
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
  result$pdbID <- clumps_all_fiter0$pdbid[i]
  result <- select(result, orf, pdbID, position, ref, alt, Freq)

}

#note
# YMR207C ref exist error which is due to the protein seq is not consitant between uniprot and sgd
        