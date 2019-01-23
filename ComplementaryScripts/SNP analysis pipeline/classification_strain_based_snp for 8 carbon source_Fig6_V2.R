#----------------note
# this main script is used to classify the strain based on residues mutation information obtained
# from CLUMPS analysis pipeline.
# firstly we obtain the residues mutation of one gene for all strains
# then we can classify the strains based on the specific residue mutations
source("Analysis_result_SNP_pipeline_hotspot.R")


# function to preprocess the SNP information
# this function will obtain all the mutated residue information for one protein
printAllMutationWithGene <- function(ss){
  #ss <- 'YMR246W'
  mutated_gene0 <- preprocessSNP(ss)
  mutated_gene1 <- mutated_gene0[which(mutated_gene0$strain %in% strain_select1$Standardized_name), ]
  gene_snp <- getGeneCoordinate(gene_name = ss, genesum = gene_feature_GEM)
  # input the protein coordinate ID
  p1 <- 1
  p2 <- length(gene_snp[['protein']])
  p3 <- paste(p1, p2, sep = "-")
  # produce the mutaed residue and its position
  pos_residue1 <- list()
  for (j in seq_along(mutated_gene1$strain)) {
    print(j)
    pos_residue1[[j]] <- PositionResidueSNP(mutated_gene1$Pos[j], mutated_gene1$Alt[j], ss)
    mutated_gene1$mutation_result[j] <- pos_residue1[[j]][1]
  }
  
  result <- mutated_gene1[!is.na(mutated_gene1$mutation_result),]
  result$gene_mutation <- paste(result$Gene2,result$mutation_result, sep = '@@')
  result0 <- select(result,strain,gene_mutation)
  return(result0)
}



# part 3
# datainput
# choose samples that need to be analyzed
growth_phenotype <- read_excel("data/phenoMatrix_35ConditionsNormalizedByYPD.xlsx") 
growth_phenotype <- select(growth_phenotype, strain_name, YPGLYCEROL)
colnames(growth_phenotype) <- c('strain_name','growth')

# here we choose strain with the related specific phenotype but remove the strains
# which have been used for finding interesting proteins list
infile1 = paste("data/strain_", substrate0,"_classification.txt", sep = "")
strain_choosed <- read.table(infile1, header = TRUE, stringsAsFactors = FALSE)
strain_classification <- read_excel("data/strain_classification.xls")
strain_select1 <- strain_classification[which(strain_classification$Standardized_name %in% strain_choosed$strain_name ==FALSE),]


# part 4
# choose the interested protein and the related mutation datasets
gene_mutation  <- read.table("result/gene_mutation_related_to_growth_hotspot.txt", header = TRUE, stringsAsFactors = FALSE)
gene_mutation0 <- gene_mutation  %>% unite(., mutation_result, orf, alt, position, sep = "@@") %>%
   filter(.,Freq >=1)
gene_mutation1 <- gene_mutation0$mutation_result
ll0 <- length(gene_mutation1)

# produce the gene mutation matrix for all the strains
# YHL032C
p <- printAllMutationWithGene(gene_mutation$orf[1])
# YAR035W 
#p_YAR035W <- printAllMutationWithGene('YAR035W')
# YML100W
#p_YML100W <- printAllMutationWithGene('YML100W')
# YJL068C
#p_YJL068C <- printAllMutationWithGene('YJL068C')

#merge the above data
strain_mutation <- rbind.data.frame(p)
mutation_summary <- data.frame(strain_name=unique(strain_mutation$strain), stringsAsFactors = FALSE) 
mutation_summary$mutation_result <- getMultipleReactionFormula(strain_mutation$gene_mutation,strain_mutation$strain,mutation_summary$strain_name)


#establish a list contains the strain name and mutation 
#here we should produce a matrix to check whether each mutation contains in each strain,
#function to check whether a mutation existed
checkExist <- function(mutate00, mutation_list=mutation_summary){
s1 <- str_split(mutation_list$mutation_result,';')
exist0 <- vector()
for(i in seq_along(s1)){
  s11 <- s1[[i]]
   if(mutate00 %in% s11){
      exist0[i] <- 1
   } else {
     exist0[i] <- 0
  }
 }
  return(exist0)
} 


for(i in gene_mutation1){
  print(i)
  checkExist(i)
  mutation_summary[,i] <- checkExist(i, mutation_list=mutation_summary)
}

mutation_matrix <- select(mutation_summary, c('strain_name',gene_mutation1 ))
#strain classification based on several genes with a series mutations
#here we can classify the strain based on cluster 
#and the total numer of mutation occured in each strain


#method2 cluster the strain based on the mutation number
##mutation number analysis
mutation_matrix0 <- mutation_matrix
mutation_matrix0$sum <- rowSums(mutation_matrix0[,2:ncol(mutation_matrix0)]) #calculate the number of aimed SNP
table(mutation_matrix0$sum)
plot(density(mutation_matrix0$sum))

# part5
# classification based on number
strain1 <- mutation_matrix0[mutation_matrix0$sum >= 1,] %>% select(.,strain_name)
cluster1 <-  growth_phenotype[growth_phenotype$strain_name %in% strain1$strain_name,]
cluster1$type <- 'g1a'

strain2 <- mutation_matrix0[mutation_matrix0$sum >= 0 & mutation_matrix0$sum <= 0,] %>% select(.,strain_name)
cluster2 <-  growth_phenotype[growth_phenotype$strain_name %in% strain2$strain_name,]
cluster2$type <- 'g2a'



# strain classification based on the mutation from non-hotspot zone
# choose the interested protein and the related mutation datasets
gene_mutation  <- read.table("result/gene_mutation_related_to_growth_nonhotspot.txt", header = TRUE, stringsAsFactors = FALSE)
gene_mutation0 <- gene_mutation  %>% unite(., mutation_result, orf, alt, position, sep = "@@") %>%
  filter(.,Freq >=1)
gene_mutation1 <- gene_mutation0$mutation_result
ll0 <- length(gene_mutation1)

# produce the gene mutation matrix for all the strains
# YHL032C
p <- printAllMutationWithGene(gene_mutation$orf[1])
# YAR035W 
#p_YAR035W <- printAllMutationWithGene('YAR035W')
# YML100W
#p_YML100W <- printAllMutationWithGene('YML100W')
# YJL068C
#p_YJL068C <- printAllMutationWithGene('YJL068C')

#merge the above data
strain_mutation <- rbind.data.frame(p)
mutation_summary <- data.frame(strain_name=unique(strain_mutation$strain), stringsAsFactors = FALSE) 
mutation_summary$mutation_result <- getMultipleReactionFormula(strain_mutation$gene_mutation,strain_mutation$strain,mutation_summary$strain_name)


#establish a list contains the strain name and mutation 
#here we should produce a matrix to check whether each mutation contains in each strain,
#function to check whether a mutation existed
checkExist <- function(mutate00, mutation_list=mutation_summary){
  s1 <- str_split(mutation_list$mutation_result,';')
  exist0 <- vector()
  for(i in seq_along(s1)){
    s11 <- s1[[i]]
    if(mutate00 %in% s11){
      exist0[i] <- 1
    } else {
      exist0[i] <- 0
    }
  }
  return(exist0)
} 


for(i in gene_mutation1){
  print(i)
  checkExist(i)
  mutation_summary[,i] <- checkExist(i, mutation_list=mutation_summary)
}

mutation_matrix <- select(mutation_summary, c('strain_name',gene_mutation1 ))
#strain classification based on several genes with a series mutations
#here we can classify the strain based on cluster 
#and the total numer of mutation occured in each strain


#method2 cluster the strain based on the mutation number
##mutation number analysis
mutation_matrix0 <- mutation_matrix
mutation_matrix0$sum <- rowSums(mutation_matrix0[,2:ncol(mutation_matrix0)]) #calculate the number of aimed SNP
table(mutation_matrix0$sum)
plot(density(mutation_matrix0$sum))

# classification based on number
strain1b <- mutation_matrix0[mutation_matrix0$sum >= 2,] %>% select(.,strain_name)
cluster1b <-  growth_phenotype[growth_phenotype$strain_name %in% strain1b$strain_name,]
cluster1b$type <- 'g1b'

strain2b <- mutation_matrix0[mutation_matrix0$sum >= 0 & mutation_matrix0$sum <= 1,] %>% select(.,strain_name)
cluster2b <-  growth_phenotype[growth_phenotype$strain_name %in% strain2b$strain_name,]
cluster2b$type <- 'g2b'




# part 6
# final analysis
# statistical analysis
x1 <- cluster1$growth
x2 <- cluster2$growth
x3 <- cluster1b$growth
x4 <- cluster2b$growth
wilcox.test(x1, x2, alternative = "two.sided")
wilcox.test(x3, x2, alternative = "two.sided")
wilcox.test(x2, x4, alternative = "two.sided")
# plot the result
strain_classification0 <- rbind.data.frame(cluster1,cluster1b,cluster2, cluster2b)
boxplot(growth~type,data=strain_classification0, main="", col=terrain.colors(4),
        xlab="Strain classification", ylab="Relative growth rate",
        ylim=c(0,1), frame=FALSE, cex.lab=1.5, cex.axis=1.2)
