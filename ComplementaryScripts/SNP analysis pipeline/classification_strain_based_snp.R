
#----------------note
# this main script is used to classify the strain based on residues mutation information obtained
# from CLUMPS analysis pipeline.
# firstly we obtain the residues mutation of one gene for all strains
# then we can classify the strains based on the specific residue mutations





source("preprocess_1011_project_function.R")
source("genomics annotation summary.R")

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

# function to preprocess the SNP information
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



# datainput
# choose samples that need to be analyzed
strain_classification <- read_excel("data/strain_classification.xls")
strain_type <- "all_strain"
strain_select1 <- chooseStrain(type = strain_type)
# choose the interested protein and the related mutation datasets
gene_mutation  <- read_excel("data/gene_mutation_related_to_heat.xlsx")
gene_mutation0 <- gene_mutation  %>% unite(., mutation_result, gene, alt, position, sep = "@@") %>%
   filter(.,Freq >=3)
gene_mutation1 <- gene_mutation0$mutation_result
ll0 <- length(gene_mutation1)



# produce the gene mutation matrix for all the strains
# YMR246W
p_YMR246W <- printAllMutationWithGene('YMR246W')
# YAR035W 
p_YAR035W <- printAllMutationWithGene('YAR035W')
# YML100W
p_YML100W <- printAllMutationWithGene('YML100W')
# YJL068C
p_YJL068C <- printAllMutationWithGene('YJL068C')
  
#merge the above data
strain_mutation <- rbind.data.frame(p_YMR246W,p_YAR035W,p_YML100W,p_YJL068C)
mutation_summary <- data.frame(strain_name=unique(strain_mutation$strain), stringsAsFactors = FALSE) 
mutation_summary$mutation_result <- getMultipleReactionFormula(strain_mutation$gene_mutation,strain_mutation$strain,mutation_summary$strain_name)
 
#establish a list contains the strain name and mutation 
#here we should produce a matrix to check whether each mutation contains in each strain,

#function to check whether a mutation existed
checkExist <- function(mutate00){
s1 <- str_split(mutation_summary$mutation_result,';')
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
  mutation_summary[,i] <- checkExist(i)
}





mutation_matrix <- select(mutation_summary, c('strain_name',gene_mutation1 ))
#strain classification based on several genes with a series mutations
#here we can classify the strain based on cluster 
#and the total numer of mutation occured in each strain


#method1 cluster the strain based on the maxrix
#cluster analysis
library(ggfortify)
mutation_summary$type <- getSingleReactionFormula(strain_classification$Clades,strain_classification$Standardized_name,mutation_summary$strain_name)
df <- select(mutation_summary, c(gene_mutation1))
name0 <- paste('s',1:ll0,sep = "")
colnames(df) <- name0
autoplot(prcomp(df), data = mutation_summary, colour = 'type')
#obtain the cluster
km <- kmeans(df, 3)
mutation_matrix$cluster <- km$cluster
growth_all <- growth_phenotype$growth

strain1 <- mutation_matrix[mutation_matrix$cluster==1,] %>% select(.,strain_name)
cluster1 <-  growth_phenotype[growth_phenotype$strain_name %in% strain1$strain_name,]

strain2 <- mutation_matrix[mutation_matrix$cluster==2,] %>% select(.,strain_name)
cluster2 <-  growth_phenotype[growth_phenotype$strain_name %in% strain2$strain_name,]

strain3 <- mutation_matrix[mutation_matrix$cluster==3,] %>% select(.,strain_name)
cluster3 <-  growth_phenotype[growth_phenotype$strain_name %in% strain3$strain_name,]
#plot the result
x1 <- cluster1$growth
x2 <- cluster2$growth
x3 <- cluster3$growth
boxplot(x1,x2,x3) 
t.test(x1,x3)
t.test(x1,x2)
t.test(x2,x3)



#method2 cluster the strain based on the mutation number
##based on number
mutation_matrix0 <- mutation_matrix
mutation_matrix0$sum <- rowSums(mutation_matrix0[,2:13])
table(mutation_matrix0$sum)
plot(density(mutation_matrix0$sum))
strain1 <- mutation_matrix0[mutation_matrix0$sum >= 5,] %>% select(.,strain_name)
cluster1 <-  growth_phenotype[growth_phenotype$strain_name %in% strain1$strain_name,]


strain2 <- mutation_matrix0[mutation_matrix0$sum >= 2 & mutation_matrix0$sum <= 4,] %>% select(.,strain_name)
cluster2 <-  growth_phenotype[growth_phenotype$strain_name %in% strain2$strain_name,]


strain3 <- mutation_matrix0[mutation_matrix0$sum <=1,] %>% select(.,strain_name)
cluster3 <-  growth_phenotype[growth_phenotype$strain_name %in% strain3$strain_name,]

x1 <- cluster1$growth
x2 <- cluster2$growth
x3 <- cluster3$growth

##The plots
boxplot(x1,x2,x3)
t.test(x1, x3)  
t.test(x1, x2) 
t.test(x2, x3) 


# strain classification based on two mutations of  'YML100W'
gene_mutation2 <- gene_mutation1[str_detect(gene_mutation1,'YML100W')]
mutation_matrix02 <- select(mutation_summary, c('strain_name',gene_mutation2 ))
mutation_matrix02$sum <- rowSums(mutation_matrix02[,2:3])

strain1 <- mutation_matrix02[mutation_matrix02$sum >= 2,] %>% select(.,strain_name)
cluster1 <-  growth_phenotype[growth_phenotype$strain_name %in% strain1$strain_name,]

strain2 <- mutation_matrix02[mutation_matrix02$sum ==1,] %>% select(.,strain_name)
cluster2 <-  growth_phenotype[growth_phenotype$strain_name %in% strain2$strain_name,]

strain3 <- mutation_matrix02[mutation_matrix02$sum <=0,] %>% select(.,strain_name)
cluster3 <-  growth_phenotype[growth_phenotype$strain_name %in% strain3$strain_name,]

x1 <- cluster1$growth
x2 <- cluster2$growth
summary(x2)
x3 <- cluster3$growth


##The plots
boxplot(x1,x2, x3)
t.test(x1,x3) 
t.test(x2,x3) 
t.test(x1,x2) 




  
  