
#----------------note
# this main script is used to classify the strain based on residues mutation information obtained
# from CLUMPS analysis pipeline.
# firstly we obtain the residues mutation of one gene for all strains
# then we can classify the strains based on the specific residue mutations

source("preprocess_1011_project_function.R")
source("genomics annotation summary.R")


# function to preprocess the SNP information
# this function will obtain all the mutated residue information for one protein
printAllMutationWithGene <- function(ss){
  #example
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
# input the relative growht rate under 42°C
growth_phenotype <- read.table('data/strain_classification based on relative growth under high temperature.txt',header = TRUE, stringsAsFactors = FALSE) 
strain_type <- "all_strain"
strain_select1 <- strain_classification[,'Standardized_name']

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
#-------------------------------------------------
#check the result
#p_YML100W$type <- getSingleReactionFormula(strain_classification$Ecological_origins,strain_classification$Standardized_name,p_YML100W$strain)
#library(tidyverse)
#snp415 <- filter(p_YML100W, gene_mutation =='YML100W@@V@@415'& type=='Bioethanol')
#snp422 <- filter(p_YML100W, gene_mutation =='YML100W@@E@@422'& type=='Bioethanol')
#snp_bioethonal <- filter(p_YML100W, type=='Bioethanol')
#length(unique(snp_bioethonal$strain))
#-------------------------------------------------
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
mutation_matrix0 <- mutation_matrix
mutation_matrix0$sum <- rowSums(mutation_matrix0[,2:13])

# statical analysis for multiple mutations from several genes 
table(mutation_matrix0$sum)
#mutation_matrix0$sum <- as.factor(mutation_matrix0$sum)
#ggplot(mutation_matrix0, aes(sum)) + geom_bar(fill = "#FF6666") +
#  xlab('Mutation number') + ylab('Strain number') +
#  theme_minimal()


# classification all strains in three type
strain1 <- mutation_matrix0[mutation_matrix0$sum >= 6,] %>% select(.,strain_name)
cluster1 <-  growth_phenotype[growth_phenotype$strain_name %in% strain1$strain_name,]
cluster1$type <- 'g1'

strain2 <- mutation_matrix0[mutation_matrix0$sum >= 2 & mutation_matrix0$sum <= 5,] %>% select(.,strain_name)
cluster2 <-  growth_phenotype[growth_phenotype$strain_name %in% strain2$strain_name,]
cluster2$type <- 'g2'

strain3 <- mutation_matrix0[mutation_matrix0$sum <=1,] %>% select(.,strain_name)
cluster3 <-  growth_phenotype[growth_phenotype$strain_name %in% strain3$strain_name,]
cluster3$type <- 'g3'

# statistical analysis
x1 <- cluster1$growth
x2 <- cluster2$growth
x3 <- cluster3$growth

t.test(x1, x3)  
t.test(x1, x2) 
t.test(x2, x3) 

wilcox.test(x1, x3, alternative = "two.sided")
wilcox.test(x1, x2, alternative = "two.sided")
wilcox.test(x2, x3, alternative = "two.sided")

#plot the result
strain_classification0 <- rbind.data.frame(cluster1, cluster2, cluster3)
par(tcl=-0.2)
boxplot(growth~type,data=strain_classification0, main="", col=terrain.colors(4),
        xlab="Strain classification", ylab="Relative growth rate",
        ylim=c(0,1.2), frame=TRUE,cex.lab=1.4, cex.axis=1.2)



# strain classification based on two mutations of  'YML100W'
gene_mutation2 <- gene_mutation1[str_detect(gene_mutation1,'YML100W')]
mutation_matrix02 <- select(mutation_summary, c('strain_name',gene_mutation2 ))
mutation_matrix02$sum <- rowSums(mutation_matrix02[,2:3])

# classification all strains in three type
strain1 <- mutation_matrix02[mutation_matrix02$sum >= 2,] %>% select(.,strain_name)
cluster1 <-  growth_phenotype[growth_phenotype$strain_name %in% strain1$strain_name,]
cluster1$type <- 's1'
strain2 <- mutation_matrix02[mutation_matrix02$sum ==1,] %>% select(.,strain_name)
cluster2 <-  growth_phenotype[growth_phenotype$strain_name %in% strain2$strain_name,]
cluster2$type <- 's2'
strain3 <- mutation_matrix02[mutation_matrix02$sum <=0,] %>% select(.,strain_name)
cluster3 <-  growth_phenotype[growth_phenotype$strain_name %in% strain3$strain_name,]
cluster3$type <- 's3'


# statistical analysis
x1 <- cluster1$growth
x2 <- cluster2$growth
x3 <- cluster3$growth

t.test(x1, x3)  
t.test(x1, x2) 
t.test(x2, x3) 

wilcox.test(x1, x3, alternative = "two.sided")
wilcox.test(x1, x2, alternative = "two.sided")
wilcox.test(x2, x3, alternative = "two.sided")



#plot the result
strain_classification0 <- rbind.data.frame(cluster1, cluster2, cluster3)
par(tcl=-0.2)
boxplot(growth~type,data=strain_classification0, main="", col=terrain.colors(4),
        xlab="Strain classification", ylab="Relative growth rate",
        ylim=c(0,1.2), frame=TRUE, cex.lab=1.4, cex.axis=1.2)


# further check the growth phenotype for the bioethonal strains
growth_phenotype$clade <- getSingleReactionFormula(strain_classification$Clades,strain_classification$Standardized_name,growth_phenotype$strain_name)
# strain name belong to 'g1'(have a higher relative growth in a higher temperature) and bioethonal as well
strain_need_test <- filter(growth_phenotype, str_detect(clade,'bioethanol'))
average_growth <- mean(strain_need_test$growth)

table(mutation_matrix02$sum)


mutation_matrix02$type <- getSingleReactionFormula(strain_classification$Ecological_origins,strain_classification$Standardized_name,mutation_matrix02$strain_name)


