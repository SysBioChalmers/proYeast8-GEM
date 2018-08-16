source('Preprocess data from swiss database.R')
source('some_function_for_quality_analysis.R')


# input the gene information in yeastGEM
gene_all <- read_excel("data/gene_list_yeastGEM.xlsx",  sheet = "Sheet1")
gene_all$geneNames <- str_trim(gene_all$geneNames, side = "both")


## experimental protein structure and homology modelling structure
## firsr run the script:  
## Preprocess data from swiss database.R
index_experiment <- which(model_EXP$locus %in% gene_all$geneNames ==TRUE)
pdb_EX <- model_EXP[index_experiment,]
gene_EX <- unique(pdb_EX$locus)

gene_no_pdb <- setdiff(gene_all$geneNames, gene_EX)
index_homology <- which(model_homo$locus %in% gene_no_pdb ==TRUE)
pdb_HOMO <-model_homo[index_homology,]
gene_HOMO <- unique(pdb_HOMO$locus)



## evalution of all homology(HOMO) protein structure
#PDB number analysis
pdb_number <- data.frame(locus = gene_HOMO, stringsAsFactors = FALSE)
pdb_number$pdb <- getMultipleReactionFormula(pdb_HOMO$template,pdb_HOMO$locus,pdb_number$locus)
pdb_number$number <- str_count(pdb_number$pdb, ";") +1
pdb_number <- filter(pdb_number, !is.na(number))


##bar based on group
plotPDBnumber(pdb_number$number)


#qmean
mean_qmean <- mean(as.numeric(pdb_HOMO$qmean))
pdb_HOMO$qmean <- as.numeric(pdb_HOMO$qmean)
sd_qmean <- sd(as.numeric(pdb_HOMO$qmean))

dens0 <- density(pdb_HOMO$qmean)
plot(dens0, frame = FALSE, col = "steelblue",
     main = "Density of QMEAN",
     xlim=c(-15, 5),
     ylim=c(0,0.25),
     xlab="qmean",
     ylab="Density",     
     lty=1, 
     lwd=3)
abline(v=-6.98, col="red",lty=2)
abline(v=-4, col="red")
abline(h=c(0.05,0.1,0.15,0.2,0.25), col="lightgray", lty=2)


pnorm(-1.5,mean_qmean,sd_qmean) #calculate the the probability smaller than -1.5
pnorm(-4,mean_qmean,sd_qmean) #calculate the the probability smaller than -1.5
qnorm(0.2,mean_qmean,sd_qmean)  #qnorm is that you give it a probability, and it returns the number whose cumulative distribution matches the probability
qnorm(0.1,mean_qmean,sd_qmean)

#resolution
pdb_HOMO$Resolution <- str_replace_all(pdb_HOMO$Resolution,'Å','')
pdb_HOMO$Resolution[which(is.na(pdb_HOMO$Resolution)==TRUE)] <- 'NA'
pdb_HOMO_refine<- filter(pdb_HOMO, Resolution !='NA')
pdb_HOMO_refine$Resolution <- as.numeric(pdb_HOMO_refine$Resolution)

max(as.numeric(pdb_HOMO_refine$Resolution))
mean_Resolution <- mean(as.numeric(pdb_HOMO_refine$Resolution))
pdb_HOMO_refine$Resolution <- as.numeric(pdb_HOMO_refine$Resolution)
sd_Resolution <- sd(as.numeric(pdb_HOMO_refine$Resolution))

dens0 <- density(pdb_HOMO_refine$Resolution)
plot(dens0, frame = FALSE, col = "steelblue",
     main = "Density of Resolution",
     xlab="Resolution",
     ylab="Density",
     lty=1,
     lwd=3)
abline(v=3.8, col="red",lty=2)
abline(v=3.4, col="red")
abline(h=c(0.1,0.2,0.4,0.6,0.8), col="lightgray", lty=2)

pnorm(3.4,mean_Resolution,sd_Resolution,lower.tail = FALSE) #calculate the the probability smaller than -1.5
qnorm(0.1,mean_Resolution,sd_Resolution,lower.tail = FALSE)  #qnorm is that you give it a probability, and it returns the number whose cumulative distribution matches the probability
qnorm(0.19,mean_Resolution,sd_Resolution,lower.tail = FALSE)  #qnorm is that you give it a probability, and it returns the number whose cumulative distribution matches the probability


#Seq-identity
pdb_HOMO_refine<- filter(pdb_HOMO, Seq_Identity !='NA')
pdb_HOMO_refine$Seq_Identity <- as.numeric(pdb_HOMO_refine$Seq_Identity)

mean_identity <- mean(as.numeric(pdb_HOMO_refine$Seq_Identity))
sd_identity <- sd(as.numeric(pdb_HOMO_refine$Seq_Identity))

dens0 <- density(pdb_HOMO_refine$Seq_Identity)
plot(dens0, frame = FALSE, col = "steelblue",
     main = "Density of Seq-identity",
     xlab="Seq-identity(%)",
     ylab="Density",
     lty=1, 
     lwd=3)
abline(v=8.82, col="red",lty=2) # P_value = 0.1
abline(v=25, col="red") # P_value = 0.3586
abline(h=c(0.01,0.02,0.03,0.04), col="lightgray", lty=2)

pnorm(25,mean_identity,sd_identity) #calculate the the probability smaller than 25
qnorm(0.1,mean_identity,sd_identity)  #qnorm is that you give it a probability, and it returns the number whose cumulative distribution matches the probability
qnorm(0.3586,mean_identity,sd_identity)  #qnorm is that you give it a probability, and it returns the number whose cumulative distribution matches the probability



#seq-similarity
pdb_HOMO_refine<- filter(pdb_HOMO, Seq_similarity !='NA')
pdb_HOMO_refine$Seq_similarity <- as.numeric(pdb_HOMO_refine$Seq_similarity)

mean_similarity <- mean(as.numeric(pdb_HOMO_refine$Seq_similarity))
sd_similarity <- sd(as.numeric(pdb_HOMO_refine$Seq_similarity))

dens0 <- density(pdb_HOMO_refine$Seq_similarity)
plot(dens0, frame = FALSE, col = "steelblue",
     main = "Density of Seq_similarity",
     xlab="Seq_similarity",
     ylab="Density",
     lty=1, 
     lwd=3)
abline(v=0.25, col="red",lty=2) # P_value = 0.1
abline(v=0.31, col="red") # P_value = 0.2985
abline(h=c(2,4,6,8), col="lightgray", lty=2)


pnorm(0.31,mean_similarity,sd_similarity) #calculate the the probability smaller than -1.5
qnorm(0.1,mean_similarity,sd_similarity)  #qnorm is that you give it a probability, and it returns the number whose cumulative distribution matches the probability
qnorm(0.2985,mean_similarity,sd_similarity)  #qnorm is that you give it a probability, and it returns the number whose cumulative distribution matches the probability


## choose homology of high quality
## Refine the Homology PDB files based on QMEAN, Resolution, similarity and identity
#  filter 1 with a higher standard
pdb_homo1 <- filter(pdb_HOMO, 
                    qmean >= -4 & 
                    Seq_Identity >=25 &
                    Seq_similarity >=0.31 &
                    Resolution <= 3.4  )
gene_homo1 <- unique(pdb_homo1$locus)


#filter 2 with a lower standard
pdb_homo2 <- filter(pdb_HOMO, 
                    qmean >= -6.98 & 
                      Seq_Identity >= 17.58 &
                      Seq_similarity >= 0.25 &
                      Resolution <= 3.8  )
gene_homo2 <- unique(pdb_homo2$locus)









