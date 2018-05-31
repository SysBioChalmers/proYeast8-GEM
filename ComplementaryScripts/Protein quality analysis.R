library(readxl)
library(readr)
library(tidyverse)
library(stringr)
getMultipleReactionFormula <- function(description, reaction_ko, ko) {###description can be any charater of metabolite
  index <- vector()
  result <- vector()
  tt <- vector()
  for (i in 1:length(ko)){
    if(length( which (reaction_ko %in%  ko[i]))){
      index <- which (reaction_ko %in%  ko[i])
      tt <- description[index]
      result[i] <- paste0(tt, collapse = ";")
    } else{
      
      result[i] <- NA
    }
  }
  return(result)
}

getSingleReactionFormula <- function(description, reaction_ko, ko) {###description can be any charater of metabolite
  index <- vector()
  result <- vector()
  tt <- vector()
  for (i in 1:length(ko)){
    if(length(match(ko[i],reaction_ko))){
      index <- match(ko[i],reaction_ko)
      tt <- description[index]
      result[i] <- paste0(tt, collapse = ";")
    } else{
      
      result[i] <- NA
    }
  }
  return(result)
}

splitAndCombine <- function(gene,rxn,sep0) { ##one rxn has several genes, this function was used to splite the genes
  
  gene <- str_split(gene, sep0)
  tt<- length(gene)
  gene0 <- list()
  for (i in 1:tt){
    gene0[[i]] <- paste(rxn[i], gene[[i]], sep = "@@@")
    
  }
  
  gene1 <- unique(unlist(gene0))
  gene2 <- str_split(gene1, "@@@" )
  rxnGene <- data.frame(v1=character(length(gene2)),stringsAsFactors = FALSE)
  tt1 <- length(gene2)
  for (j in 1:tt1){
    rxnGene$v1[j] <- gene2[[j]][2]
    rxnGene$v2[j] <- gene2[[j]][1]
  }
  
  return(rxnGene)
}

# calculate the number between a range (number1, number2)
number.counter <-function(ss, number1, number2){  
  counts <- 0
  for(i in 1:length(ss)){
    if (number1 < ss[i] && ss[i] <= number2 ){
      counts <- counts + 1
    } else{
      counts <- counts + 0
    }
  }
  return(counts)
}


# input the gene information in yeastGEM
gene_all <- read_excel("gene_list_yeastGEM.xlsx",  sheet = "Sheet1")
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


### evalution of all experimental protein structure
# PDB number for each genes

pdb_number <- data.frame(locus = gene_EX, stringsAsFactors = FALSE)
pdb_number$pdb <- getMultipleReactionFormula(pdb_EX$template,pdb_EX$locus,pdb_number$locus)
pdb_number$number <- str_count(pdb_number$pdb, ";") +1
pdb_number <- filter(pdb_number, !is.na(number))

hist(pdb_number$number)##hist
d <- density(pdb_number$number) ##density
plot(d, main="Density of pdb number") ##density

      ##bar based on group
group <- c('With one','Between 2 and 10','Between 10 and 20','Between 20 and 50','Between 50 and 200')
pdb_analysis <- data.frame(group=group,stringsAsFactors = FALSE)

pdb_analysis$num[1] <- number.counter(pdb_number$number,1,2)
pdb_analysis$num[2] <- number.counter(pdb_number$number,2,10)
pdb_analysis$num[3] <- number.counter(pdb_number$number,10,20)
pdb_analysis$num[4] <- number.counter(pdb_number$number,20,50)
pdb_analysis$num[5] <- number.counter(pdb_number$number,50,200)
ggplot(data=pdb_analysis, aes(x=reorder(group,-num), y=num,fill=group)) +
  geom_bar(stat="identity") +  # reorder: adjust the order
  theme(legend.title = element_blank(), legend.position = "right")+
  labs(y="Number", x = "PDB_exp number of protein") +
  theme(legend.title=element_text(size=10))


#coverage analysis
min(as.numeric(pdb_EX$coverage))
mean_coverage <- mean(as.numeric(pdb_EX$coverage))
pdb_EX$coverage <- as.numeric(pdb_EX$coverage)
sd_coverage <- sd(as.numeric(pdb_EX$coverage))


dens0 <- density(pdb_EX$coverage)
plot(dens0, frame = FALSE, col = "steelblue",
     main = "Density of coverage",
     xlim=c(0, 1.0),
     ylim=c(0,350),
     xlab="coverage",
     ylab="Density")

plot(ecdf(pdb_EX$coverage),
     main = "Cumulative density of coverage",
     xlim=c(0, 1),
     ylim=c(0,1),
     xlab="Coverage",
     ylab="Cumulative")

#resolution analysis
# preprocess
pdb_EX$Resolution <- str_replace_all(pdb_EX$Resolution,'Å','')
pdb_EX$Resolution[which(is.na(pdb_EX$Resolution)==TRUE)] <- 'NA'

pdb_EX1<- filter(pdb_EX, Resolution !='NA')
pdb_EX1$Resolution <- as.numeric(pdb_EX1$Resolution)

mean_Resolution <- mean(as.numeric(pdb_EX1$Resolution))
pdb_EX1$Resolution <- as.numeric(pdb_EX1$Resolution)
sd_Resolution <- sd(as.numeric(pdb_EX1$Resolution))


dens0 <- density(pdb_EX1$Resolution)
plot(dens0, frame = FALSE, col = "steelblue",
     main = "Density of Resolution",
     xlim=c(0, 10),
     ylim=c(0,0.5),
     xlab="Resolution",
     ylab="Density")
abline(v=3, col="red")


pnorm(0.5,mean_Resolution,sd_Resolution) #calculate the the probability smaller than -1.5
qnorm(0.4,mean_Resolution,sd_Resolution,lower.tail = FALSE)  #qnorm is that you give it a probability, and it returns the number whose cumulative distribution matches the probability


#completeness analysis
pdb_EX <-pdb_EX[which(pdb_EX$mismatch !="NA"),] 
pdb_EX$mismatch <- as.numeric(pdb_EX$mismatch)
dens0 <- density(pdb_EX$mismatch)
plot(dens0, frame = FALSE, col = "steelblue",
     main = "Density of mismatch number",
     xlim=c(0, 20),
     ylim=c(0,0.1),
     xlab="Mismatch number",
     ylab="Density")
abline(v=1, col="red")



## Refine the PDB files based on pident (obtained by blast) and Resolution
#######################################################################
    # filtering using mutation
pdb_EX_filter1 <- filter(pdb_EX,pident == 100)
gene_EX_filter1 <- unique(pdb_EX_filter1$locus)


    # fitering using resolution cut off: 3.4
pdb_EX_filter2 <- filter(pdb_EX_filter1, Resolution <= 3.4 | Resolution =="NA")
gene_EX_filter2 <- unique(pdb_EX_filter2$locus)




pdb_number <- data.frame(locus = gene_EX_filter2, stringsAsFactors = FALSE)
pdb_number$pdb <- getMultipleReactionFormula(pdb_EX_filter2$template,pdb_EX_filter2$locus,pdb_number$locus)
pdb_number$number <- str_count(pdb_number$pdb, ";") +1
pdb_number <- filter(pdb_number, !is.na(number))

hist(pdb_number$number)##hist
d <- density(pdb_number$number) ##density
plot(d, main="Density of pdb number") ##density

##bar based on group
group <- c('With one','Between 2 and 10','Between 10 and 20','Between 20 and 50','Between 50 and 200')
pdb_analysis <- data.frame(group=group,stringsAsFactors = FALSE)

pdb_analysis$num[1] <- number.counter(pdb_number$number,1,2)
pdb_analysis$num[2] <- number.counter(pdb_number$number,2,10)
pdb_analysis$num[3] <- number.counter(pdb_number$number,10,20)
pdb_analysis$num[4] <- number.counter(pdb_number$number,20,50)
pdb_analysis$num[5] <- number.counter(pdb_number$number,50,200)
ggplot(data=pdb_analysis, aes(x=reorder(group,-num), y=num,fill=group)) +
  geom_bar(stat="identity") +  # reorder: adjust the order
  theme(legend.title = element_blank(), legend.position = "right")+
  labs(y="Number", x = "PDB_exp number of protein") +
  theme(legend.title=element_text(size=10))


d <- density(pdb_number$number) ##density
plot(d, main="Density of pdb number") ##density



############################################################################################
## obtain the Homology PDB file for another 49 proteins with mutation in PDB exp
gene_ex_homo <- setdiff(gene_EX, gene_EX_filter1)
index49 <- which(model_homo$locus %in% gene_ex_homo ==TRUE)
pdb_ex_HOMO <- model_homo[index49,]
write.table(pdb_ex_HOMO,"pdb_homo for PDB structure with mutation.txt", row.names = FALSE, sep = "\t")






##########################################################################################
## evalution of all homology(HOMO) protein structure

#PDB number analysis
pdb_number <- data.frame(locus = gene_HOMO, stringsAsFactors = FALSE)
pdb_number$pdb <- getMultipleReactionFormula(pdb_HOMO$template,pdb_HOMO$locus,pdb_number$locus)
pdb_number$number <- str_count(pdb_number$pdb, ";") +1
pdb_number <- filter(pdb_number, !is.na(number))

hist(pdb_number$number)##hist
d <- density(pdb_number$number) ##density
plot(d, main="Density of pdb number") ##density

    ##bar based on group
group <- c('With one','Between 2 and 10','Between 10 and 20','Between 20 and 50','Between 50 and 200')
pdb_analysis <- data.frame(group=group,stringsAsFactors = FALSE)

pdb_analysis$num[1] <- number.counter(pdb_number$number,1,2)
pdb_analysis$num[2] <- number.counter(pdb_number$number,2,10)
pdb_analysis$num[3] <- number.counter(pdb_number$number,10,20)
pdb_analysis$num[4] <- number.counter(pdb_number$number,20,50)
pdb_analysis$num[5] <- number.counter(pdb_number$number,50,200)
ggplot(data=pdb_analysis, aes(x=reorder(group,-num), y=num,fill=group)) +
  geom_bar(stat="identity") +  # reorder: adjust the order
  theme(legend.title = element_blank(), legend.position = "right")+
  labs(y="Number", x = "PDB_homo number of protein") +
  theme(legend.title=element_text(size=10))


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

plot(ecdf(pdb_HOMO$qmean),
     main = "Cumulative density of QMEN",
     xlim=c(-15, 5),
     ylim=c(0,1),
     xlab="qmean",
     ylab="Cumulative")
abline(v=-6.98, col="red")
abline(h=c(0.1,0.2,0.4,0.6,0.8), col="lightgray", lty=2)


pnorm(-1.5,mean_qmean,sd_qmean) #calculate the the probability smaller than -1.5
pnorm(-4,mean_qmean,sd_qmean) #calculate the the probability smaller than -1.5
qnorm(0.2,mean_qmean,sd_qmean)  #qnorm is that you give it a probability, and it returns the number whose cumulative distribution matches the probability
qnorm(0.1,mean_qmean,sd_qmean)
qnorm(0.5,mean_qmean,sd_qmean)

#resolution
  # preprocess
pdb_HOMO$Resolution <- str_replace_all(pdb_HOMO$Resolution,'Å','')
pdb_HOMO$Resolution[which(is.na(pdb_HOMO$Resolution)==TRUE)] <- 'NA'
pdb_HOMO1<- filter(pdb_HOMO, Resolution !='NA')
pdb_HOMO1$Resolution <- as.numeric(pdb_HOMO1$Resolution)

max(as.numeric(pdb_HOMO1$Resolution))
mean_Resolution <- mean(as.numeric(pdb_HOMO1$Resolution))
pdb_HOMO1$Resolution <- as.numeric(pdb_HOMO1$Resolution)
sd_Resolution <- sd(as.numeric(pdb_HOMO1$Resolution))

dens0 <- density(pdb_HOMO1$Resolution)
plot(dens0, frame = FALSE, col = "steelblue",
     main = "Density of Resolution",
     xlim=c(0, 8),
     ylim=c(0,0.8),
     xlab="Resolution",
     ylab="Density",
     lty=1,
     lwd=3)
abline(v=3.8, col="red",lty=2)
abline(v=3.4, col="red")
abline(h=c(0.1,0.2,0.4,0.6,0.8), col="lightgray", lty=2)

plot(ecdf(pdb_HOMO1$Resolution),
     main = "Cumulative density of Resolution",
     xlim=c(0, 8),
     ylim=c(0,1),
     xlab="Resolution",
     ylab="Cumulative")
abline(v=3.8, col="red")
abline(h=c(0.1,0.2,0.4,0.6,0.8), col="lightgray", lty=2)

pnorm(3.4,mean_Resolution,sd_Resolution,lower.tail = FALSE) #calculate the the probability smaller than -1.5
qnorm(0.1,mean_Resolution,sd_Resolution,lower.tail = FALSE)  #qnorm is that you give it a probability, and it returns the number whose cumulative distribution matches the probability
qnorm(0.19,mean_Resolution,sd_Resolution,lower.tail = FALSE)  #qnorm is that you give it a probability, and it returns the number whose cumulative distribution matches the probability

#Seq-identity

# preprocess
pdb_HOMO1<- filter(pdb_HOMO, Seq_Identity !='NA')
pdb_HOMO1$Seq_Identity <- as.numeric(pdb_HOMO1$Seq_Identity)

mean_identity <- mean(as.numeric(pdb_HOMO1$Seq_Identity))
sd_identity <- sd(as.numeric(pdb_HOMO1$Seq_Identity))

dens0 <- density(pdb_HOMO1$Seq_Identity)
plot(dens0, frame = FALSE, col = "steelblue",
     main = "Density of Seq-identity",
     xlim=c(0, 100),
     ylim=c(0,0.04),
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

# preprocess
pdb_HOMO1<- filter(pdb_HOMO, Seq_similarity !='NA')
pdb_HOMO1$Seq_similarity <- as.numeric(pdb_HOMO1$Seq_similarity)

mean_similarity <- mean(as.numeric(pdb_HOMO1$Seq_similarity))
sd_similarity <- sd(as.numeric(pdb_HOMO1$Seq_similarity))

dens0 <- density(pdb_HOMO1$Seq_similarity)
plot(dens0, frame = FALSE, col = "steelblue",
     main = "Density of Seq_similarity",
     xlim=c(0.1, 0.7),
     ylim=c(0,8),
     xlab="Seq_similarity",
     ylab="Density",
     lty=1, 
     lwd=3)
abline(v=0.25, col="red",lty=2) # P_value = 0.1
abline(v=0.31, col="red") # P_value = 0.2985
abline(h=c(2,4,6,8), col="lightgray", lty=2)

plot(ecdf(pdb_HOMO1$Seq_similarity),
     main = "Cumulative density of Seq_similarity",
     xlim=c(0, 8),
     ylim=c(0,1),
     xlab="Seq_similarity",
     ylab="Cumulative")
abline(v=3.8, col="red")
abline(h=c(0.1,0.2,0.4,0.6,0.8), col="lightgray", lty=2)

pnorm(0.31,mean_similarity,sd_similarity) #calculate the the probability smaller than -1.5
qnorm(0.1,mean_similarity,sd_similarity)  #qnorm is that you give it a probability, and it returns the number whose cumulative distribution matches the probability
qnorm(0.2985,mean_similarity,sd_similarity)  #qnorm is that you give it a probability, and it returns the number whose cumulative distribution matches the probability




## Refine the Homology PDB files based on QMEAN, Resolution, similarity and identity
#######################################################################
#number of pdb in higher quality based on higher standard
pdb_HOMO$locus_template <- paste(pdb_HOMO$locus, pdb_HOMO$template, sep = "@")
pdb_homo_h1 <- filter(pdb_HOMO, qmean >= -4) # filter based on qmean
gene_homo_h1 <- unique(pdb_homo_h1$locus_template)

pdb_homo_h2 <- filter(pdb_HOMO,Seq_Identity >=25) # filter based on seq-identity
gene_homo_h2 <- unique(pdb_homo_h2$locus_template)

pdb_homo_h3 <- filter(pdb_HOMO,Seq_similarity >=0.31) # filter based on seq-identity
gene_homo_h3 <- unique(pdb_homo_h3$locus_template)

pdb_homo_h4 <- filter(pdb_HOMO,Resolution <= 3.4) # filter based on seq-identity
gene_homo_h4 <- unique(pdb_homo_h4$locus_template)

gene_homo_h <- intersect(gene_homo_h1,gene_homo_h2)
gene_homo_h <- intersect(gene_homo_h,gene_homo_h3)
gene_homo_h <- intersect(gene_homo_h,gene_homo_h4)

index_high <- which(pdb_HOMO$locus_template %in% gene_homo_h ==TRUE)
pdb_homo_h <- pdb_HOMO[index_high,]
gene_h <- unique(pdb_homo_h$locus)



#number of pdb in higher quality based on lower standard
pdb_HOMO$locus_template <- paste(pdb_HOMO$locus, pdb_HOMO$template, sep = "@")
pdb_homo_h1 <- filter(pdb_HOMO, qmean >= -6.98) # filter based on qmean
gene_homo_h1 <- unique(pdb_homo_h1$locus_template)

pdb_homo_h2 <- filter(pdb_HOMO,Seq_Identity >=17.58) # filter based on seq-identity
gene_homo_h2 <- unique(pdb_homo_h2$locus_template)

pdb_homo_h3 <- filter(pdb_HOMO,Seq_similarity >=0.25) # filter based on seq-identity
gene_homo_h3 <- unique(pdb_homo_h3$locus_template)

pdb_homo_h4 <- filter(pdb_HOMO,Resolution <= 3.8) # filter based on seq-identity
gene_homo_h4 <- unique(pdb_homo_h4$locus_template)

gene_homo_h <- intersect(gene_homo_h1,gene_homo_h2)
gene_homo_h <- intersect(gene_homo_h,gene_homo_h3)
gene_homo_h <- intersect(gene_homo_h,gene_homo_h4)

index_high <- which(pdb_HOMO$locus_template %in% gene_homo_h ==TRUE)
pdb_homo_h <- pdb_HOMO[index_high,]
gene_h <- unique(pdb_homo_h$locus)
























#pie graph
#EXAMPLE: 
#dt = data.frame(A = c(2, 7, 4, 10, 1), B = c('B','A','C','D','E'))

dt = select(pdb_analysis,num,group)
colnames(dt) <- c('A','B')
dt = dt[order(dt$A, decreasing = TRUE),]
myLabel = as.vector(dt$B)   
myLabel = paste(myLabel, "(", round(dt$A / sum(dt$A) * 100, 2), "%)", sep = "")   

ggplot(dt, aes(x = "", y = A, fill = B)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  scale_fill_discrete(breaks = dt$B, labels = myLabel) + 
  theme(axis.text.x = element_blank())







##coverage of each pdb files
group1 <- c('0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1.0')
coverage_analysis <- data.frame(group=group1,stringsAsFactors = FALSE)

coverage_analysis$num[1] <- number.counter(proteinStructure0$coverage,0,0.2)
coverage_analysis$num[2] <- number.counter(proteinStructure0$coverage,0.2,0.4)
coverage_analysis$num[3] <- number.counter(proteinStructure0$coverage,0.4,0.6)
coverage_analysis$num[4] <- number.counter(proteinStructure0$coverage,0.6,0.8)
coverage_analysis$num[5] <- number.counter(proteinStructure0$coverage,0.8,2.0)

#bar graph
ggplot(data=coverage_analysis, aes(x=group1, y=num,fill=group)) +
  geom_bar(stat="identity") + 
  theme(legend.title = element_blank(), legend.position = "right") +
  labs(y="Number", x = "coverage")



#pie graph
#EXAMPLE: 
#dt = data.frame(A = c(2, 7, 4, 10, 1), B = c('B','A','C','D','E'))

dt = select(coverage_analysis,num,group)
colnames(dt) <- c('A','B')
dt = dt[order(dt$A, decreasing = TRUE),]
myLabel = as.vector(dt$B)   
myLabel = paste(myLabel, "(", round(dt$A / sum(dt$A) * 100, 2), "%)", sep = "")   

ggplot(dt, aes(x = "", y = A, fill = B)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "right") + 
  scale_fill_discrete(breaks = dt$B, labels = myLabel) + 
  theme(axis.text.x = element_blank()) +
  geom_text(aes(y = A/2 + c(0, cumsum(A)[-length(A)]), x = sum(A)/5000, label = myLabel), size = 2) 















