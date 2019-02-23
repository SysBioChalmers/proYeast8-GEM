source('Preprocess data from swiss database.R')
source('some_function_for_quality_analysis.R')


# input the gene information in yeastGEM
#initial version
gene_all <- read_excel("data/gene_list_yeastGEM.xlsx",  sheet = "Sheet1")
gene_all$geneNames <- str_trim(gene_all$geneNames, side = "both")

#check the newly added gene in latest yeastGEM version
proYeast_DataFrame_from_yeastGEM_november <- read_excel("data/proYeast_DataFrame_from_yeastGEM_november.xlsx")
new_gene <- setdiff(proYeast_DataFrame_from_yeastGEM_november$gene, gene_all$geneNames)
new_gene <- new_gene[!is.na(new_gene)]
new_gene0 <- data.frame(geneNames=new_gene, stringsAsFactors = FALSE)
new_gene0$gene_source <- 'new_gene_yeastGEM'
new_gene0$note <- NA

#merge the new gene from yeastGEM update into gene list file-'gene_all'
gene_all <- rbind.data.frame(gene_all, new_gene0)





## experimental protein structure and homology modelling structure
## firsr run the script:  
## Preprocess data from swiss database.R
index_experiment <- which(model_EXP$locus %in% gene_all$geneNames ==TRUE)
pdb_EX <- model_EXP[index_experiment,]
gene_EX <- unique(pdb_EX$locus)

gene_no_pdb <- setdiff(gene_all$geneNames, gene_EX)
index_homology <- which(model_homo$locus %in% gene_no_pdb ==TRUE)
pdb_HOMO <- model_homo[index_homology,]
gene_HOMO <- unique(pdb_HOMO$locus)



## evalution of all homology(HOMO) protein structure
#PDB number analysis
pdb_number <- data.frame(locus = gene_HOMO, stringsAsFactors = FALSE)
pdb_number$pdb <- getMultipleReactionFormula(pdb_HOMO$template,pdb_HOMO$locus,pdb_number$locus)
pdb_number$number <- str_count(pdb_number$pdb, ";") +1
pdb_number <- filter(pdb_number, !is.na(number))


##bar based on group
plotPDBnumber(pdb_number$number)
ggplot(pdb_number, aes(number)) +
  geom_density(fill="lightblue") +
  xlim(-1, 10) +
  labs(x='PDB number per protein')

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
    ##ggplot
ggplot(pdb_HOMO, aes(qmean)) +
  geom_density(fill="blue",alpha=.2) +
  labs(x='QMEAN') +
  labs(y='Density') +
  theme(axis.text=element_text(size=20,face="bold", family="Arial"),
        axis.title=element_text(size=24,face="bold", family="Arial") ) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(legend.position="none") +
  geom_vline(xintercept = -6.98,col="black",lty=2) +
  geom_vline(xintercept = -4,col="red",lty=2)
ggsave(out <- paste('result/','QMEAN','.png', sep = ""), width=8, height=6, dpi=300)






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
      ##ggplot
ggplot(pdb_HOMO_refine, aes(Resolution)) +
  geom_density(fill="blue",alpha=.2) +
  labs(x='Resolution') +
  labs(y='Density') +
  theme(axis.text=element_text(size=20,face="bold", family="Arial"),
        axis.title=element_text(size=24,face="bold", family="Arial") ) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(legend.position="none") +
  geom_vline(xintercept = 3.8,col="black",lty=2) +
  geom_vline(xintercept = 3.4,col="red",lty=2) +
  xlim(0, 8)
ggsave(out <- paste('result/','Resolution','.png', sep = ""), width=8, height=6, dpi=300)







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
##ggplot
ggplot(pdb_HOMO_refine, aes(Seq_Identity)) +
  geom_density(fill="blue",alpha=.2) +
  labs(x='Seq-identity(%)') +
  labs(y='Density') +
  theme(axis.text=element_text(size=20,face="bold", family="Arial"),
        axis.title=element_text(size=24,face="bold", family="Arial") ) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(legend.position="none") +
  geom_vline(xintercept = 8.82,col="black",lty=2) +
  geom_vline(xintercept = 25,col="red",lty=2) +
  xlim(0, 100)
ggsave(out <- paste('result/','Seq-identity','.png', sep = ""), width=8, height=6, dpi=300)







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
##ggplot
ggplot(pdb_HOMO_refine, aes(Seq_similarity)) +
  geom_density(fill="blue",alpha=.2) +
  labs(y='Density') +
  labs(x='Seq_similarity') +
  theme(axis.text=element_text(size=20,face="bold", family="Arial"),
        axis.title=element_text(size=24,face="bold", family="Arial") ) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(legend.position="none") +
  geom_vline(xintercept = 0.25,col="black",lty=2) +
  geom_vline(xintercept = 0.31,col="red",lty=2) +
  xlim(0, 1)
ggsave(out <- paste('result/','Seq_similarity','.png', sep = ""), width=8, height=6, dpi=300)












## choose homology of high quality
## Refine the Homology PDB files based on QMEAN, Resolution, similarity and identity
#  filter 1 with a higher standard
pdb_homo1 <- filter(pdb_HOMO, 
                    qmean >= -4 & 
                    Seq_Identity >=25 &
                    Seq_similarity >=0.31 &
                    Resolution <= 3.4  )
gene_homo1 <- unique(pdb_homo1$locus)



# filter 2 with a lower standard
pdb_homo2 <- filter(pdb_HOMO, 
                    qmean >= -6.98 & 
                      Seq_Identity >= 17.58 &
                      Seq_similarity >= 0.25 &
                      Resolution <= 3.8  )
gene_homo2 <- unique(pdb_homo2$locus)
length(gene_homo2)

# save the result for the pdb_homo with high quality
# the file can be used as the input to calculate the 
write.table(pdb_homo1,"result/pdb_homo for PDB structure without experiment pdb.txt", row.names = FALSE, sep = "\t")




#check the newly proteins with PDB-homo in high quality
library(readxl)
pdb_homo_summary <- read_excel("result/pdb_homo summary for manual check.xlsx")
setdiff(gene_homo1, pdb_homo_summary$locus)




library(tidyverse)
# systematically evaluate the effect of different cut-off on number of proteins
df <- function(s0,s1){
  s <- data.frame(x=s0,y=s1,stringsAsFactors = FALSE)
  return(s)
}

# check the influence of resolution on number of pdb file
calculatePDBhomo_num <- function(resolution0=3.4, qmean0=-4, SI=25, SS=0.31){
  pdb_homo1 <- filter(pdb_HOMO, 
                      qmean >= qmean0 & 
                        Seq_Identity >=SI &
                        Seq_similarity >=SS &
                        Resolution <= resolution0  )
  gene_homo <- unique(pdb_homo1$locus)
  return(length(gene_homo))
}



# resolution

protein_number <- vector()
resolution_range <- seq(0, 9, by=0.1)
for(i in seq_along(resolution_range)) {
  s1 <- calculatePDBhomo_num(resolution0 = resolution_range[i])
  protein_number[i] <- s1
  
}


df1 <- df(resolution_range,protein_number )
ln <- ggplot(df1, aes(x, y)) + geom_line() +
  geom_vline(xintercept = 3.4, linetype="dashed", color = "red", size=1) +
  labs(x ="resolution_range", y = "protein_number") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


# qmean
protein_number1 <- vector()
qmean_range <- seq(-10, 9, by=0.1)
for(i in seq_along(qmean_range)) {
  s1 <- calculatePDBhomo_num(qmean0 = qmean_range[i])
  protein_number1[i] <- s1
  
}

df2 <- df(qmean_range,protein_number1 )
lm <- ggplot(df2, aes(x, y)) + geom_line() +
  geom_vline(xintercept = -4, linetype="dashed", color = "red", size=1) +
  labs(x ="qmean_range", y = "protein_number") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


# SI
protein_number2 <- vector()
SI_range <- seq(10, 99, by=1)
for(i in seq_along(SI_range)) {
  s1 <- calculatePDBhomo_num(SI = SI_range[i])
  protein_number2[i] <- s1
  
}
df3 <- df(SI_range,protein_number2 )
lq <- ggplot(df3, aes(x, y)) + geom_line() +
  geom_vline(xintercept = 25, linetype="dashed", color = "red", size=1) +
  labs(x ="SI_range", y = "protein_number") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


# SS
protein_number3 <- vector()
SS_range <- seq(0, 1, by=0.01)
for(i in seq_along(SS_range)) {
  s1 <- calculatePDBhomo_num(SS = SS_range[i])
  protein_number3[i] <- s1
  
}
df4 <- df(SS_range,protein_number3 )
lp <- ggplot(df4, aes(x, y)) + geom_line() +
  geom_vline(xintercept = 0.31, linetype="dashed", color = "red", size=1) +
  labs(x ="SS_range", y = "protein_number") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))




library(ggpubr)
figure <- ggarrange(ln, lm, lq,lp,
                    ncol = 2, nrow = 2)
figure


# analysis the qmean for all the pdb_homo from swiss model database
mean_qmean <- mean(as.numeric(model_homo$qmean))
model_homo$qmean <- as.numeric(model_homo$qmean)
sd_qmean <- sd(as.numeric(model_homo$qmean))

ggplot(model_homo, aes(qmean)) +
  geom_density(fill="lightblue") +
  labs(x='QMEAN4') +
  geom_vline(aes(xintercept= -4), color="red", linetype="dashed", size=1) +
  geom_vline(aes(xintercept= -5), color="blue", linetype="dashed", size=1)

model_homo1 <- filter(model_homo, qmean >=-4.5)

library(fitdistrplus)
descdist(model_homo$qmean, discrete =  FALSE)
fit.norm1 <- fitdist(model_homo$qmean,'norm')
plot(fit.norm1)


