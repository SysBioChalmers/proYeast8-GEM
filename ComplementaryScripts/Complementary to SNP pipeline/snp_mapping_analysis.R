library(readxl)
library(tidyverse)
library(fitdistrplus)
library(logspline)
library(hongR)

#------part 1 snp number analysis based on strain
# input the snp number for each strain
# this data is coming from 1011 project article
SNP_num_in_each_strain <- read_excel("data/SNP_num_in_each_strain.xlsx")
df1 <- SNP_num_in_each_strain$homozygous_SNPs
summary(unlist(df1))
dens0 <- density(SNP_num_in_each_strain$homozygous_SNPs)

#save the graph in pdf
pdf("result/SNP_distribution.pdf")
plot(dens0, frame = FALSE, col = "steelblue",
     main = "Density of SNP num",
     xlim = c(0,100000),
     xlab="SNP number",
     ylab="Density",     
     lty=1, 
     lwd=3)
dev.off()

#save the result in tiff
tiff("result/SNP_distribution.tiff", width = 6, height = 4, units = 'in', res = 300, compression = 'rle')
plot(dens0, frame = FALSE, col = "steelblue",
     main = "Density of SNP num",
     xlim = c(0,100000),
     xlab="SNP number",
     ylab="Density",     
     lty=1, 
     lwd=3)
dev.off()

tiff("result/SNP_distribution2.tiff", width = 4, height = 6, units = 'in', res = 300, compression = 'rle')
boxplot(df1, ylab="SNP number")
dev.off()


#------part2 snp number analysis based on the metabolic gene
# input SNP for each gene for all strains
# this dataset is produced based on python codes
# so why we exclude the gene without SNP?
# There are two reasons for gene without SNP, first is that the gene is not analysed as it not in reference
# genome , second is that the gene is really don't have snp data across all these strains, but this condition 
# is quite rare
# thus based on the above two reasons, in the enrichment analysis of genes with least and largest SNP number
# the gene with no SNP is not considered.
geneGEM <- read_excel("data/geneGEM with SNP number.xlsx", sheet = "Remove gene with no SNP")
df2 <- geneGEM$SNP_NUM
summary(unlist(df2))

dens1 <- density(geneGEM$SNP_NUM)

tiff("result/SNP_distribution based on gene.tiff", width = 6, height = 4, units = 'in', res = 300, compression = 'rle')
plot(dens1, frame = FALSE, col = "steelblue",
     main = "Density of SNP num",
     xlim = c(0,45000),
     xlab="SNP number",
     ylab="Density",     
     lty=1, 
     lwd=3)
dev.off()

tiff("result/SNP_distribution based on gene2.tiff", width = 4, height = 6, units = 'in', res = 300, compression = 'rle')
boxplot(df2, ylab="SNP number")
dev.off()

#get the distribution for the above dataset
x <- geneGEM$SNP_NUM
descdist(x, discrete =  FALSE)
#fit a weibull distribution and a normal distribution
fit.weibull <- fitdist(x,'weibull')
fit.norm <- fitdist(x,'norm')
plot(fit.norm)
plot(fit.weibull)
#parameter for the weibull distribution
print(fit.weibull) #obtain the parameters
shape0  <- 1.174417
scale0 <- 5137.219642
s = dweibull(x=1:40000, shape=shape0, scale = scale0, log = FALSE) #density (密度函数)
plot(s)
pweibull(q=10000, shape=shape0, scale = scale0, lower.tail = TRUE, log.p = FALSE) #分布函数
qweibull(p=0.01, shape=shape0, scale = scale0, lower.tail = TRUE, log.p = FALSE) #分位函数 start from left side
qweibull(p=0.02,shape=shape0, scale = scale0, lower.tail = FALSE, log.p = FALSE) #分位函数 start from right side
rweibull(n, shape=shape0, scale = scale0)#随机函数


#choose the gene with snp number in two groups
#group 1, pvalue =0.01, left-tailed
gene_with_lessSNP <- filter(geneGEM, SNP_NUM <= 102.238)
gene_with_moreSNP <- filter(geneGEM, SNP_NUM >= 16411.56)

gene_list_lessSNP <- paste0(gene_with_lessSNP$geneNames, collapse = ",")
gene_list_moreSNP <- paste0(gene_with_moreSNP$geneNames, collapse = ",")


#------part3 nsSNP number analysis based on the metabolic gene
# input nsSNP for each gene for all strains
# this dataset is produced based on R codes
nsSNP <- read.table("data/geneGEM_with_SNP_and_nsSNP_number.txt", header = TRUE, sep = "\t")
y <- nsSNP$nsSNP
y <- y[!is.na(y)]
y <- y[y>0]
descdist(y, discrete =  FALSE)

#fit a weibull distribution and a normal distribution
fit.weibull <- fitdist(y,'weibull')
fit.norm <- fitdist(y,'norm')
plot(fit.norm)

plot(fit.weibull)

shape0  <- 0.85744
scale0 <- 1227.97488

s = dweibull(x=1:10000, shape=shape0, scale = scale0, log = FALSE) #density (密度函数)
plot(s)
pweibull(q=10000, shape=shape0, scale = scale0, lower.tail = TRUE, log.p = FALSE) #分布函数
qweibull(p=0.015, shape=shape0, scale = scale0, lower.tail = TRUE, log.p = FALSE) #分位函数 start from left side
qweibull(p=0.025,shape=shape0, scale = scale0, lower.tail = FALSE, log.p = FALSE) #分位函数 start from right side

#choose the gene with snp number in two groups
#group 1, pvalue =0.01, left-tailed
gene_with_less_nsSNP <- filter(nsSNP, nsSNP <= 9.243918)
gene_with_more_nsSNP <- filter(nsSNP, nsSNP >= 5627.781)

#obtain the gene list with less nsSNP and more nsSNP
gene_list_less_nsSNP <- paste0(gene_with_less_nsSNP$geneNames, collapse = ",")
gene_list_more_nsSNP <- paste0(gene_with_more_nsSNP$geneNames, collapse = ",")

# comparison between SNP and nsSNP
# density 
plot(dens1,xlab="",col = "steelblue",
     main = "Comparing between SNP and nsSNP distribution",
     ylim = c(0,0.0006),
     lty=1,
     lwd=3)
lines(dens2, lty=1,lwd=3)
title( xlab="Number of SNP (nsSNP) number")
# vioplot
library(vioplot)
x1 <- geneGEM$SNP_NUM[!is.na(geneGEM$SNP_NUM)]
x2 <- nsSNP$nsSNP[!is.na(nsSNP$nsSNP)]
vioplot(x1, x2, names=c("SNP", "nsSNP"), col="gold")
title("Violin Plots of Miles Per Gallon")
#correlation plot
plot(nsSNP$SNP_NUM,nsSNP$nsSNP, xlab="", ylab="")
lines(lowess(nsSNP$SNP_NUM,nsSNP$nsSNP), col=2)
title( xlab="total SNP number", ylab="nsSNP number")




#------part 4 snp number analysis based on the all gene
#we want to look into the function enrichment of gene with least and largest SNP number in a whole genome
#which can be compared with the above metabolic gene analysis
all_gene_with_SNP_number <- read_excel("data/all_gene_with SNP number.xlsx")
#remove the 113 genes with no SNP
all_gene_with_SNP0 <- filter(all_gene_with_SNP_number, SNP_NUM > 0)
summary(unlist(all_gene_with_SNP0$SNP_NUM))

dens1 <- density(all_gene_with_SNP0$SNP_NUM)
plot(dens1, frame = FALSE, col = "steelblue",
     main = "Density of SNP num",
     xlim = c(0,60000),
     xlab="SNP number",
     ylab="Density",     
     lty=1, 
     lwd=3)

#get the distribution for the above dataset
x <- all_gene_with_SNP0$SNP_NUM
descdist(x, discrete =  FALSE)
#fit a weibull distribution and a normal distribution
fit.weibull <- fitdist(x,'weibull')
fit.norm <- fitdist(x,'norm')
plot(fit.norm)
plot(fit.weibull)
#parameter for the weibull distribution
print(fit.weibull) #obtain the parameters
shape0  <- 1.08
scale0 <- 5182.922
s = dweibull(x=1:40000, shape=shape0, scale = scale0, log = FALSE) #density (密度函数)
plot(s)
pweibull(q=10000, shape=shape0, scale = scale0, lower.tail = TRUE, log.p = FALSE) #分布函数
qweibull(p=0.01, shape=shape0, scale = scale0, lower.tail = TRUE, log.p = FALSE) #分位函数 start from left side
qweibull(p=0.01,shape=shape0, scale = scale0, lower.tail = FALSE, log.p = FALSE) #分位函数 start from right side
rweibull(n, shape=shape0, scale = scale0)#随机函数

#choose the gene with snp number in two groups
#group 1, pvalue =0.01, left-tailed
#choose these genes to do the gene enrichment analysis
lessSNP <- filter(all_gene_with_SNP0, SNP_NUM <= 73.23891)
moreSNP <- filter(all_gene_with_SNP0, SNP_NUM >= 21315.28)
gene_fromGenome_list_lessSNP <- paste0(lessSNP$locus_tag, collapse = ",")
gene_fromGenome_list_moreSNP <- paste0(moreSNP$locus_tag, collapse = ",")





#------part5 nsSNP number analysis based on genome
# input nsSNP for each gene for all strains
# this dataset is produced based on R codes
nsSNP <- read_excel("data/all_gene_with SNP number.xlsx")
nsSNP$nsSNP <- as.numeric(nsSNP$nsSNP )
y <- nsSNP$nsSNP
y <- y[!is.na(y)]
y <- y[y>0]
descdist(y, discrete =  FALSE)

#fit a weibull distribution and a normal distribution
fit.weibull <- fitdist(y,'weibull')
fit.norm <- fitdist(y,'norm')
plot(fit.norm)

plot(fit.weibull)

shape0  <- 0.8320273
scale0 <- 1581.2802268

s = dweibull(x=1:10000, shape=shape0, scale = scale0, log = FALSE) #density (密度函数)
plot(s)
pweibull(q=10000, shape=shape0, scale = scale0, lower.tail = TRUE, log.p = FALSE) #分布函数
qweibull(p=0.01, shape=shape0, scale = scale0, lower.tail = TRUE, log.p = FALSE) #分位函数 start from left side
qweibull(p=0.01,shape=shape0, scale = scale0, lower.tail = FALSE, log.p = FALSE) #分位函数 start from right side

#choose the gene with snp number in two groups
#group 1, pvalue =0.01, left-tailed
gene_with_less_nsSNP <- filter(nsSNP, nsSNP <= 6.278593)
gene_with_more_nsSNP <- filter(nsSNP, nsSNP >= 9911.81)

#obtain the gene list with less nsSNP and more nsSNP
gene_list_less_nsSNP <- paste0(gene_with_less_nsSNP$locus_tag, collapse = ",")
gene_list_more_nsSNP <- paste0(gene_with_more_nsSNP$locus_tag, collapse = ",")

# comparison between SNP and nsSNP
# density 
dens2 <- density(y)
plot(dens1,xlab="",col = "steelblue",
     main = "Comparing between SNP and nsSNP distribution",
     ylim = c(0,0.0006),
     lty=1,
     lwd=3)
lines(dens2, lty=1,lwd=3)
title( xlab="Number of SNP (nsSNP) number")

# vioplot
library(vioplot)
x1 <- x
x2 <- y

vioplot(x1, x2, names=c("SNP", "nsSNP"), col="gold")
title("Comparison between SNP and nsSNP")

# correlation plot
plot(nsSNP$SNP_NUM,nsSNP$nsSNP, xlab="", ylab="")
lines(lowess(nsSNP$SNP_NUM,nsSNP$nsSNP), col=2)
title( xlab="total SNP number", ylab="nsSNP number")




























