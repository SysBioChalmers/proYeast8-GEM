library(readxl)
library(tidyverse)
library(fitdistrplus)
library(logspline)
library(hongR)

#input the genome annotation of gene from SGD database
s288_sgd_annotation <- read_tsv('data/s288_genome.tsv')

# function to standardize the snp number for each gene
calc_relative_snp <- function(geneName, snp_number, annotation = s288_sgd_annotation) {
  #gene_list <- geneGEM$geneNames
  #snp_list <- geneGEM$SNP_NUM
  gene_list <- geneName
  snp_list <- snp_number
  snp_list <- as.numeric(snp_list)
  protein_length <- getSingleReactionFormula(annotation$protein_length, annotation$systematic_name, gene_list)
  protein_length <- as.numeric(protein_length)
  snp_relative <- snp_list / protein_length
  return(snp_relative)
}


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
geneGEM <- read_excel("data/geneGEM with SNP number.xlsx")
#standard the SNP number based on the protein lenght
geneGEM$SNP_NUM <- calc_relative_snp(geneGEM$locus_tag,geneGEM$SNP_NUM)
geneGEM$nsSNP <- calc_relative_snp(geneGEM$locus_tag,geneGEM$nsSNP)
geneGEM0 <- geneGEM[geneGEM$SNP_NUM > 0 & !is.na(geneGEM$SNP_NUM),]

x <- geneGEM0$SNP_NUM
dens1 <- density(x)

#get the distribution for the above dataset
descdist(x, discrete =  FALSE)
#fit a weibull distribution and a normal distribution
fit.weibull <- fitdist(x,'weibull')
fit.norm <- fitdist(x,'norm')
plot(fit.norm)
plot(fit.weibull)
#parameter for the weibull distribution
print(fit.weibull) #obtain the parameters
shape0  <- fit.weibull[['estimate']][1]
scale0 <- fit.weibull[['estimate']][2]
s = dweibull(x=1:60, shape=shape0, scale = scale0, log = FALSE) #density (密度函数)
plot(s)
pweibull(q=50, shape=shape0, scale = scale0, lower.tail = TRUE, log.p = FALSE) #分布函数
n1 <- qweibull(p=0.015, shape=shape0, scale = scale0, lower.tail = TRUE, log.p = FALSE) #分位函数 start from left side
n2 <- qweibull(p=0.02,shape=shape0, scale = scale0, lower.tail = FALSE, log.p = FALSE) #分位函数 start from right side
#rweibull(n, shape=shape0, scale = scale0)#随机函数
#choose the gene with snp number in two groups
#group 1, pvalue =0.01, left-tailed
gene_with_lessSNP <- filter(geneGEM0, SNP_NUM <= n1)
gene_with_moreSNP <- filter(geneGEM0, SNP_NUM >= n2)

gene_list_lessSNP <- paste0(gene_with_lessSNP$locus_tag, collapse = ",")
gene_list_moreSNP <- paste0(gene_with_moreSNP$locus_tag, collapse = ",")


#------part3 nsSNP number analysis based on the metabolic gene
# input nsSNP for each gene for all strains
y <- geneGEM0$nsSNP
y[y==0] <- 0.0001 #in weibull distributoíon, the value should be larger than 0
descdist(y, discrete =  FALSE)

#fit a weibull distribution and a normal distribution
fit.norm <- fitdist(y,'norm')
fit.weibull <- fitdist(y,'weibull')
plot(fit.norm)
plot(fit.weibull)

#parameter for the weibull distribution
print(fit.weibull) #obtain the parameters
shape0  <- fit.weibull[['estimate']][1]
scale0 <- fit.weibull[['estimate']][2]

s = dweibull(x=1:60, shape=shape0, scale = scale0, log = FALSE) #density (密度函数)
plot(s)
pweibull(q=50, shape=shape0, scale = scale0, lower.tail = TRUE, log.p = FALSE) #分布函数
n1 <- qweibull(p=0.025, shape=shape0, scale = scale0, lower.tail = TRUE, log.p = FALSE) #分位函数 start from left side
n2 <- qweibull(p=0.025,shape=shape0, scale = scale0, lower.tail = FALSE, log.p = FALSE) #分位函数 start from right side

#choose the gene with snp number in two groups
#group 1, pvalue =0.01, left-tailed
gene_with_less_nsSNP <- filter( geneGEM0, nsSNP <= n1)
gene_with_more_nsSNP <- filter( geneGEM0, nsSNP >= n2)

#obtain the gene list with less nsSNP and more nsSNP
gene_list_less_nsSNP <- paste0(gene_with_less_nsSNP$locus_tag, collapse = ",")
gene_list_more_nsSNP <- paste0(gene_with_more_nsSNP$locus_tag, collapse = ",")

# comparison between SNP and nsSNP
# density 
dens1 <- density(x)
dens2 <- density(y)
plot(dens1,xlab="",col = "steelblue",
     main = "Comparing between SNP and nsSNP distribution",
     ylim = c(0, 0.3),
     lty=1,
     lwd=3)
lines(dens2, lty=1,lwd=3)
title( xlab="Relative number of SNP (nsSNP)")
# vioplot
library(vioplot)
vioplot(x, y, names=c("SNP", "nsSNP"), col="gold")
title("Comparing between SNP and nsSNP distribution")
# correlation plot
ggplot(geneGEM0, aes(SNP_NUM, nsSNP)) +
  geom_bin2d(bins = 60) + 
  labs(x="total SNP number", y="nsSNP number")



#------part 4 snp number analysis based on the all gene
#we want to look into the function enrichment of gene with least and largest SNP number in a whole genome
#which can be compared with the above metabolic gene analysis
all_gene <- read_excel("data/all_gene_with SNP number.xlsx")

#standardize the snp number based on the protein length
all_gene$SNP_NUM <- calc_relative_snp(all_gene$locus_tag,all_gene$SNP_NUM)
all_gene$nsSNP <- calc_relative_snp(all_gene$locus_tag,all_gene$nsSNP)
all_gene0 <- all_gene[all_gene$SNP_NUM > 0 & !is.na(all_gene$SNP_NUM),]

x <- all_gene0$SNP_NUM
dens1 <- density(x)
plot(dens1, frame = FALSE, col = "steelblue",
     main = "Density of SNP num",
     xlab="SNP number",
     ylab="Density",     
     lty=1, 
     lwd=3)

#get the distribution for the above dataset
descdist(x, discrete =  FALSE)
#fit a weibull distribution and a normal distribution
fit.norm <- fitdist(x,'norm')
fit.weibull <- fitdist(x,'weibull')
plot(fit.norm)
plot(fit.weibull)
#parameter for the weibull distribution
print(fit.weibull) #obtain the parameters
shape0  <- fit.weibull[['estimate']][1]
scale0 <- fit.weibull[['estimate']][2]
s = dweibull(x=1:95, shape=shape0, scale = scale0, log = FALSE) #density (密度函数)
plot(s)
pweibull(q=10000, shape=shape0, scale = scale0, lower.tail = TRUE, log.p = FALSE) #分布函数
n1 <- qweibull(p=0.01, shape=shape0, scale = scale0, lower.tail = TRUE, log.p = FALSE) #分位函数 start from left side
n2 <- qweibull(p=0.01,shape=shape0, scale = scale0, lower.tail = FALSE, log.p = FALSE) #分位函数 start from right side

#choose the gene with snp number in two groups
#group 1, pvalue =0.01, left-tailed
#choose these genes to do the gene enrichment analysis
lessSNP <- filter(all_gene0, SNP_NUM <= n1)
moreSNP <- filter(all_gene0, SNP_NUM >= n2)
gene_fromGenome_list_lessSNP <- paste0(lessSNP$locus_tag, collapse = ",")
gene_fromGenome_list_moreSNP <- paste0(moreSNP$locus_tag, collapse = ",")



#------part5 nsSNP number analysis based on genome
# input nsSNP for each gene for all strains
# this dataset is produced based on R codes
y <- all_gene0$nsSNP
y[y==0] <- 0.0001 #in weibull distributoíon, the value should be larger than 0
descdist(y, discrete =  FALSE)

#fit a weibull distribution and a normal distribution
fit.weibull <- fitdist(y,'weibull')
fit.norm <- fitdist(y,'norm')
plot(fit.norm)
plot(fit.weibull)
print(fit.weibull) #obtain the parameters
shape0  <- fit.weibull[['estimate']][1]
scale0 <- fit.weibull[['estimate']][2]

s = dweibull(x=1:100, shape=shape0, scale = scale0, log = FALSE) #density (密度函数)
plot(s)
pweibull(q=90, shape=shape0, scale = scale0, lower.tail = TRUE, log.p = FALSE) #分布函数
n1 <- qweibull(p=0.01, shape=shape0, scale = scale0, lower.tail = TRUE, log.p = FALSE) #分位函数 start from left side
n2 <- qweibull(p=0.01,shape=shape0, scale = scale0, lower.tail = FALSE, log.p = FALSE) #分位函数 start from right side

#choose the gene with snp number in two groups
#group 1, pvalue =0.01, left-tailed
gene_with_less_nsSNP <- filter(all_gene0, nsSNP <= n1)
gene_with_more_nsSNP <- filter(all_gene0, nsSNP >= n2)

#obtain the gene list with less nsSNP and more nsSNP
gene_list_less_nsSNP <- paste0(gene_with_less_nsSNP$locus_tag, collapse = ",")
gene_list_more_nsSNP <- paste0(gene_with_more_nsSNP$locus_tag, collapse = ",")

# comparison between SNP and nsSNP
# density 
dens2 <- density(y)
plot(dens1,xlab="",col = "steelblue",
     main = "Comparing between SNP and nsSNP distribution",
     ylim = c(0,0.3),
     lty=1,
     lwd=3)
lines(dens2, lty=1,lwd=3)
title( xlab="Relative number of SNP (nsSNP) ")

# vioplot
library(vioplot)
x1 <- x
x2 <- y

vioplot(x1, x2, names=c("SNP", "nsSNP"), col="gold")
title("Comparison between SNP and nsSNP")

# correlation plot
ggplot(all_gene0, aes(SNP_NUM, nsSNP)) +
  geom_bin2d(bins = 80) + 
  labs(x="total SNP number", y="nsSNP number")


















