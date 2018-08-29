##----------note
##----------the domain analysis will be mainly based on SGD database, not using data from Gang

library(readxl)
library(readr)
library(tidyverse)
library(stringr)
library(hongR) # personal packages contain single and multiple mapping


gene_all <- read_excel("data/gene_list_yeastGEM.xlsx",  sheet = "Sheet1")
gene_all$geneNames <- str_trim(gene_all$geneNames, side = "both")

# yeast gene domain from swiss and SGD
# yeast_domain_swiss0 <- read_excel("data/yeast_domain_Pfam_April.xls", sheet = "Domain")
yeast_domain_swiss <- read_excel("data/20180808_yeast_3D_information.xls", sheet = "Domain")


#id-mapping in swisss database
uniprotID_gene <- read_excel("data/uniprotGeneID_mapping.xlsx")
yeast_domain_swiss$gene <- getSingleReactionFormula(uniprotID_gene$GeneName,uniprotID_gene$Entry,yeast_domain_swiss$`seq id`)

#domain information from pfam fro these metabolic genes
index2 <- which(yeast_domain_swiss$gene %in% gene_all$geneNames ==TRUE)
domain_pfam <- yeast_domain_swiss[index2,]
gene_pfam_withDomain <- unique(domain_pfam$gene)

#simplify the information
domain_pfam0 <- select(domain_pfam,
                       gene,
                       `alignment start`,
                       `alignment end`,
                       `hmm acc`,
                       `hmm name`,
                       type,
                       `E-value`,
                       pdb_id)
colnames(domain_pfam0) <- c('gene','domain_start','domain_end','domain_name','domain_decription','type','e_value','pdb_id')

#domain statistics analysis from pfam database
index3 <- which(duplicated(domain_pfam0$`hmm acc`) ==FALSE)
unique_domain_pfam <- domain_pfam0[index3,]

domain_pfam0$gene_domain <- paste(domain_pfam0$gene, domain_pfam0$`hmm acc`, sep = '_')
index4 <- which(duplicated(domain_pfam0$gene_domain)==FALSE)
unique_gene_domain <- domain_pfam0[index4,]

#domain number of each metabolic genes
unique_gene_pfam <- unique(domain_pfam0$gene)
pfam_domain_number <- data.frame(gene=unique_gene_pfam, stringsAsFactors = FALSE)
pfam_domain_number$domain <- getMultipleReactionFormula(domain_pfam0$`hmm acc`,domain_pfam0$gene,unique_gene_pfam)
pfam_domain_number$number <- str_count(pfam_domain_number$domain, ";") +1

unique(pfam_domain_number$number)
# calculate the number between a range (number1, number2)
number.counter <-function(ss, number1, number2){  
  counts <- 0
  for(i in 1:length(ss)){
    if (number1 <= ss[i] && ss[i] < number2 ){
      counts <- counts + 1
    } else{
      counts <- counts + 0
    }
  }
  return(counts)
}


##bar based on group
group <- c('With one','Between 2 and 6','Between 6 and 20','Between 20 and 100')
pfam_domain_analysis <- data.frame(group=group,stringsAsFactors = FALSE)

pfam_domain_analysis$num[1] <- number.counter(pfam_domain_number$number,1,2)
pfam_domain_analysis$num[2] <- number.counter(pfam_domain_number$number,2,6)
pfam_domain_analysis$num[3] <- number.counter(pfam_domain_number$number,6,20)
pfam_domain_analysis$num[4] <- number.counter(pfam_domain_number$number,20,100)


ggplot(data=pfam_domain_analysis, aes(x=reorder(group,-num), y=num,fill=group)) +
  geom_bar(stat="identity") + # reorder: adjust the order
  theme(legend.title = element_blank(), legend.position = "right") +
  theme(axis.text=element_text(size=10),
         axis.title=element_text(size=12,face="bold")) +
  ggtitle("Number analysis of domain for all genes") +
  labs(x="Group",y="Number") 

## density analysis
##density
d <- density(pfam_domain_number$number)
plot(d, main="Density of domain number",
    xlab="Domain number",
    ylab="Density")

## enrichment analysis for these genes more than 5 domains
gene_over_5_domain <- filter(pfam_domain_number, number <= 200 & number >= 5)
write.table(gene_over_5_domain,'result/gene_over_5_domain.txt', row.names = FALSE, sep='\t')
write.table(domain_pfam0,'result/domain_pfam0_for_SNP_pipeline.txt', row.names = FALSE, sep='\t')





# part 2
yeast_domain_SGD <- read_excel("data/yeast_domain_SGD.xlsx",  sheet = "Sheet4")

# preprocess domain information from SGD database
yeast_domain_SGD <- filter(yeast_domain_SGD, `Protein Domains Method` == 'Pfam' )

# re_name the column names of SGD domain annotation
colnames(yeast_domain_SGD) <- c('SGDid','gene','method','domain_name','domain_decription','e_value','domain_start','domain_end','run_data')

# domain information from sgd for these metabolic genes
index1 <- which(yeast_domain_SGD$gene %in% gene_all$geneNames ==TRUE)
domain_SGD <- yeast_domain_SGD[index1,]
gene_SGD_withDomain <- unique(domain_SGD$gene)

# compare the domain information from SGD and pfam
# it can be found that the latest information from gang is more comprehensive
domain_SGD_YGR240C <- filter(domain_SGD, gene == 'YGR240C')
