library(tidyverse)
library(stringr)


## parse the genomic gbff file
s288 <- scan("data/GCF_000146045.2_R64_genomic.gbff", sep = "\n", what = "complex")
s288_n <- s288[72:length(s288)]
s288_n <- str_replace_all(s288_n, "     ", "")


# gene feature summary
gene <- which(str_detect(s288_n, "gene  ") ==TRUE)
gene_name <- gene+1
gene_orf <- gene+2

gene_feature <- data.frame(gene_name=character(length = length(gene_name)), stringsAsFactors = FALSE)
gene_feature$gene_name <- s288_n[gene_name]
gene_feature$gene_orf <- s288_n[gene_orf]
gene_feature$complement <- s288_n[gene]

gene_feature$gene_name <- str_replace_all(gene_feature$gene_name, "/", "")
gene_feature$gene_orf <- str_replace_all(gene_feature$gene_orf, "/", "")

for (i in seq(length(gene_feature$gene_name))){ # seq in R and range in python
  if (str_detect(gene_feature$gene_name[i], "locus_tag")){
    gene_feature$gene_orf[i] <- gene_feature$gene_name[i]
    gene_feature$gene_name[i] <- ""
  } else{
    gene_feature$gene_orf[i] <- gene_feature$gene_orf[i]
    gene_feature$gene_name[i] <- gene_feature$gene_name[i]
    
  }
}

gene_feature0 <- select(gene_feature, gene_orf, complement)
colnames(gene_feature0) <- c('locus_tag','location')

gene_feature0$locus_tag <- str_replace_all(gene_feature0$locus_tag, "locus_tag=","") %>%
  str_replace_all(.,"\"","") %>%
  str_trim(.,side = "both")

gene_feature0$location <- str_replace_all(gene_feature0$location, "gene ","")



#mRNA feature summary
mRNA1 <- which(str_detect(s288_n, "mRNA  ") ==TRUE)
mRNA_feature <- data.frame(location=s288_n[mRNA1], stringsAsFactors = FALSE)
mRNA_feature$gene <- s288_n[mRNA1+1]
mRNA_feature$locus_tag <- s288_n[mRNA1+2]
mRNA_feature$product <- s288_n[mRNA1+3]

for (i in seq(length(mRNA1))){
  if(str_detect(mRNA_feature$gene[i],"locus_tag")){
   mRNA_feature$locus_tag[i] <- mRNA_feature$gene[i]
   mRNA_feature$gene[i] <- ""
  } else{
   mRNA_feature$locus_tag[i] <- mRNA_feature$locus_tag[i]
   mRNA_feature$gene[i] <- mRNA_feature$gene[i]
  }
  
  if(str_detect(mRNA_feature$locus_tag[i],"gene")){
     mRNA_feature$locus_tag[i] <- mRNA_feature$product[i]
     mRNA_feature$location[i] <- paste0(mRNA_feature$location[1], mRNA_feature$gene[1])
  } else{
   mRNA_feature$locus_tag[i] <- mRNA_feature$locus_tag[i]
   mRNA_feature$location[i] <- mRNA_feature$location[i]
  }

}
mRNA_feature0 <- select(mRNA_feature, locus_tag, location)
mRNA_feature0$locus_tag <- str_replace_all(mRNA_feature0$locus_tag, "\\/locus_tag=","") %>%
  str_replace_all(.,"\"","") %>%
  str_trim(.,side = "both")

mRNA_feature0$location <- str_replace_all(mRNA_feature0$location, "mRNA ","")


#CDS feature summary
CDS <- which(str_detect(s288_n, "CDS  ") ==TRUE)
CDS_feature <- data.frame(location=s288_n[CDS], stringsAsFactors = FALSE)
CDS_feature$gene <- s288_n[CDS+1]
CDS_feature$locus_tag <- s288_n[CDS+2]
CDS_feature$note <- s288_n[CDS+3]

for (i in seq(length(CDS))){
  if(str_detect(CDS_feature$gene[i],"locus_tag")){
    CDS_feature$locus_tag[i] <- CDS_feature$gene[i]
    CDS_feature$gene[i] <- ""
  } else{
    CDS_feature$locus_tag[i] <- CDS_feature$locus_tag[i]
    CDS_feature$gene[i] <- CDS_feature$gene[i]
  }
  if(str_detect(CDS_feature$locus_tag[i],"gene")){
    CDS_feature$locus_tag[i] <- CDS_feature$note[i]
    CDS_feature$location[i] <- paste0(CDS_feature$location[1], CDS_feature$gene[1])
  } else{
    CDS_feature$locus_tag[i] <- CDS_feature$locus_tag[i]
    CDS_feature$location[i] <- CDS_feature$location[i]
  }

}

CDS_feature0 <- select(CDS_feature,location,locus_tag)

CDS_feature0$locus_tag <- str_replace_all(CDS_feature0$locus_tag, "\\/locus_tag=","") %>%
  str_replace_all(.,"\"","") %>%
  str_trim(.,side = "both")

CDS_feature0$location <- str_replace_all(CDS_feature0$location, "CDS ","")


## cds fna file analysis
s288_cds <- scan("data/GCF_000146045.2_R64_cds_from_genomic.fna", sep = "\n", what = "complex")
s288_cds[1:200]
index1 <- which(str_detect(s288_cds,">"))

cds <- s288_cds[index1]
cds <- str_replace_all(cds, ">","")
cds <- str_replace_all(cds, " \\[","@")
cds <- str_replace_all(cds, "\\]","")
cds0 <- str_split(cds,"@")

locus <- vector()
gene <- vector()
location_index <- vector()
location <- vector()

for (i in seq(length(cds0))){
  locus[i] <- which(str_detect(cds0[[i]],"locus_tag="))
  gene[i] <- cds0[[i]][locus[i]]
  location_index[i] <- which(str_detect(cds0[[i]],"location="))
  location[i] <- cds0[[i]][location_index[i]]
}

cds_fna <- data.frame(gene=gene, location= location, stringsAsFactors = FALSE)


seq_cds <- list()
for (i in seq(length(index1)-1)){
  seq_cds[[i]] <- s288_cds[(index1[i]+1):(index1[i+1]-1)]
}
seq_cds[[6008]] <- s288_cds[(index1[length(index1)]+1):length(s288_cds)]

nchar(paste(seq_cds[[6008]], sep = "", collapse = ""))


for (i in 1:6008){
  cds_fna$cds[i] <- paste(seq_cds[[i]], sep = "", collapse = "")
  cds_fna$length_cds[i] <- nchar(cds_fna$cds[i])
}

cds_fna$gene <- str_replace_all(cds_fna$gene, "locus_tag=", "")




##gene location summary based on the above three files:
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

##give the mRNA seq, mRNA length, amino acid sequece and amino acid length
library(readr)
s288_SGD <- read_tsv('data/s288_genome.tsv')
gene_feature0$cds_location <- getSingleReactionFormula(cds_fna$location,cds_fna$gene,gene_feature0$locus_tag)
gene_feature0$cds_seq <- getSingleReactionFormula(cds_fna$cds,cds_fna$gene,gene_feature0$locus_tag)
gene_feature0$cds_length <- getSingleReactionFormula(cds_fna$length_cds,cds_fna$gene,gene_feature0$locus_tag)
gene_feature0$aa_seq <- getSingleReactionFormula(s288_SGD$protein_residue,s288_SGD$systematic_name,gene_feature0$locus_tag)
gene_feature0$aa_length <- getSingleReactionFormula(s288_SGD$protein_length,s288_SGD$systematic_name,gene_feature0$locus_tag)
gene_feature0$chromosome <- getSingleReactionFormula(s288_SGD$chromosome,s288_SGD$systematic_name,gene_feature0$locus_tag)
gene_feature0$start <- getSingleReactionFormula(s288_SGD$locations_start,s288_SGD$systematic_name,gene_feature0$locus_tag)
gene_feature0$end <- getSingleReactionFormula(s288_SGD$locations_end,s288_SGD$systematic_name,gene_feature0$locus_tag)
#gene_feature0$DNA_SGD <- getSingleReactionFormula(s288_SGD$sequence,s288_SGD$systematic_name,gene_feature0$locus_tag) It should be note the gene sequence from SGD seems not right

gene_feature0$start <- as.numeric(gene_feature0$start)
gene_feature0$end <- as.numeric(gene_feature0$end)


##choose the metabolic genes
library(readxl)
gene_list_yeastGEM <- read_excel("data/gene_list_yeastGEM.xlsx")
index1 <-  which (gene_feature0$locus_tag %in% gene_list_yeastGEM$geneNames ==TRUE)
gene_feature_GEM <- gene_feature0[index1,]


##evaluate the quality
gene_feature_GEM$check <- ((as.numeric(gene_feature_GEM$cds_length))/3-1) == as.numeric(gene_feature_GEM$aa_length)
