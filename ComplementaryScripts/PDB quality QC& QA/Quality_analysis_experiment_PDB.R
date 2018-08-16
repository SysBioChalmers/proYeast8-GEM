source('Preprocess data from swiss database.R')
source('some_function_for_quality_analysis.R')
source('SIFT_yeast.R')

# function used in this part
# calculate the number between a range (number1, number2)
# input the gene information in yeastGEM
gene_all <- read_excel("data/gene_list_yeastGEM.xlsx",  sheet = "Sheet1")
gene_all$geneNames <- str_trim(gene_all$geneNames, side = "both")


index_experiment <- which(model_EXP$locus %in% gene_all$geneNames ==TRUE)
pdb_EX <- model_EXP[index_experiment,]
gene_EX <- unique(pdb_EX$locus)

#The data from swiss database only contains the mapping between the geneID and pdbID, but no chainid
#obtain the chainID for pdb_EX
pdb_EX$chain <- getMultipleReactionFormula(pdb_chain_sce$CHAIN, pdb_chain_sce$id_mapping, pdb_EX$id_mapping)

#we can find for each pair of geneID and pdbID there could be several chainID
#so then we futher establish the relation between geneID, pdbID, chainID
pdb_EX0 <- splitAndCombine(pdb_EX$chain,pdb_EX$id_mapping,sep0 = ";")
colnames(pdb_EX0) <- c('chain_new','id_mapping')

#combine pdb_EX0 and pdb_EX
df_merged <- right_join(pdb_EX, pdb_EX0, by = 'id_mapping')

#then we should bring in the blast results based on sequence using diamond
#it should be noted that there two different types of sequence for each PDB structure
#PDB sequence type1, which could be downloaded from PDB database directly
#PDB sequence type2, which is extracted from each PDB structure
#we should compare these two types of sequence
df_merged$id_mapping_chain <- paste(df_merged$id_mapping, df_merged$chain_new,sep = "@")

#input the blast results using PDB sequence type1
blast1 <- read_excel("data/map_exppdb_prot_id_July.xlsx")
blast1$Entry <- getSingleReactionFormula(geneIDmapping$Entry,geneIDmapping$GeneName,blast1$orf)
blast1$id <- paste(blast1$Entry,str_to_lower(blast1$PDBid),blast1$ChainID, sep = "@")

df_merged$pident1 <- getMultipleReactionFormula(blast1$pident,blast1$id,df_merged$id_mapping_chain)
df_merged$mismatch1 <- getMultipleReactionFormula(blast1$mismatch,blast1$id,df_merged$id_mapping_chain)


#input the blast results using PDB sequence type2
blast2 <- read_excel("data/pdb_seq_blast_sgd.xlsx")
blast2$Entry <- getSingleReactionFormula(geneIDmapping$Entry,geneIDmapping$GeneName,blast2$orf)
blast2 <- blast2 %>% separate(PDBid, into = c('PDBid','ChainID'), sep = "\\.")
blast2$id <- paste(blast2$Entry,str_to_lower(blast2$PDBid),blast2$ChainID, sep = "@")

df_merged$pident2 <- getMultipleReactionFormula(blast2$pident,blast2$id,df_merged$id_mapping_chain)
df_merged$mismatch2 <- getMultipleReactionFormula(blast2$mismatch,blast2$id,df_merged$id_mapping_chain)

#obtain the coordinate of the mapped residue from pdb strucutre
df_merged$qstart2 <- getMultipleReactionFormula(blast2$qstart,blast2$id,df_merged$id_mapping_chain)
df_merged$qend2 <- getMultipleReactionFormula(blast2$qend,blast2$id,df_merged$id_mapping_chain)

#obtain the coordinate of the mapped residue from the original protein sequence
df_merged$sstart2 <- getMultipleReactionFormula(blast2$sstart,blast2$id,df_merged$id_mapping_chain)
df_merged$send2 <- getMultipleReactionFormula(blast2$send,blast2$id,df_merged$id_mapping_chain)


## quality evalution
## evalution of all experimental protein structure
## PDB number for each genes
pdb_number <- data.frame(locus = gene_EX, stringsAsFactors = FALSE)
pdb_number$pdb <- getMultipleReactionFormula(pdb_EX$template,pdb_EX$locus,pdb_number$locus)
pdb_number$number <- str_count(pdb_number$pdb, ";") +1
pdb_number <- filter(pdb_number, !is.na(number))

##bar based on group
plotPDBnumber(pdb_number$number)

#coverage analysis
#the coverage for experimental pdb files is close to 1.
min(as.numeric(pdb_EX$coverage))
mean_coverage <- mean(as.numeric(pdb_EX$coverage))
pdb_EX$coverage <- as.numeric(pdb_EX$coverage)
sd_coverage <- sd(as.numeric(pdb_EX$coverage))


#resolution analysis
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
     xlab="Resolution",
     ylab="Density")
abline(v=3, col="red")
pnorm(0.5,mean_Resolution,sd_Resolution) #calculate the the probability smaller than -1.5
qnorm(0.4,mean_Resolution,sd_Resolution,lower.tail = FALSE)  #qnorm is that you give it a probability, and it returns the number whose cumulative distribution matches the probability


#completeness analysis
pdb_EX2 <- df_merged[which(!is.na(df_merged$mismatch2) ==TRUE),] 
pdb_EX2$mismatch2 <- as.numeric(pdb_EX2$mismatch2)
dens0 <- density(pdb_EX2$mismatch2)
plot(dens0, frame = FALSE, col = "steelblue",
     main = "Density of mismatch number",
     xlab="Mismatch number",
     xlim=c(0,10),
     ylab="Density")
abline(v=1, col="red")


## using the above code, we merge the blast result with each gene-pdbid relation
## next the filtration is used to choose the experimental PDB files without mutation
## for the pdb structure with mutation the homology model will be used
## Refine the PDB files based on pident (obtained by blast) and Resolution
    # filtering using mutation
pdb_EX_filter1 <- filter(df_merged, pident2 == 100)
gene_EX_filter1 <- unique(pdb_EX_filter1$locus)
    # fitering using resolution cut off: 3.4
pdb_EX_filter2 <- filter(pdb_EX_filter1, Resolution <= 3.4 | Resolution =="NA")
gene_EX_filter2 <- unique(pdb_EX_filter2$locus)


pdb_number_refine <- data.frame(locus = gene_EX_filter2, stringsAsFactors = FALSE)
pdb_number_refine$pdb <- getMultipleReactionFormula(pdb_EX_filter2$template,pdb_EX_filter2$locus,pdb_number_refine$locus)
pdb_number_refine$number <- str_count(pdb_number_refine$pdb, ";") +1
pdb_number_refine <- filter(pdb_number_refine, !is.na(number))

##bar based on group
plotPDBnumber(pdb_number_refine$number)


## save the PDB-ex after two round of filtration
write.table(pdb_EX_filter2,"result/pdb_Ex refine for final residue distance calculation.txt", row.names = FALSE, sep = "\t")
write.table(pdb_EX,"result/pdb_EX for PDB structure.txt", row.names = FALSE, sep = "\t")


## obtain the Homology PDB file for another 129 proteins with mutation in PDB exp
gene_ex_homo <- setdiff(gene_EX, gene_EX_filter1)
index2 <- which(model_homo$locus %in% gene_ex_homo ==TRUE)
pdb_ex_HOMO <- model_homo[index2,]
write.table(pdb_ex_HOMO,"result/pdb_homo for PDB structure with mutation.txt", row.names = FALSE, sep = "\t")

## obtain the Homology PDB file for another 28 proteins with smaller resolution
gene_ex_lowR <- setdiff(gene_EX_filter1, gene_EX_filter2)
index3 <- which(model_homo$locus %in% gene_ex_lowR ==TRUE)
pdb_ex_HOMO2 <- model_homo[index3,]
write.table(pdb_ex_HOMO2,"result/pdb_homo for PDB structure with low resolution.txt", row.names = FALSE, sep = "\t")




