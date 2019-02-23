library(readxl)
library(tidyverse)
library(hongR)

# main function
parseSingleSite <- function(site_dataframe, site_type0){
  site_type <- site_type0
  colnames(site_dataframe) <- c('Entry',site_type)
  active_site <- site_dataframe   
  active_site[[site_type]] <- str_replace_all(active_site[[site_type]], "\\{.*?\\}", "") %>%
    str_replace_all(.,"\\. \\.","") %>%
    str_replace_all(.,"\\.;",";")
  
  
  active_site1 <- splitAndCombine(active_site[[site_type]], active_site$Entry, sep0 = ";")
  
  ss <- active_site1$v1
  ss1 <- str_split(ss," ")
  singleSite <- vector()
  for (i in seq_along(ss1)) {
    if (site_type %in% ss1[[i]]) {
      ss <- which(ss1[[i]] == site_type)
      singleSite[i] <- ss1[[i]][ss + 1]
    } else {
      singleSite[i] <- "NA"
    }
  }
  
  active_site1$v3 <- singleSite
  active_site1$v4 <- site_type
  
  return(active_site1)
}

parseMutipleSite <- function(site_dataframe, site_type0, single_site = TRUE) {
  site_type <- site_type0
  colnames(site_dataframe) <- c("Entry", site_type)
  active_site <- site_dataframe
  active_site[[site_type]] <- str_replace_all(active_site[[site_type]], "\\{.*?\\}", "") %>%
    str_replace_all(., "\\. \\.", "") %>%
    str_replace_all(., "\\.;", ";")


  active_site1 <- splitAndCombine(active_site[[site_type]], active_site$Entry, sep0 = ";")

  ss <- active_site1$v1
  ss1 <- str_split(ss, " ")
  singleSite <- vector()
  for (i in seq_along(ss1)) {
    if (site_type %in% ss1[[i]]) {
      ss <- which(ss1[[i]] == site_type)
      if (single_site) {
        singleSite[i] <- ss1[[i]][ss + 1]
      } else {
        singleSite[i] <- paste(ss1[[i]][ss + 1], ss1[[i]][ss + 2], sep = "-")
      }
    } else {
      singleSite[i] <- "NA"
    }
  }

  active_site1$v3 <- singleSite
  active_site1$v4 <- site_type

  return(active_site1)
}

# input the data
sce_site00 <- read_excel("data/sce_active site.xlsx")

# parse the active site
active_site0 <- select(sce_site00,Entry,Active_site) %>%
  filter(.,!is.na(Active_site))

active_site1 <- parseSingleSite(active_site0, site_type0 = 'ACT_SITE')

# parse the metal binding
metal_site <- select(sce_site00,Entry, Metal_binding) %>%
  filter(.,!is.na(Metal_binding))

metal_site1 <- parseSingleSite(metal_site, site_type0 = 'METAL')

# parse the binding site 
binding_site <- select(sce_site00,Entry, Binding_site) %>%
  filter(.,!is.na(Binding_site))

binding_site1 <- parseSingleSite(binding_site, site_type0 = 'BINDING')

# parse the modified residue
modified_site <- select(sce_site00,Entry, Modified_residue) %>%
  filter(.,!is.na(Modified_residue))

modified_site1 <- parseSingleSite(modified_site, site_type0 = 'MOD_RES')


# parse the PTM
# will do based on the need
PTM_site <- select(sce_site00,Entry, Post_translational_modification) %>%
  filter(.,!is.na(Post_translational_modification))


# parse the DNA binding
Nucleotide_binding_site <- select(sce_site00,Entry, Nucleotide_binding) %>%
  filter(.,!is.na(Nucleotide_binding))

Nucleotide_binding_site1 <- parseMutipleSite(Nucleotide_binding_site, site_type0 = "NP_BIND", single_site = FALSE)


# merge the above information
sce_site_refine <- rbind.data.frame(active_site1,metal_site1,binding_site1,modified_site1,Nucleotide_binding_site1,Nucleotide_binding_site1)
colnames(sce_site_refine) <- c('description','Entry','coordinate','site_type')
sce_site_refine <- filter(sce_site_refine, coordinate != "NA")

# find the gene orf name
ID_mapping <- read_excel("data/uniprotGeneID_mapping.xlsx")
sce_site_refine$orf <- getMultipleReactionFormula(ID_mapping$GeneName,ID_mapping$Entry,sce_site_refine$Entry)
which(str_detect(sce_site_refine$orf, ";"))

# refine the dataformat
# one entry id from uniprot database could mapping onto gene orf ID
ss <- sce_site_refine %>% unite(.,inf,c("description", "Entry", "coordinate",  "site_type"), sep = "@@")
tt <- str_split(ss$orf,";")
ss0 <- list()
for (i in seq_along(ss$inf)){
  ss0[[i]] <- paste(ss$inf[i],tt[[i]], sep = "@@")
}
ss1 <- unlist(ss0)

sce_site_refine0 <- data.frame(inf = ss1, stringsAsFactors = FALSE)
sce_site_refine0 <- sce_site_refine0 %>% separate(., inf, into = c("description", "Entry", "coordinate",  "site_type","orf"), sep = "@@")
sce_site_refine0$id_p <- paste(sce_site_refine0$orf, sce_site_refine0$coordinate, sep = "@")


write.table(sce_site_refine0, "result/coordinate for the interesting site of sce genome.txt", row.names = FALSE, sep = "\t")












