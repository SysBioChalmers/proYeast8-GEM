# some general function used for the protein quality analysis
# 12th, November, 2018
# hongzhong Lu

library(ggplot2)
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

# function to compare whether two lists contain the same elements
SameElements <- function(a, b) return(identical(sort(a), sort(b)))
a=c('a',1,'b')
b=c('b',1,'a')
SameElements(a,b)


#function to plot the pdb file number distribution

##bar based on group
plotPDBnumber <- function(pdb_number0) {
  #input pdb_number0 a vector contains the numbers of pdb files for each gene
  group <- c("(1)", "(2,10)", "(10,20)", "(20,50)", "(50,200)")
  pdb_analysis <- data.frame(group = group, stringsAsFactors = FALSE)
  pdb_analysis$group <- factor(pdb_analysis$group, levels = pdb_analysis$group)
  pdb_analysis$num[1] <- number.counter(pdb_number0, 1, 2)
  pdb_analysis$num[2] <- number.counter(pdb_number0, 2, 10)
  pdb_analysis$num[3] <- number.counter(pdb_number0, 10, 20)
  pdb_analysis$num[4] <- number.counter(pdb_number0, 20, 50)
  pdb_analysis$num[5] <- number.counter(pdb_number0, 50, 200)
  ggplot(data = pdb_analysis, aes(x =group , y = num, fill = group)) +
    geom_bar(stat = "identity") + # reorder: adjust the order
    theme(legend.title = element_blank(), legend.position = "right") +
    labs(y = "Number of protein", x = "PDB_EX number") +
    theme(legend.title = element_text(size = 10)) +
    theme(axis.text=element_text(size=20,face="bold", family="Arial"),
          axis.title=element_text(size=24,face="bold", family="Arial") ) +
    ggtitle('') +
    theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
    theme(legend.position="none")
}


# function remove the duplicated number in a column of dataframe
calculateUnique <- function (s1){
  #s1 <- pfam_domain_number$domain
  s2 <- str_split(s1, ";")
  s3 <- sapply(s2, unique)
  s4 <- sapply(s3, length)
  return(s4)
}

