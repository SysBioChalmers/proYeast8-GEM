#evalue the seq leng

#first example with intron
gene_test1 <- "YER003C.fasta"
df_list <- processFasta(gene_test = gene_test1)
seq_len <- vector()

for (i in seq_along(df_list)){
  seq_len[i] <- length(df_list[[i]])
}

#based on SGD this gene with intron has 1383 
right_seq <- length(which(seq_len == 1383))
wrong_seq <- length(df_list) - right_seq





#second example without intron
gene_test1 <- "YAL012W.fasta"
df_list <- processFasta(gene_test = gene_test1)
seq_len <- vector()

for (i in seq_along(df_list)){
  seq_len[i] <- length(df_list[[i]])
}

#based on SGD this gene with intron has 1185
right_seq <- length(which(seq_len == 1185))
wrong_seq <- length(df_list) - right_seq






#third example without intron
gene_test1 <- "YCL018W.fasta"
df_list <- processFasta(gene_test = gene_test1)
seq_len <- vector()

for (i in seq_along(df_list)){
  seq_len[i] <- length(df_list[[i]])
}

#based on SGD this gene with intron has 1095
right_seq <- length(which(seq_len == 1095))
wrong_seq <- length(df_list) - right_seq




#four example without intron
gene_test1 <- "YBR006W.fasta"
df_list <- processFasta(gene_test = gene_test1)
seq_len <- vector()

for (i in seq_along(df_list)){
  seq_len[i] <- length(df_list[[i]])
}

#based on SGD this gene with intron has 1494
right_seq <- length(which(seq_len == 1494))
wrong_seq <- length(df_list) - right_seq


which(seq_len != 1494)



