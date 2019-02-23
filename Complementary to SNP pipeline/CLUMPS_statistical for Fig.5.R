library(tidyverse)
a <- c(35,54,390)
b<- c('Bioethonal','Wild', 'Wine')
strain_choosed <- data.frame(strain = b, number = a, stringsAsFactors = FALSE)
ggplot(strain_choosed, aes(strain, number, fill=strain)) +
  geom_col() +
  xlab('Strain type') +
  ylab('Number of strain') +
  theme(axis.text=element_text(size=20,face="bold", family="Arial"),
        axis.title=element_text(size=24,face="bold", family="Arial"),
        legend.text = element_text(size=20,face="bold", family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(legend.position = '') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(out <- paste('result/','strain number for the clumps analysis','.eps', sep = ""), width=5, height=5, dpi=300)


c <- c(22,65,41)
b <- c('Bioethonal','Wild', 'Wine')
strain_choosed <- data.frame(strain = b, number = c, stringsAsFactors = FALSE)
ggplot(strain_choosed, aes(strain, number, fill=strain)) +
  geom_col() +
  xlab('Strain type') +
  ylab('Number of protein \n with p_value < 0.05') +
  theme(axis.text=element_text(size=20,face="bold", family="Arial"),
        axis.title=element_text(size=24,face="bold", family="Arial"),
        legend.text = element_text(size=20,face="bold", family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(legend.position = '') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(out <- paste('result/','protein number for the clumps analysis','.eps', sep = ""), width=5, height=5, dpi=300)