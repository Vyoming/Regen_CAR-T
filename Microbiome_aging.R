#Vyom Shah - Spring 2024
#CAR-T Senescence Project
#Microbiome plotting

#permANOVA test 
library(ggvegan)
micro_matrix <- read_csv("data/Aging_microbiome_beta_diversity.csv", col_names = TRUE) 
micro_matrix <- read_csv("data/Aging_microbiome_beta_diversity_p16.csv", col_names = TRUE) 
micro_matrix <- read_csv("data/Preventative_Revision_Beta_diversity.csv", col_names = TRUE) 
micro_matrix <- read_csv("data/BMT_revision_beta_diversity.csv", col_names = TRUE) 
micro_matrix <- read_csv("data/DQ_revision_beta_diversity.csv", col_names = TRUE) 

micro_matrix1 <- as.dist(micro_matrix)
micro_matrix$sample <- gsub("\\..*", "",colnames(micro_matrix))

micro.div<-adonis2(micro_matrix1~as.factor(micro_matrix$sample), data=micro_matrix, permutations=9999)
micro.div
write.table(micro.div, 'permanova_stat.txt')

set.seed(11) #reproducible results
#write.table(micro.div, 'Micro_permanova_fig4.txt')

dispersion<-betadisper(micro_matrix1, group=micro_matrix$sample)
dispersion

plot(dispersion, hull=FALSE, ellipse=TRUE)

