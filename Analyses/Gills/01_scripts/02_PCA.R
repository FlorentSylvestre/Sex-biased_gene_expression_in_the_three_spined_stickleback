##Take a filtered count matrix
##Perform a VST normalisation with deseq2 packages
##Run a PCA and color it according to a groupe file
##Rscript PCA_counts.R count_matrix group_file
##Group file contain two column: samples names (matching those in count_matrix)
#and grouping variable
##########################################################################################################

######################
##Library:############
######################
#Normalisation
library(DESeq2)

#Data formating/plots
library(ggplot2)
library(magrittr)

#####################
##Parse User input:##
#####################
#args = commandArgs(trailingOnly = T)
#count_path = args[1]
#grouping_path = args[2]
#grouping_name = args[3]

count_path = "04_filtered/gene_count_filtered"
grouping_path = "02_metadata/clean_sample_sex"
grouping_name = "Sex"

#####################
##fixed variable#####
#####################
out_path = "05_PCA"

#####################
##Reading input######
#####################
count_matrix <- read.table(count_path, h=T)
grouping <- read.table(grouping_path, h=F)
colnames(grouping) <- c("sample","group")

####################
##ordering input####
####################
count_matrix_ord <- count_matrix[,match(grouping$sample,colnames(count_matrix))]
count_matrix_ord_norm <- vst(as.matrix(count_matrix_ord))

####################
##PCA###############
####################
pca <- prcomp(t(count_matrix_ord_norm), scale = F)

p <-  pca$x %>% data.frame(group = grouping$group, .) %>%
  ggplot(., aes(x= PC1, y = PC2, shape = group, group = group)) +
        geom_point() +
        stat_ellipse()+
        labs(shape = grouping_name) +
        theme_bw()

ggsave(filename = paste(out_path, "PCA.png", sep = "/"),
       plot = p,
       device = "png",
       dpi = 300)
ggsave(filename = paste(out_path, "PCA.pdf", sep = "/"),
       plot = p,
       device = "pdf",
       dpi = 300)
##Last one is to save a plot object
##I will later use it to make a combined plot
##Using the patchwork package

saveRDS(file = paste(out_path, "PCA.RDS", sep = "/"),
       object = p)
