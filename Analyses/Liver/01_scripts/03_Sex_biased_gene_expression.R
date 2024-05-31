##Take a filtered count matrix
##Perform a median of ratio normalisation following DESEq2 practice
##estimate log fold change and associated wilcoxon ranksumtest for significativity
##Rscript diff_expr_analysis.R count_matrix sex_file qval_thresh, LFC_thresh
##Group file contain two column: samples names (matching those in count_matrix)
#and sex, tab separated
##########################################################################################################

######################
##Library:############
######################
#Data formating/plots
library(ggplot2)
library(magrittr)
library(tidyverse)
#####################
##Function###########
#####################
Normalize<- function(count_matrix){
        log_count <- log2(count_matrix)

        pseudo_sample <- rowMeans(log_count)
        pseudo_sample_noinf <- pseudo_sample[!is.infinite(pseudo_sample)]
        logcount_noinf <- log_count[!is.infinite(pseudo_sample),]

        ratio_to_ref <- logcount_noinf - pseudo_sample_noinf
        scaling_factor <- apply(ratio_to_ref,2, median)
        print(length(scaling_factor))
        Normalized_matrix <- t(t(count_matrix)/(2^scaling_factor))
        return(list("Counts" = Normalized_matrix, "norm_factor" = scaling_factor))
}

#####################
##Parse User input:##
#####################
#args = commandArgs(trailingOnly = T)
#count_path = args[1]
#grouping_path = args[2]
#qval_thresh = as.numeric(args[3])
#LFC_thresh = as.numeric(args[4])

count_path = "04_filtered/gene_count_filtered"
grouping_path = "02_metadata/clean_sample_sex"
qval_thresh = 0.05
LFC_thresh = 1

#####################
##fixed variable#####
#####################
out_path = "06_differential_expression"
colors = c("grey", viridisLite::viridis(n = 2, option = 'H', begin = 0.13, end = 0.85))
names(colors) = c("FALSE", "Males", "Females")

#####################
##Reading input######
#####################
count_matrix <- read.table(count_path, h=T)
grouping <- read.table(grouping_path, h=F)
colnames(grouping) <- c("sample","sex")

####################
##ordering input####
####################
count_matrix <- count_matrix[,match(grouping$sample,colnames(count_matrix))]

####################
##Normalization#####
####################
count_matrix <-rbind(count_matrix)
count_norm <- Normalize(count_matrix)
log2count_norm <- log2(count_norm$Counts + 0.5)


####################
##Analyses##########
####################
wilcoxon <- unlist(apply(count_norm$Counts,
                         1,
                         function(x,sex){wilcox.test(x[sex == "M"], x[sex == "F"])$p.value},
                         grouping$sex))

LFC <- apply(log2count_norm[,grouping$sex == "M"], 1, mean) -
        apply(log2count_norm[,grouping$sex == "F"], 1, mean)

report <- data.frame(gene = rownames(count_matrix),
                     Mean = (apply(count_norm$Counts,1, mean)+1),
                     LFC = LFC,
                     p.value = wilcoxon,
                     fdr = p.adjust(wilcoxon, "BH"))



###################
###plot Autosome###
###################
report$signif <- report$fdr <= qval_thresh
report$signif[report$fdr <= qval_thresh & report$LFC >0] <- "Males"
report$signif[report$fdr <= qval_thresh & report$LFC <0] <- "Females"

auto_annot <- rbind(head(report[order(report$LFC),], n = 10),
               head(report[order(report$LFC, decreasing =T), ], n = 10),
               head(report[order(report$fdr), ], n = 20)) %>%
        tibble() %>%
        distinct()


auto_volcano <- report %>%
        ggplot(aes(x = LFC, y= -log10(fdr), color = signif)) +
        geom_point(size = 1.2) +
        scale_color_manual(values= colors ,breaks  = c('Males', 'Females'))+
        labs(x = "log2 fold change",
             y = "-log10(fdr)",
             main = "",
             color = "") +
        geom_vline(xintercept = c(-1,1)* LFC_thresh, linetype = "longdash") +
        theme_bw()


auto_MA_plot <- report %>%

        ggplot(aes(x = (Mean), y = LFC, col = signif))+
        geom_point(size = 1.2) +
        scale_color_manual(values = colors) +
        labs(x = "Mean gene expression",
             y = "Log2 fold change",
             main ="",
             color = "")+
        scale_x_continuous(trans='log10') +
        geom_hline(yintercept = c(-1,1)*LFC_thresh, linetype = "longdash") +
  theme_bw()


##################################
###### Export plots and files#####
##################################

ggsave(paste(out_path,"auto_volcano.pdf",sep = "/"),
       auto_volcano,
       device = "pdf")

ggsave(paste(out_path,"auto_volcano.png",sep = "/"),
       auto_volcano,
       device = "png")

ggsave(paste(out_path,"auto_MA_plot.pdf",sep = "/"),
       auto_MA_plot,
       device = "pdf")

ggsave(paste(out_path,"auto_MA_plot.png",sep = "/"),
       auto_MA_plot,
       device = "png")

saveRDS(auto_MA_plot,
        paste(out_path,"auto_MA_plot.RDS",sep = "/"))

saveRDS(auto_volcano,
        paste(out_path,"auto_volcano.RDS",sep = "/"))

write.table(report,
            paste(out_path,"report.txt",sep = "/"),
            row.names = F,
            col.names = T,
            quote = F)

write.table(count_norm$norm_factor,
            paste(out_path,"norm_factor", sep ="/"),
            row.names = T,
            col.names = F,
            quote = F)

write.table(auto_annot,
            paste(out_path, "autosomal_main_gene", sep = "/"),
            row.names = F,
            col.names = T,
            quote = F)
