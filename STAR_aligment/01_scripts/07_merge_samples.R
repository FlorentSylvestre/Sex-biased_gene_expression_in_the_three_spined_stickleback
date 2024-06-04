#!/usr/bin/env Rscript
##Merge counts from HTSEQ
#script <prefix>
#prefix is a common prefix to all individual to be merged
#if empty, merge all individual in 09_gene_counts

#######################
##Variables
#######################

PATH = "09_gene_counts"
PREFIX = "htseq-count_gene_"
OUTPUT = "10_summarised_gene_counts/"
tissue = commandArgs(trailingOnly=TRUE)[1]


#######################
##Script
#######################

path_files_to_merge = list.files(PATH,
                 pattern = paste(PREFIX, tissue,"*", sep =""),
                 full.names =TRUE)

print(path_files_to_merge)
sample_names = substr(path_files_to_merge,
                     33,
                     (nchar(path_files_to_merge)-4))

infos_files = data.frame(path = unlist(path_files_to_merge),
                         names = sample_names)

files_to_merge = apply(infos_files,
                       1,
                       function(x){
                           read.table(x[1], h=F,
                           col.names = c("Gene", x[2]))}
                       )


dataset = Reduce(function(x, y) merge(x, y, all=TRUE, by = "Gene"), files_to_merge)


write.table(x = dataset,
            file = paste(OUTPUT, tissue,"_counts.txt", sep =""),
            col.names = TRUE,
            row.names = FALSE,
            sep = "\t",
            quote = FALSE)
