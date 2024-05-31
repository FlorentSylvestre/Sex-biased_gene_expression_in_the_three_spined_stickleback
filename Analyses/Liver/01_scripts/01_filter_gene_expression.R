############
#Rscript that filter bad individuals based on number of reads in fastq file
#And unexpressed genes based of number of sample with less than X reads
#Usage: Rscript Filtering_HTseq.R ht_count_matrix read_counts min_Nreads cpm_thresh N_samples blacklist_genes
#with
#ht_cout_matrix a txt files with count for each genes (row) and samples (column)
#read_counts a txt files with one row for each samples with sample names and his readsd count in fastq files
#min_Nreads an integer specifying the minimum read numbr required to keep a sample
#cpm_thresh and N_samples value representing genes to keep. Kepping genes with at least cpm_thresh cpm in N_samples
#Also remove any genes in the blacklist_genes list if any

##############################
##Functions###################
##############################
cpm <- function(count_vec){

  return(10^6*count_vec/sum(count_vec))
}

#############################
##Fixed variables############
#############################
Output_path="04_filtered"

#############################
##Parsing user input#########
#############################
args = commandArgs(trailingOnly = T)
count_path <- args[1]
read_counts_path <- args[2]
min_Nreads <- as.numeric(args[3])
cpm_thresh <- as.numeric(args[4])
N_samples <- as.numeric(args[5])
sex_chr_gene_path <- args[6]

##These default corresponds to value used in our paper
if(is.na(min_Nreads)){min_Nreads <- 5*10^6}
if(is.na(cpm_thresh)){cpm_thresh <- 1}
if(is.na(N_samples)){N_samples <- 10}

#############################
##Importing / handling data##
#############################
Raw_counts <- read.table(count_path, h =T)
read_counts <- read.table(read_counts_path, h =F)
if(!is.na(sex_chr_gene_path)){
        blacklist <- read.table(sex_chr_gene_path, h=F)
}

#Removing lines corresponding to unmaped reads:
flags <- c("__alignment_not_unique",
           "__ambiguous",
           "__no_feature",
           "__not_aligned",
           "__too_low_aQual")

Raw_counts <- Raw_counts[!(Raw_counts[,1] %in% flags),]

#Removing genes names from first column
rownames(Raw_counts) <- Raw_counts[,1]
Raw_counts <- Raw_counts[,-1]

#Colnames
colnames(read_counts) <- c("Sample","Coverage")


#############################
##Filtering individuals######
#############################
bad_samples <- read_counts$Sample[read_counts$Coverage < min_Nreads]
Raw_counts_good_samples <- Raw_counts[,!(colnames(Raw_counts) %in% bad_samples)]


#############################
##Filtering genes############
#############################
##Filtering based on low cpm
cpm_counts_good_sample <- data.frame(apply(Raw_counts_good_samples,
                                2,
                                cpm))
expressed_genes <- Raw_counts_good_samples[rowSums(cpm_counts_good_sample>=cpm_thresh) >= N_samples,]
##genes representing more than 15% of total cpm
percent_total_expr <- 100 * as.matrix(expressed_genes) %*% diag(1/colSums(expressed_genes))
overexpressed_genes <- unique(unlist(apply(percent_total_expr, 2 , function(x){ which(x >= 15)})))
overexpressed_genes <- expressed_genes[overexpressed_genes,]

##Removing unwanted genes
if(!is.na(sex_chr_gene_path)){
auto_counts_good_samples <- expressed_genes[!rownames(expressed_genes) %in% blacklist$V1,]
sex_counts_good_samples <- expressed_genes[rownames(expressed_genes) %in% blacklist$V1,]

}



############################
##Outputing files###########
############################

#Expressed gene list
write.table(data.frame(gene = rownames(auto_counts_good_samples)),
            paste(Output_path, "/expressed_genes", sep = ""),
            row.names = F, quote = F, col.names = F)

#Bad samples
write.table(data.frame(bad_samples),
            paste(Output_path,"/bad_samples", sep = ""),
            row.names = F, quote = F, col.names = F)
#Overexpressed_genes
write.table(data.frame(overexpressed_genes),
            paste(Output_path,"/overexpressed_genes", sep = ""),
            row.names = T, quote = F, col.names = T)

#Clean counts
write.table(data.frame(auto_counts_good_samples),
            paste(Output_path,"/gene_count_filtered", sep = ""),
            row.names = T, quote = F, col.names = T)

write.table(data.frame(sex_counts_good_samples),
            paste(Output_path, "/sexchr_gene_counts", sep = ""),
            row.names = T, quote = F, col.names = T)
