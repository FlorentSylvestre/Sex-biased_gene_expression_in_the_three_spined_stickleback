#Packages
library(tidyverse)
library(data.table)
library(ggvenn)
library(ComplexHeatmap)

#functions
parse_gene_count_per_chr <- function(repport, tissue){
  table_s = table(repport$Chrom[!startsWith(repport$gene, "trna") & repport$signif != "FALSE"])
  table_ns = table(repport$Chrom[!startsWith(repport$gene, "trna")])
  table_sm = table(repport$Chrom[!startsWith(repport$gene, "trna") & repport$signif == "Males"])
  table_sf = table(repport$Chrom[!startsWith(repport$gene, "trna") & repport$signif == "Females"])
  res = Reduce(function(x, y) merge(x, y, by = "chr", all = T),
               list(S = data.frame(chr = names(table_s), count_s = as.vector(unname(table_s))),
                    NS = data.frame(chr = names(table_ns), count_ns = as.vector(table_ns)),
                    SM = data.frame(chr = names(table_sm), count_sm = as.vector(table_sm)),
                    SF = data.frame(chr = names(table_sf), count_sf = as.vector(table_sf))))
  return(data.frame(Tissue = tissue,
                    chr = res$chr,
                    All_significant = res$count_s,
                    Males = res$count_sm,
                    Females = res$count_sf,
                    Non_Significant = res$count_ns))
}

conv_num_chr <- function(listchromnum){
  convtable <-   c("chrI", "chrII","chrIII","chrIV",
                   "chrV","chrVI", "chrVII", "chrVIII",
                   "chrIX","chrX", "chrXI", "chrXII", "chrXIII",
                   "chrXIV","chrXV","chrXVI","chrXVII", "chrXVIII",
                   "chrXIX", "chrXX","chrXXI")
  return(convtable[listchromnum])
}


#setup color
colors = viridisLite::viridis(n = 2, option = 'H', begin = 0.13, end = 0.85)
names(colors) = c("Males", "Females")


#input files:
Liver_counts <- read.table("../Liver/04_filtered/gene_count_filtered")
Brain_counts <- read.table("../Brain/04_filtered/gene_count_filtered")
Gills_counts <- read.table("../Gills/04_filtered/gene_count_filtered")

Liver_LFC <- fread("../Liver/06_differential_expression/report.txt") %>% as.data.frame()
Brain_LFC <- fread("../Brain/06_differential_expression/report.txt") %>% as.data.frame()
Gills_LFC <- fread("../Gills/06_differential_expression/report.txt") %>% as.data.frame()

Liver_popmap <- fread("../Liver/02_metadata/clean_sample_sex", h= F, col.names = c("Ind", "Sex")) %>% as.data.frame()
Brain_popmap <- fread("../Brain/02_metadata/clean_sample_sex", h= F, col.names = c("Ind", "Sex")) %>% as.data.frame()
Gills_popmap <- fread("../Gills/02_metadata/clean_sample_sex", h = F, col.names = c("Ind", "Sex")) %>% as.data.frame()

gene_pos <- fread("01_data/gene_position", col.names = c("Chrom","Start","End","gene"))

##manipulating input file:
#adding pos to LFC:
Liver_LFC <- merge(Liver_LFC, gene_pos, by = "gene")
Brain_LFC <- merge(Brain_LFC, gene_pos, by = "gene")
Gills_LFC <- merge(Gills_LFC, gene_pos, by = "gene")


##Shared genes datasets
all_shared_genes <- Reduce(intersect, list(Liver_LFC$gene[Liver_LFC$signif != "FALSE" ],
                                           Brain_LFC$gene[Brain_LFC$signif != "FALSE"],
                                           Gills_LFC$gene[Gills_LFC$signif != "FALSE"]))
BL_shared_genes <- Reduce(intersect, list(Liver_LFC$gene[Liver_LFC$signif != "FALSE" ],
                                          Brain_LFC$gene[Brain_LFC$signif != "FALSE"]))
LG_shared_genes <- Reduce(intersect, list(Liver_LFC$gene[Liver_LFC$signif != "FALSE" ],
                                          Gills_LFC$gene[Gills_LFC$signif != "FALSE"]))
BG_shared_genes <- Reduce(intersect, list(Brain_LFC$gene[Brain_LFC$signif != "FALSE"],
                                          Gills_LFC$gene[Gills_LFC$signif != "FALSE"]))

all_shared <- data.frame(rbind(cbind(Liver_LFC[Liver_LFC$gene %in% all_shared_genes,],Tissue = "Liver"),
                               cbind(Brain_LFC[Brain_LFC$gene %in% all_shared_genes,],Tissue = "Brain"),
                               cbind(Gills_LFC[Gills_LFC$gene %in% all_shared_genes,],Tissue = "Gills")))
BL_shared <- data.frame(rbind(cbind(Liver_LFC[Liver_LFC$gene %in% BL_shared_genes,],Tissue = "Liver"),
                              cbind(Brain_LFC[Brain_LFC$gene %in% BL_shared_genes,],Tissue = "Brain")))

LG_shared <- data.frame(rbind(cbind(Liver_LFC[Liver_LFC$gene %in% LG_shared_genes,],Tissue = "Liver"),
                              cbind(Gills_LFC[Gills_LFC$gene %in% LG_shared_genes,],Tissue = "Gills")))
BG_shared <- data.frame(rbind(cbind(Brain_LFC[Brain_LFC$gene %in% BG_shared_genes,],Tissue = "Brain"),
                              cbind(Gills_LFC[Gills_LFC$gene %in% BG_shared_genes,],Tissue = "Gills")))

##plots
##Heatmaps
Figure_3_data %<>%
  select(LFC, Tissue, gene) %>%
  pivot_wider(names_from = Tissue, values_from = LFC) %>%
  as.data.frame()%>%
  mutate(order = rowSums(across(c("Liver","Brain","Gills")) >0))%>%
  arrange(order) %>%
  select(-gene, -order)%>%
  as.matrix() %>%
  t()
all_shared <- all_shared[c("Brain", "Gills", "Liver"),]

png("figure/99_figure/Figure_3.png")
tmp_dev <- dev.cur()

pdf("figure/99_figure/Figure_3.pdf",
    width = 4,
    height = 2)
dev.control("enable")

Heatmap(Figure_3_data,
        cluster_rows = F,
        cluster_columns = F,
        col = circlize::colorRamp2(breaks = sort(c(0,range(all_shared))),colors = c(colors[2], "white", colors[1])),
        top_annotation = HeatmapAnnotation(df = data.frame(shared = c(rep("Females", 24), rep("Ambiguous",2), rep("Males", 16))),
                                           which = "col",col = list(shared=c(Females = "black",Ambiguous= "orange",Males ="grey"))))

dev.copy(which = tmp_dev)
dev.off()
dev.off()

Figure_S1_data %<>%
  select(LFC, Tissue, gene) %>%
  pivot_wider(names_from = Tissue, values_from = LFC) %>%
  as.data.frame()%>%
  mutate(order = rowSums(across(c("Brain","Gills")) >0))%>%
  arrange(order) %>%
  select(-gene, -order)%>%
  as.matrix() %>%
  t()


png("figure/99_figure/Figure_S1.png")
tmp_dev <- dev.cur()
pdf("figure/99_figure/Figure_S1.pdf")
dev.control("enable")

Heatmap(Figure_S1_data,
        cluster_rows = F,
        cluster_columns = F,
        col = circlize::colorRamp2(breaks = sort(c(0,range(BG_shared))),colors = c(colors[2], "white", colors[1])),
        top_annotation = HeatmapAnnotation(df = data.frame(shared = c(rep("Females", 35), rep("Ambiguous",0), rep("Males", 21))),
                                           which = "col",col = list(shared=c(Females = "black",Ambiguous= "orange",Males="grey"))))

dev.copy(which = tmp_dev)
dev.off()
dev.off()

Figure_S2_data %<>%
  select(LFC, Tissue, gene) %>%
  pivot_wider(names_from = Tissue, values_from = LFC) %>%
  as.data.frame()%>%
  mutate(order = rowSums(across(c("Brain","Liver")) >0))%>%
  arrange(order) %>%
  select(-gene, -order)%>%
  as.matrix() %>%
  t()

png("figure/99_figure/Figure_S2.png")
tmp_dev <- dev.cur()
pdf("figure/99_figure/Figure_S2.pdf")
dev.control("enable")

Heatmap(Figure_S2_data,
        cluster_rows = F,
        cluster_columns = F,
        col = circlize::colorRamp2(breaks = sort(c(0,range(BL_shared))),colors = c(colors[2], "white", colors[1])),
        top_annotation = HeatmapAnnotation(df = data.frame(shared = c(rep("Males", 35), rep("Ambiguous",9), rep("F", 31))),
                                           which = "col",col = list(shared=c(Males = "black",Ambiguous= "orange",F="grey"))))

dev.copy(which = tmp_dev)
dev.off()
dev.off()

Figure_S3_data %<>%
  select(LFC, Tissue, gene) %>%
  pivot_wider(names_from = Tissue, values_from = LFC) %>%
  as.data.frame()%>%
  mutate(order = rowSums(across(c("Gills","Liver")) >0))%>%
  arrange(order) %>%
  select(-gene, -order)%>%
  as.matrix() %>%
  t()

png("figure/99_figure/Figure_S3.png")
tmp_dev <- dev.cur()
pdf("figure/99_figure/Figure_S3.pdf")
dev.control("enable")

Heatmap(Figure_S3_data,
        cluster_rows = F,
        cluster_columns = F,
        col = circlize::colorRamp2(breaks = sort(c(0,range(LG_shared))),colors = c(colors[1], "white", colors[2])),
        top_annotation = HeatmapAnnotation(df = data.frame(shared = c(rep("Females", 115), rep("Ambiguous",80), rep("Males", 130))),
                                           which = "col",col = list(shared=c(Females = "black",Ambiguous= "orange",Males="grey"))))
dev.copy(which = tmp_dev)
dev.off()
dev.off()



##Distribution of DEG across chromosomes
Distrib_DEG <- data.frame(rbind(parse_gene_count_per_chr(Brain_LFC, "Brain"),
                                parse_gene_count_per_chr(Liver_LFC, "Liver"),
                                parse_gene_count_per_chr(Gills_LFC, "Gills")))

res <- list()
for(tissue in unique(Distrib_DEG$Tissue)){
  subdat = Distrib_DEG %>% filter(Tissue == tissue)
  subdat[is.na(subdat)] <- 0
  dist_glob <- c(sum(subdat$Males),
                 sum(subdat$Females),
                 sum(subdat$Non_Significant-subdat$All_significant))
  for(chr in unique(subdat$chr)){
    dist_chr = c(sum(subdat$Males[subdat$chr == chr]),
             sum(subdat$Females[subdat$chr == chr]),
             sum(subdat$Non_Significant[subdat$chr == chr]-subdat$All_significant[subdat$chr == chr]))

  Fisher_DEG = fisher.test(rbind(c(sum(dist_chr[1:2]), dist_chr[3]),
                                 c(sum(dist_glob[1:2])-sum(dist_chr[1:2]), dist_glob[3]-dist_chr[3])))
  Fisher_sex= fisher.test(rbind(dist_chr[1:2],
                                (dist_glob[1:2] - dist_chr[1:2])))


      res[[paste(tissue,chr)]] <- c(tissue, chr,unlist(dist_chr),Fisher_DEG$p.value,Fisher_sex$p.value)
  }
}
res <- as.data.frame(do.call(rbind,res))
colnames(res) <- c("Tissue", "Chr", "m","f","ns","p_n","p_sex")
res$q_n <- p.adjust(res$p_n, "BH")
res$q_sex <- p.adjust(res$p_sex, "BH")


Distrib_DEG$chr <- factor(Distrib_DEG$chr, levels = c("chrI", "chrII","chrIII","chrIV",
                                                      "chrV","chrVI", "chrVII", "chrVIII",
                                                      "chrIX","chrX", "chrXI", "chrXII", "chrXIII",
                                                      "chrXIV","chrXV","chrXVI","chrXVII", "chrXVIII",
                                                      "chrXX","chrXXI","M","chrXIX", "chrY"))

Figure_2 <- Distrib_DEG %>%
  select(chr,Tissue, Males, Females, Non_Significant) %>%
  pivot_longer(cols = c(Males, Females),
               names_to = "Sex",
               values_to = "prop",
               names_prefix = "count_") %>%
  ggplot(aes(x = chr, y = prop/Non_Significant, fill = Sex)) +
  geom_col() +
  scale_fill_manual(values = colors) +
  facet_grid(Tissue~., scales = "free_y") +
  theme_bw() +
  ylab("Proportion of significant genes") +
  xlab("Chromosomes") +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45,vjust = 0.5))

ggsave("figure/99_figure/Figure_2.png",
       plot = Figure_3,
       device = "png",
       dpi = 300)


ggsave("figure/99_figure/Figure_2.pdf",
       plot = Figure_3,
       device = "pdf",
       dpi = 300)

write.table("figure/99_figure/Figure_2_p_values.txt",
            x=res, row.names =F,col.names = T, quote = F)



#Plot Overexpressed_genes
overexpressed = data.frame(t(read.table("Foie/04_filtered/overexpressed_genes")))
popmap_liver = read.table("Foie/02_metadata/clean_sample_sex")
norm_factor_liver = read.table("Foie/06_differential_expression/norm_factor")
overexpressed = overexpressed[order(rownames(overexpressed)),]
norm_factor_liver = norm_factor_liver[order(norm_factor_liver$V1),]
popmap_liver = popmap_liver[order(popmap_liver$V1),]
norm_overexpressed = overexpressed/2^norm_factor_liver$V2
norm_overexpressed$sample = rownames(norm_overexpressed)
norm_overexpressed= merge(norm_overexpressed, popmap_liver, by.x = "sample", by.y = "V1", sort = F)
norm_overexpressed = norm_overexpressed %>%
pivot_longer(cols = colnames(norm_overexpressed)[-c(1,7)],
             names_to ="gene",
            values_to = "count") %>%
  rename(Sex = V2)

Figure_4 = norm_overexpressed %>%
  ggplot(.,aes(x= gene, y = log2(count), col = Sex)) +
  geom_point(position = position_dodge(width = 0.75)) +
  geom_boxplot(alpha =0) +
  theme_bw()+
  scale_color_manual(values = unname(colors[c(2,1)]))

ggsave(plot = Figure_4,
      filename = "figure/99_figure/Overexpressed.png",
      device = "png",
      dpi = 300)

ggsave(plot = Figure_4,
       filename = "figure/99_figure/Overexpressed.pdf",
       device = "pdf",
       dpi = 300)
