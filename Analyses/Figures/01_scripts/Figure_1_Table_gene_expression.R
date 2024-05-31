library(tidyverse)
library(patchwork)
library(VennDiagram)
library(ggvenn)


#Constants
colors = viridisLite::viridis(n = 2, option = 'H', begin = 0.13, end = 0.85)
names(colors) = c("Males", "Females")

#Importing subplots:
PCA <- readRDS("../all_tissue/05_PCA/PCA.RDS")
Liver_auto_volc <- readRDS("../Foie/06_differential_expression/auto_volcano.RDS")
Brain_auto_volc <- readRDS("../Cerveau/06_differential_expression/auto_volcano.RDS")
Gills_auto_volc <- readRDS("../Branchie/06_differential_expression/auto_volcano.RDS")

#Importing differential expression reports:
Liver_rep <- read.table("../Foie/06_differential_expression/report.txt", h = T)
Brain_rep <- read.table("../Cerveau/06_differential_expression/report.txt", h = T)
Gills_rep <- read.table("../Branchie/06_differential_expression/report.txt", h = T)

#Venn Diagram
signif_Liver <- Liver_rep %>% filter(signif != FALSE) %>% pull(gene)
signif_Brain <- Brain_rep %>% filter(signif != FALSE) %>% pull(gene)
signif_Gills <- Gills_rep %>% filter(signif != FALSE) %>% pull(gene)

all_tissue <- signif_Liver[signif_Liver %in% signif_Brain & signif_Liver %in% signif_Gills]


venn_signif <- ggvenn(list(Liver = signif_Liver,
                    Brain = signif_Brain,
                    Gills = signif_Gills),
               fill_alpha = 0,text_size = 3,set_name_size = 3,show_percentage = F
               )
##Stacked plot for sex-bias count
Liver_stackp <- Liver_rep %>%
  filter(signif != FALSE) %>%
  pull(signif) %>%
  table() %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(group = row.names(.), percent = V1/sum(V1)) %>%
  ggplot(.,aes(y = percent, x ="", fill = group)) +
  scale_fill_manual(values = colors)+
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.01f%%",percent*100)),
            position = position_stack(vjust = 0.5),
            size = 1.4) +
  theme(line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "none")

brain_stackp <- Brain_rep %>%
  filter(signif != FALSE) %>%
  pull(signif) %>%
  table() %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(group = row.names(.), percent = V1/sum(V1)) %>%
  ggplot(.,aes(y = percent, x ="", fill = group)) +
  scale_fill_manual(values = colors) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.01f%%",percent*100)),
            position = position_stack(vjust = 0.5),
            size = 1.4) +
  theme(line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "none")

gills_stackp <- Gills_rep %>%
  filter(signif != FALSE) %>%
  pull(signif) %>%
  table() %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(group = row.names(.), percent = V1/sum(V1)) %>%
  ggplot(.,aes(x = percent, y ="", fill = group)) +
  scale_fill_manual(values = colors)+
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.01f%%",percent*100)),
            position = position_stack(vjust = 0.5),
            size = 1.4) +
  theme(line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "none")


w = 0.22
L = 0.35
venn_signif <- venn_signif +
  inset_element(brain_stackp, 0.74, 0.47, 0.74+w, 0.47+L) +
  inset_element(Liver_stackp, 0.02, 0.47, 0.02+w,  0.47+L) +
  inset_element(gills_stackp, 0.31, 0.08, 0.31+L , 0.08+w)



##Groupped plot
A <- PCA+theme(legend.position = "none")|venn_signif
B <- (Liver_auto_volc + ggtitle(label = "Liver") | Gills_auto_volc + ylab("")+ ggtitle(label = "Gills") | Brain_auto_volc+ ggtitle(label = "Brain") + ylab("")) +
        plot_layout(guides = 'collect')

Fig1 <- (A/B)

ggsave(file = "figure/99_figure/Fig1.png",
       plot = Fig1,
       device = "png",
       width = 17,
       height = 14,
      units = "cm",
       dpi = 700)

ggsave(file = "figure/99_figure/Fig1.pdf",
       plot = Fig1,
       device = "pdf",
       width = 17,
       height = 14,
       units = "cm",
       dpi = 700)

write.table(all_tissue,
            "figure/98_Gene_list//Table_gene_signif_all.txt",
            row.names = F,
            col.names = F,
            quote = F)
