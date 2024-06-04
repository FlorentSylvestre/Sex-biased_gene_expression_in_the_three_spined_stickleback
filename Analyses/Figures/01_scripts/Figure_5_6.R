#library
library(tidyverse)
library(patchwork)

#############
##functions##
#############

auto_sex_ratio <- function(dataset){
  return(median(dataset[dataset[,2] == 2, 1]) - median(dataset[dataset[,2] == 1, 1]))
}


##Bootstrap is Adapted from 
##https://stats.stackexchange.com/questions/437477/calculate-accelerated-bootstrap-interval-in-r

bootstrap <- function(data, N, func){
  bootmetric <- rep(NA, N)
  for(i in 1:N){
    samp1 <- sample(which(data[,2]==1), nrow(data[data[,2]==1,]), replace = TRUE)
    samp2 <- sample(which(data[,2]==2), nrow(data[data[,2]==2,]), replace = TRUE)
    bootmetric[i] <- func(data[c(samp1,samp2),])
  }
  return(bootmetric)
}

ACB <- function(data, N.boot, func, alpha = 0.05){
  theta_hat <- func(data)
  theta_boot <- bootstrap(data, N.boot, func)
  a <- 0
  z0 <- qnorm(mean(theta_boot <= theta_hat))
  zu <- qnorm(c(alpha/2, 1- alpha/2))
  u = quantile(theta_boot, pnorm(z0 + (z0+zu)/(1-a*(z0+zu))))
  return(list("IC" = data.frame("Value" = theta_hat, IC_min = u[1], IC_UP = u[2]),
              "bootstrap" = theta_boot,
              "Function" = func,
              "N_sexchr" = sum(data[,2] == 2),
              "N_auto" = sum(data[,2] == 1)))
}

add_strata <- function(gene_pos, brk_X, brk_Y){
  gene_pos$strata <- NA
  for(i in 1:nrow(gene_pos)){
    mid = (gene_pos$end[i] + gene_pos$start[i] )/(2*10^6)  
    if( gene_pos$Chr[i] == "chrXIX"){
      gene_pos$strata[i] = brk_X$Region[which.max(brk_X$pos[mid > brk_X$pos])]
    }
    if( gene_pos$Chr[i] == "chrY"){
      gene_pos$strata[i] = brk_Y$Region[which.max(brk_Y$pos[mid > brk_Y$pos])]
    }
  }
  return(gene_pos)
}

merge_counts <- function(sexcounts, gene_pos, relationship){
  compiled_counts = list()
  female.XX <- list()
  male.X <- list()
  male.Y <- list()
  count = 0
  meta_coll_n <- c(-1,0) + ncol((sexcounts))+2
  
  
  for(gene in 1:nrow(sexcounts)){
    gene_name <- rownames(sexcounts)[gene]
    gene_pos_line <- which(gene_pos$Gene == paste("gene",gene_name, sep ="-"))
    
    #Removing transfer rna from dataset as all copies are merge by star so not exploitable  
    if(substr(gene_name,1,4) == "trna"){
      next() 
    }
    
    ##If gene in pseudo-autosomal region, counts are already correct
    if(gene_pos$strata[gene_pos_line] == "PAR"){
      compiled_counts[[gene_name]] <- unlist(c(sexcounts[gene,], "PAR", "one2one"))
      next()
    }
    
    if(gene_pos$Chr[gene_pos_line] == "chrXIX"){
      line_relationship = which(relationship$Gene_A == paste("gene",gene_name, sep ="-"))
      link = relationship$type[line_relationship]
      
    }else if(gene_pos$Chr[gene_pos_line] == "chrY"){
      line_relationship = which(relationship$Gene_B == paste("gene",gene_name, sep ="-"))
      link = relationship$type[line_relationship]
    }
    
    if(length(line_relationship) >1){
      next()
    }
    
    #if gene not in relationship, skipping (relation could not be established)
    if(identical(line_relationship, integer(0))){
      next()
    }
    
    gene_X = substring(relationship$Gene_A[line_relationship],6,
                       nchar(relationship$Gene_A[line_relationship]))
    gene_Y = substring(relationship$Gene_B[line_relationship],6,
                       nchar(relationship$Gene_B[line_relationship]))
    
    if( link == "one2one" & !(paste(gene_X, gene_Y, sep = "_") %in% names(compiled_counts))){
      X = sexcounts[rownames(sexcounts) == gene_X,]
      Y = sexcounts[rownames(sexcounts) == gene_Y,]
      strata = gene_pos$strata[gene_pos$Gene == relationship$Gene_B[line_relationship]]
      
      if(dim(Y)[1] == 0){
        compiled_counts[[paste(gene_X, gene_Y, sep = "_")]] <- unlist(c(X, strata, "one2one"))
        female.XX[[paste(gene_X, gene_Y, sep = "_")]]<- unlist(c(gene_X, gene_Y, strata, X[,popmap$V2 == "F"]))
        male.X[[paste(gene_X, gene_Y, sep = "_")]]<- unlist(c(gene_X, gene_Y, strata, X[,popmap$V2 == "M"]))
        male.Y[[paste(gene_X, gene_Y, sep = "_")]]<- unlist(c(gene_X, gene_Y, strata, rep(0, sum(popmap$V2 == "M"))))
        next()    
      }
      if(dim(X)[1] == 0){
        compiled_counts[[paste(gene_X, gene_Y, sep = "_")]] <- unlist(c(Y, strata, "one2one"))
        female.XX[[paste(gene_X, gene_Y, sep = "_")]]<- unlist(c(gene_X, gene_Y, strata, rep(0, sum(popmap$V2 == "F"))))
        male.X[[paste(gene_X, gene_Y, sep = "_")]]<- unlist(c(gene_X, gene_Y,strata, rep(0, sum(popmap$V2 == "M"))))
        male.Y[[paste(gene_X, gene_Y, sep = "_")]]<- unlist(c(gene_X, gene_Y,strata, Y[,popmap$V2 == "M"]))
        next()
      }
      
      female.XX[[paste(gene_X, gene_Y, sep = "_")]]<- unlist(c(gene_X, gene_Y, strata, X[,popmap$V2 == "F"]))
      male.X[[paste(gene_X, gene_Y, sep = "_")]]<- unlist(c(gene_X, gene_Y, strata, X[,popmap$V2 == "M"]))
      male.Y[[paste(gene_X, gene_Y, sep = "_")]]<- unlist(c(gene_X, gene_Y,strata, Y[,popmap$V2 == "M"]))
      compiled_counts[[paste(gene_X, gene_Y, sep = "_")]] <- unlist(c(X + Y, strata, link))  
      
    
    
    }else if(link == "lostY" | link == "gainX"){
      if(length(strsplit(relationship$Gene_A[line_relationship], ",")[[1]]) ==1){
        strata = gene_pos$strata[gene_pos_line]
        X = sexcounts[rownames(sexcounts) == gene_X,]
        compiled_counts[[paste(gene_X, sep = "_")]] <- unlist(c(X, strata, link))
      }
      
    }else if(link == "lostX"| link == "gainY"){
      if(length(strsplit(relationship$Gene_B[line_relationship], ",")[[1]]) == 1){
        strata = gene_pos$strata[gene_pos_line]
        Y = sexcounts[rownames(sexcounts) == gene_Y,]
        compiled_counts[[paste(gene_Y, sep = "_")]] <- unlist(c(Y, strata, link))
      }
    
    }else if(link == "No_group"){
      strata = gene_pos$strata[gene_pos_line]
      if (gene_pos$Chr[gene_pos_line] == "chrY"){
        Y = sexcounts[rownames(sexcounts) == gene_Y,]
        compiled_counts[[paste(gene_Y, sep = "_")]] <- unlist(c(Y, strata, paste(link, "Y", sep ="")))
      }else{
        X = sexcounts[rownames(sexcounts) == gene_X,]
        compiled_counts[[paste(gene_X, sep = "_")]] <- unlist(c(X, strata,  paste(link, "X", sep ="")))
      }
      
      
    }else{
      next()}
  }
  compiled_counts <-data.frame(do.call(rbind, compiled_counts))
  colnames(compiled_counts)[meta_coll_n] <- c("strata", "type")
  return(list("compiled_counts" = compiled_counts,
              "female.XX" = female.XX,
              "male.X" = male.X,
              "male.Y" = male.Y))
}

run_bootstrap <- function(data, Strata, Type, Sex, autosome_data, N){ 
  subdat <- data %>%
    filter(type %in% Type,
           strata %in% Strata) %>%
    pull(paste("mean", Sex, sep ="_"))

  if(length(subdat) == 0){return(list("IC" = data.frame("Value" = NA, IC_min = NA, IC_UP = NA),
                                    "bootstrap" = NA,
                                    "Function" = NA,
                                    "N_sexchr" = 0,
                                    "N_auto" = 0))}
  
  return(ACB(data = rbind(cbind(autosome_data, 1),
                                              cbind(subdat,2)),
                                 N.boot = N,
                                 func = auto_sex_ratio))
}

#################
### parameters###
#################
  
N <- 1000 # number of bootstrap replicate
colorscheme <- viridisLite::viridis(n = 2, option = 'H', begin = 0.13, end = 0.85)
names(colorscheme) <- c("M", "F")

###############
##Parse input##
###############
relationship <- read.table("../Gene_relationship_sexchr/05_relationship_table/relationship_table.txt", h= F)
colnames(relationship) <- c("Gene_A", "chr_A", "Gene_B", "chr_B", "type")
region_breakpoints_X <- data.frame(Region = c("PAR", "II", "III", "I"), pos = c(0, 2.5, 6.89, 12.5))
region_breakpoints_Y <- data.frame(Region = c("PAR", "I", "II", "III"), pos = c(0, 0.34, 4.67, 9.67))
gene_pos <- read.table("./02_data/gene_position", h= F)
colnames(gene_pos) <- c("Chr","start","end", "Gene")

#add strata to gene_pos:
gene_pos <- add_strata(gene_pos, region_breakpoints_X, region_breakpoints_Y)

##################
##Results tables##
##################
Res2 <- list()
allele.specific.XY <- list()
SBG_report_sex <- list()

##########################
#Merge each tissue counts#
##########################
for(tiss in c("Cerveau", "Branchie", "Foie")){
  print("running tissue:")
  print(tiss)

  #parsing tissue datasets:
  ########################
  sexcounts <- read.table(paste("./Script/Data", tiss, "04_filtered/sexchr_gene_counts", sep = "/"))
  autocounts<- read.table(paste("./Script/Data/", tiss, "04_filtered/gene_count_filtered", sep = "/"))
  norm_factor <- read.table(paste("./Script/Data/", tiss,"06_differential_expression/norm_factor", sep = "/"))
  popmap <- read.table(paste("./Script/Data/",tiss, "02_metadata/clean_sample_sex", sep = "/"))

  #Sorting datasets:
  ###################
  autocounts <- autocounts[,match(popmap$V1,colnames(autocounts))]
  sexcounts <- sexcounts[,match(popmap$V1,colnames(sexcounts))]
  meta_coll_n <- (c(-1,0) + ncol((sexcounts))+2)
  
  #building gene count for one2one orthologues:
  ############################################
  compiled_counts <- merge_counts(sexcounts, gene_pos, relationship)
  
  #Normalization
  ##############
  tmp = t(compiled_counts$compiled_counts[,-meta_coll_n])
  tmp <- apply(tmp,2,as.numeric)
  Normalized_matrixsummed_chrsex <- log2(data.frame(t(tmp/(2^norm_factor$V2)))  +0.5)
  Normalized_autosomal <- log2(data.frame(t(t(autocounts)/(2^norm_factor$V2))) + 0.5)
  rm(tmp)
  
  #Differential expression analysis:
  ##################################
  wilcoxon <- unlist(apply(Normalized_matrixsummed_chrsex,
                           1,
                           function(x,sex){
                             wilcox.test(x[sex == "M"], x[sex == "F"])$p.value},
                           popmap$V2))
  Mean_F <- apply(Normalized_matrixsummed_chrsex[,popmap$V2 == "F"], 1, mean)
  Mean_M <- apply(Normalized_matrixsummed_chrsex[,popmap$V2 == "M"], 1, mean)
     

  SBG_report_sex[[tiss]] <- data.frame(gene = rownames(Normalized_matrixsummed_chrsex),
                     Mean = (apply(2^Normalized_matrixsummed_chrsex,1, mean)),
                     LFC = Mean_M - Mean_F,
                     Mean_F = Mean_F,
                     Mean_M = Mean_M,
                     p.value = wilcoxon,
                     fdr = p.adjust(wilcoxon, "BH"),
                     strata = unlist(compiled_counts$compiled_counts$strata),
                     type = unlist(compiled_counts$compiled_counts$type),
                     tissue = tiss)
  
  #Calcul of per sex median gene expression per strata and type for dosage compensation
  #####################################################################################
  Normalized_matrixsummed_chrsex$strata= unlist(compiled_counts$compiled_counts$strata)
  Normalized_matrixsummed_chrsex$type = unlist(compiled_counts$compiled_counts$type)

  colnames(Normalized_matrixsummed_chrsex) <- colnames(compiled_counts$compiled_counts)
  #values to be used for all calculation, autosomal median of gene expression

  per_gene_autosomal_medians <- apply(Normalized_autosomal,1, mean)
  Normalized_matrixsummed_chrsex$mean_M <- apply(Normalized_matrixsummed_chrsex[,-meta_coll_n][,popmap$V2 == "M"], 1, mean)
  Normalized_matrixsummed_chrsex$mean_F <- apply(Normalized_matrixsummed_chrsex[,-c(meta_coll_n,meta_coll_n[2]+1)][,popmap$V2 == "F"], 1, mean)

  #Decomposision of allele and genotype gene expression
  #####################################################
  female.XX <- data.frame(do.call(rbind, compiled_counts$female.XX))
  colnames(female.XX)[c(1:3)] <- c("gene.X", "gene.Y", "Strata")
  female.XX[,-c(1:3)] <- apply(female.XX[,-c(1:3)], 2 , as.numeric)
  
  male.X <- data.frame(do.call(rbind, compiled_counts$male.X))
  colnames(male.X)[c(1:3)] <- c("gene.X", "gene.Y", "Strata")
  male.X[,-c(1:3)] <- apply(male.X[,-c(1:3)], 2 , as.numeric)
  
  male.Y <- data.frame(do.call(rbind, compiled_counts$male.Y))
  colnames(male.Y)[c(1:3)] <- c("gene.X", "gene.Y", "Strata")
  male.Y[,-c(1:3)] <- apply(male.Y[,-c(1:3)], 2 , as.numeric)
  
  male.XY <- male.X[,-c(1:3)] + male.Y[,-c(1:3)]
  
  Norm.female.XX <- log2(data.frame(t(t(female.XX[,-c(1,2,3)])/(2^norm_factor$V2[popmap$V2 == "F"])))  +0.5) 
  Norm.male.X <- log2(data.frame(t(t(male.X[,-c(1,2,3)])/(2^norm_factor$V2[popmap$V2 == "M"])))  +0.5) 
  Norm.male.Y <- log2(data.frame(t(t(male.Y[,-c(1,2,3)])/(2^norm_factor$V2[popmap$V2 == "M"])))  +0.5) 
  Norm.male.XY <- log2(data.frame(t(t(male.XY)/(2^norm_factor$V2[popmap$V2 == "M"])))  +0.5)
  
  
  Mean.female.XX <- female.XX[,c(1:3)]
  Mean.female.XX$fdr <- SBG_report_sex[[tiss]]$fdr[SBG_report_sex[[tiss]]$gene %in% paste(female.XX$gene.X, female.XX$gene.Y, sep= "_")|
                                                     SBG_report_sex[[tiss]]$gene %in% female.XX$gene.X]
  Mean.male.X <- male.X[,c(1:3)]
  Mean.male.Y <- male.Y[,c(1:3)]
  Mean.male.XY <- male.X[,c(1:3)]
  
  Mean.male.XY$tissue <- Mean.female.XX$tissue <- Mean.male.X$tissue <- Mean.male.Y$tissue <-  tiss
  
  Mean.female.XX$Mean.F.XX <-  apply(Norm.female.XX, 1, mean)
  Mean.male.X$Mean.M.X <-  apply(Norm.male.X, 1, mean)
  Mean.male.Y$Mean.M.Y <-  apply(Norm.male.Y, 1, mean)
  Mean.male.XY$Mean.M.XY <-  apply(Norm.male.XY, 1, mean)

  #Parsing allele_Specific gene expression
  ########################################
  allele.specific.XY[[tiss]] <- reduce(list(Mean.female.XX, Mean.male.X, Mean.male.Y, Mean.male.XY), .f = merge)
  allele.specific.XY[[tiss]]$Autosomal.mediane <- median(per_gene_autosomal_medians)

  #Dosage compensation bootstrap
  ##############################
  for(Sex in c("M", "F")){
      for(strata in c("all", "PAR", "I", "II", "III"))
        for(type in c("all", "one2one", "lostY", "lostX")){
          if (type == "lostX" & Sex == "F"){next}
          if (type ==  "all"){
            if(Sex == "F"){Type <- c("one2one", "lostY", "No_groupX")}
            if(Sex == "M"){Type <- c("one2one", "lostY", "lostX", "No_groupX", "No_groupY")}
            }else{
              if(type == "lostY"){Type <- c("lostY", "No_groupX")}
              else if(type == "lostX"){Type <- c("lostX", "No_groupY")}
              else{Type <- type}

            }
          if (strata ==  "all"){
            Strata <- c("PAR", "I", "II", "III")}else{
              Strata <- strata
            }
          boot.res <- run_bootstrap(data = Normalized_matrixsummed_chrsex,
                                    Type = Type,
                                    Strata = Strata,
                                    Sex = Sex,
                                    autosome_data = per_gene_autosomal_medians, N= 1000)
          Res2[[paste(tiss,Sex,strata,type, sep ="")]] <- data.frame(
            rbind(c(tiss,
                    type,
                    Sex,
                    strata,
                    unlist(boot.res$IC),
                    boot.res$N_sexchr)))
      }
  }
    }

##Prepare results for plot
Res2 <- do.call(rbind, Res2)
colnames(Res2) <- c("Tissue", "type", "Sex", "Strata",
                    "Median X(X/Y):AA","IC_low", "IC_high", "N")
Res2$`Median X(X/Y):AA` <- as.numeric(Res2$`Median X(X/Y):AA`)
Res2$IC_low <- as.numeric(Res2$IC_low)
Res2$IC_high <- as.numeric(Res2$IC_high)

#plots
###################### 
#Dosage_compensation##
######################
  
Res2$Strata <- factor(Res2$Strata, levels = c("all", "PAR","I", "II", "III"))
Res2$Tissue <- factor(Res2$Tissue, levels = c("Cerveau", "Branchie", "Foie"))
Res2$type <- factor(Res2$type, levels = c("all", "one2one", "lostY", "lostX"))

Figure_5 <- Res2 %>%
  na.omit(.) %>%
    ggplot(aes(x = type, y = `Median X(X/Y):AA`, color = Sex)) +
  geom_point(size = 1.5,position=position_dodge(0.5))+
  facet_grid(Tissue ~Strata, scale = "free") +
  geom_linerange(aes(ymin = IC_high, ymax = IC_low,x = type, ),
                 linewidth = 0.6,
                 position=position_dodge(0.5)) +
  geom_hline(yintercept = 0) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, vjust = 0, hjust = 0.2)) +
  scale_color_manual(values = colorscheme)

    ggsave("99_FiguresFigure_5.png",
       plot = Figure_5,
       device = "png",
       width = 17.9,
       height = 12,
       units = "cm",
       dpi = 300)

    ggsave("99_Figures/Figure_5.pdf",
       plot = Figure_5,
       device = "pdf",
       width = 17.9,
       height = 10,
       units = "cm",
       dpi = 300)

  write.table(Res2,
              "99_Figures/Table_S9",
              row.names = F,
              col.names = T,
              quote = F)
  
#########################################
#######Allele specific gene expression###
#########################################
  
#Sex-bias gene expression
allele.specific.XY <-do.call(rbind, allele.specific.XY) %>% data.frame()
allele.specific.XY$LFC_sex <- allele.specific.XY$Mean.F.XX - allele.specific.XY$Mean.M.XY
allele.specific.XY$LFC.XY <- allele.specific.XY$Mean.F.X - allele.specific.XY$Mean.M.Y

allele.specific.XY$group_LFC_sex <- factor(NA, levels = 1:3, )
for(i in 1:nrow(allele.specific.XY)){
  if(allele.specific.XY$LFC_sex[i] < -0.5 ){
    allele.specific.XY$group_LFC_sex[i] <- 1
    next}
  if(allele.specific.XY$LFC_sex[i] <= 0.5){
    allele.specific.XY$group_LFC_sex[i] <- 2
    next}
  allele.specific.XY$group_LFC_sex[i] <- 3
  
}

allele.specific.XY.long <- pivot_longer(allele.specific.XY,
                                        cols = c(Mean.M.X, Mean.M.XY, Mean.F.XX, Mean.M.Y),
                                        values_to = "Mean.read.count",
                                        names_to = "Allele")
allele.specific.XY.long$Allele <- factor(allele.specific.XY.long$Allele,
                                         c("Mean.F.XX", "Mean.M.XY", "Mean.M.X", "Mean.M.Y"))
Figure_6 <- allele.specific.XY.long %>% 
  ggplot() +
  geom_boxplot(aes(y = (Mean.read.count), x = group_LFC_sex , fill = Allele),
               lwd = 0.1,
               width = 0.5,
               position = position_dodge(0.8),
               outlier.shape = NA) + 
  theme_bw()+
  scale_x_discrete(labels = c("Male-biased", "unbiased", "Female-biased"))+
  facet_grid(.~tissue) +
  scale_fill_hue(labels = c("XX", "XY", "Male X", "Male Y")) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.text.x = element_text(angle = -45,vjust = 0, hjust = 0.2)) +
  xlab("Log2 fold change") +
  ylab("Mean expression (log)")

    ggsave("99_Figures/Figure_6.png",
       plot = Figure_6,
       device = "png",
       width = 8.9,
       height = 6,
       units = "cm",
       dpi = 300)

    ggsave("99_Figures/Figure_6.pdf",
       plot = Figure_6,
       device = "pdf",
       width = 17.9,
       height = 10,
       units = "cm",
       dpi = 300)

#####STATISTICAL TESTING
Res <- list()
for(tiss in c("Cerveau", "Branchie", "Foie")){
  Mean.M.X.MBG.soft <- allele.specific.XY.long %>%
    filter(tissue == tiss, Allele == "Mean.M.X", group_LFC_sex  == 1) %>% pull(Mean.read.count)
  Mean.M.Y.MBG.soft <- allele.specific.XY.long %>%
    filter(tissue == tiss, Allele == "Mean.M.Y", group_LFC_sex  == 1) %>% pull(Mean.read.count)
  Mean.M.XY.MBG.soft <- allele.specific.XY.long %>%
    filter(tissue == tiss, Allele == "Mean.M.XY", group_LFC_sex  == 1) %>% pull(Mean.read.count)
  Mean.F.XX.MBG.soft <- allele.specific.XY.long %>%
    filter(tissue == tiss, Allele == "Mean.F.XX", group_LFC_sex  == 1) %>% pull(Mean.read.count)

  Mean.M.X.none <- allele.specific.XY.long %>%
    filter(tissue == tiss, Allele == "Mean.M.X", group_LFC_sex  == 2) %>% pull(Mean.read.count)
  Mean.M.Y.none <- allele.specific.XY.long %>%
    filter(tissue == tiss, Allele == "Mean.M.Y", group_LFC_sex  == 2) %>% pull(Mean.read.count)
  Mean.M.XY.none <- allele.specific.XY.long %>%
    filter(tissue == tiss, Allele == "Mean.M.XY", group_LFC_sex  == 2) %>% pull(Mean.read.count)
  Mean.F.XX.none <- allele.specific.XY.long %>%
    filter(tissue == tiss, Allele == "Mean.F.XX", group_LFC_sex  == 2) %>% pull(Mean.read.count)

  Mean.M.X.FBG.soft <- allele.specific.XY.long %>%
    filter(tissue == tiss, Allele == "Mean.M.X", group_LFC_sex  == 3) %>% pull(Mean.read.count)
  Mean.M.Y.FBG.soft <- allele.specific.XY.long %>%
    filter(tissue == tiss, Allele == "Mean.M.Y", group_LFC_sex  == 3) %>% pull(Mean.read.count)
  Mean.M.XY.FBG.soft <- allele.specific.XY.long %>%
    filter(tissue == tiss, Allele == "Mean.M.XY", group_LFC_sex  == 3) %>% pull(Mean.read.count)
  Mean.F.XX.FBG.soft <- allele.specific.XY.long %>%
    filter(tissue == tiss, Allele == "Mean.F.XX", group_LFC_sex  == 3) %>% pull(Mean.read.count)

Res[[tiss]] <-rbind(
c(tiss, "MBG", "XX", wilcox.test(Mean.F.XX.MBG.soft,Mean.F.XX.none)$p.value),
c(tiss, "FBG", "XX", wilcox.test(Mean.F.XX.FBG.soft,Mean.F.XX.none)$p.value),
c(tiss, "MBG", "X", wilcox.test(Mean.M.X.MBG.soft,Mean.M.X.none)$p.value),
c(tiss, "MBG", "Y", wilcox.test(Mean.M.Y.MBG.soft,Mean.M.Y.none)$p.value),
c(tiss, "MBG", "XY", wilcox.test(Mean.M.XY.MBG.soft,Mean.M.XY.none)$p.value),
c(tiss, "FBG", "X", wilcox.test(Mean.M.X.FBG.soft,Mean.M.X.none)$p.value),
c(tiss, "FBG", "Y", wilcox.test(Mean.M.Y.FBG.soft,Mean.M.Y.none)$p.value),
c(tiss, "FBG", "XY", wilcox.test(Mean.M.XY.FBG.soft,Mean.M.XY.none)$p.value))
}
Res <- do.call(rbind, Res)
write.table(Res,
            "99_Figures/Figure_6_signif_table",
            quote = F,
            row.names =F,
            col.names =T)

