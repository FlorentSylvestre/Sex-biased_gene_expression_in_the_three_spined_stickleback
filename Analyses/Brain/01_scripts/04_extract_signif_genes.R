##Take a DEG report
##Extract signif genes toward female and male for subsequent GO ontology
##Rscript extract_signif_gene.R
##########################################################################################################

######################
##Library:############
######################
library(data.table)
library(tidyverse)

#####################
##Function###########
#####################

#####################
##Parse User input:##
#####################
report_path = "06_differential_expression/report.txt"
outpath = "06_differential_expression/DEG.txt"

#####################
##Reading input######
#####################
report = fread(report_path) %>%
                as_tibble()

#####################
##Output#############
#####################

report %>%
        dplyr::filter(signif != "FALSE") %>%
        select(gene) %>%
        write.table(.,
                    outpath,
                    quote = F,
                    col.names = F,
                    row.names = F)
