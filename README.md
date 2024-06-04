This Github repository explain how to reproduce analysis presented in 
Title: Sex-biased gene expression across tissues reveals unexpected differentiation in the gills of the threespine stickleback
Authors: Florent Sylvestre, Nadia Aubin-Horth, Louis Bernatchez,
which is currently unpublished.

Corresponding author: Florent Sylvestre, flosylv@hotmail.fr
Don't hesitate to contact me if you have questions

[![DOI](https://zenodo.org/badge/787045413.svg)](https://zenodo.org/doi/10.5281/zenodo.11477975)



It consist of 4 main folders:

1) STAR_alignment contain all the script needed to reproduce alignement to the threespine stickleback reference genome and gene expression quantification for each tissue
2) Genome annotation contain all the script used to link each annotated genes in the threespined stickleback genome to an Uniprot entry, and extract gene function and ontology
3) Gene_relationship_sexchr contain all the script used to identify gene still sharing (or not) a copy between sex-chromosome of the threespine sticleback

4)Analyses makes use of results from all three precedent folder to reproduce analyses presented in the manuscript.

Finally, the Utils forlder contain some script usefull for converting data between format or removing pseudo-autosomal region from our reference

Detailled information for each part is available in correspondant subfolder


Usefull information:
All analyses are base of the V5 version of the threespine stickleback and its Y reference sequence.
it is available for download at https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016920845.1/
Used Annotation file and transcriptome are also available there.

Dependencie:
Software:
blastx
Orthofinder v2.5.5
STAR v2.7.2b
HTSeq v0.11.3
FastQC v0.11.8
fastp v0.15.0
samtools  v1.10

Databases:
Swissprot Database
Threespine stickleback reference genome: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016920845.1/
Ninespine stickleback reference genome: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_949316345.1/

R session infos:
R version 4.3.2 (2023-10-31)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods
[8] base

other attached packages:
 [1] DESeq2_1.40.2               SummarizedExperiment_1.30.2
 [3] Biobase_2.60.0              MatrixGenerics_1.12.3
 [5] matrixStats_1.2.0           GenomicRanges_1.52.1
 [7] GenomeInfoDb_1.36.4         IRanges_2.36.0
 [9] S4Vectors_0.40.2            BiocGenerics_0.48.1
[11] data.table_1.14.10          lubridate_1.9.3
[13] forcats_1.0.0               stringr_1.5.1
[15] dplyr_1.1.4                 purrr_1.0.2
[17] readr_2.1.4                 tidyr_1.3.0
[19] tibble_3.2.1                ggplot2_3.4.4
[21] tidyverse_2.0.0

loaded via a namespace (and not attached):
 [1] utf8_1.2.4              generics_0.1.3          bitops_1.0-7
 [4] lattice_0.22-5          stringi_1.8.3           hms_1.1.3
 [7] magrittr_2.0.3          grid_4.3.2              timechange_0.2.0
[10] Matrix_1.6-3            fansi_1.0.6             scales_1.3.0
[13] codetools_0.2-19        abind_1.4-5             cli_3.6.2
[16] rlang_1.1.3             crayon_1.5.2            XVector_0.40.0
[19] munsell_0.5.0           DelayedArray_0.26.7     withr_2.5.2
[22] S4Arrays_1.0.6          parallel_4.3.2          tools_4.3.2
[25] BiocParallel_1.34.2     tzdb_0.4.0              colorspace_2.1-0
[28] locfit_1.5-9.8          GenomeInfoDbData_1.2.10 vctrs_0.6.5
[31] R6_2.5.1                lifecycle_1.0.4         zlibbioc_1.46.0
[34] pkgconfig_2.0.3         pillar_1.9.0            gtable_0.3.4
[37] Rcpp_1.0.12             glue_1.7.0              tidyselect_1.2.0
[40] compiler_4.3.2          RCurl_1.98-1.14


