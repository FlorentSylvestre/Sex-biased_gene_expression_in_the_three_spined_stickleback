#!/bin/sh

#Run HTSeq in union mode
#script <Sex [F|M]>, <gff_file> <sampleName>

########################################
## Global variables
########################################

#Global variables
INPUT="08_star_pass2"
OUTPUT="09_gene_counts"

#######################################
## User input
#######################################
sex=$1
gff_file=$2
base=$3

########################################
## Execute script ##
########################################

# for gene expression
htseq-count -f bam \
    -s no \
    -r pos \
    -t "gene" \
    --mode "union" \
    -i "Name" \
    "$INPUT"/$sex/"$base"Aligned.uniq.sorted.bam \
    "$gff_file" >"$OUTPUT"/htseq-count_gene_"$base".txt
