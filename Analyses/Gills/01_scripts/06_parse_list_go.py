#!/bin/python3

##Script that parse a list of GO term and retun all gene ID for those GO
##Using gene annotation to provide alias for generic gene names such as "LOCXXXXXXXX"
##Usage: <script.py> <go_list_file> <output>

# imports
import sys
from collections import defaultdict

#Constant
ANNOTATION_FILE = "03_data/GFF_gene_annotation2.txt"
GO_FILE = "07_GO_analysis/ALL_DEG_GO.tsv"

#user input:
GO_list = sys.argv[1]
output = sys.argv[2]

##parsing gene-alias relation
gene_alias_conv = defaultdict(str)

with open(ANNOTATION_FILE, "r") as annot:
    for line in annot:
        l = line.strip().split("\t")
        gene_alias_conv[l[0]] = l[6]


#reading GO_list
GO_list = [x.strip() for x in open(GO_list, "r").readlines()]


#outputing selected gene lists:
with open(GO_FILE, "r") as go:
    with open(output, "w") as outp:

        for line in go:
            l = line.strip().split("\t")
            if l[0] in GO_list:
                for gene in l[13].split(", "):
                    outp.write(l[0] + "\t" + l[3] + "\t" + gene + "\t" + gene_alias_conv[gene] + "\n")

