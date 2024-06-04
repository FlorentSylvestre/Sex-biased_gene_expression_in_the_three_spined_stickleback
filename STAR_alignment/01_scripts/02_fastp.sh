#!/bin/bash

#trim single end fastq file from QUANTseq SE 50bp sequencing
#02_fastp.sh <file_listing_all_fastp>


##Variable
NCPU=6
#Suggested memory: 10G
List_file=$1

while read line
do
    fastp -w $NCPU\
    -i $line \
    -o 06_fastp/$(basename $line) \
    --length_required=20 \
    -f 12\
    --trim_poly_g \
    --poly_g_min_len=6 \
    --trim_poly_x \
    --poly_x_min_len=7 \
    --json 06_fastp/$(basename ${line/.fastq.gz/}).json\
    --html 06_fastp/$(basename ${line/.fastq.gz/}).html;
done<$List_file
