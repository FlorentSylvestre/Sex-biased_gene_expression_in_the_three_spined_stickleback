#!/bin/sh

#Run star fo each individual
#made to be run by parallel
#script <Sex [M|F]>, <input.fasta>, <NCPU>

#Modules
#Manitou
module load star/2.7.2b

#variables:
GENOME="03_reference/genome"
OUTPUT="07_star_alignment"

SEX="$1"
INPUT="$2"
NCPU="$3"

#sample_file
sample="$(ls 06_fastp/*${INPUT}_*fastq.gz)"
echo "aligning sample $sample"
zcat $sample >tmp_$INPUT

#Scripts:
STAR --runThreadN ${NCPU} \
    --genomeDir 03_reference/STARindex_$SEX \
    --readFilesIn tmp_$INPUT \
    --readFilesCommand -\
    --twopassMode None \
    --outSAMmapqUnique 37 \
    --outFileNamePrefix ${OUTPUT}/$SEX/${INPUT} \
    --outSAMtype BAM Unsorted \
    --outFilterMultimapNmax 10 \
    --outFilterMismatchNmax 10 \
    --outSAMattrRGline ID:$INPUT SN:$INPUT PL:ILLUMINA

rm tmp_$INPUT
