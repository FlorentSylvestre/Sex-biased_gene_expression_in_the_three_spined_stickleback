#!/usr/bin/env bash
#Run star fo each individual
#made to be run by parallel
#script <Sex [M|F]>, <input.fasta>, <NCPU>
########################################
## Global variables
########################################

# Global variables
FIRST_PASS="07_star_alignment"
OUTPUT="08_star_pass2"

SEX=$1
INPUT="$2"
NCPUS="$3"


########################################
## Load environment ##
########################################
#manitou
module load star/2.7.2b
module load samtools/1.15
#valeria:
#module load nixpkgs/16.09 gcc/7.3.0 intel/2018.3 star/2.7.3a #samtools/1.15

########################################
## Prepare sample file
########################################
#sample_file
sample="$(ls 06_fastp/*${INPUT}_*fastq.gz)"
echo $sample
echo "aligning sample $sample"
zcat $sample >tmp_$INPUT



########################################
## Execute script ##
########################################

#remove phase1 bam files
echo "aligning $INPUT"

STAR --runThreadN ${NCPUS} \
    --genomeDir 03_reference/STARindex_$SEX/ \
    --readFilesIn tmp_$INPUT\
    --readFilesCommand - \
    --twopassMode None \
    --sjdbFileChrStartEnd $FIRST_PASS/$SEX/*SJ.out.tab \
    --limitSjdbInsertNsj=2000000 \
    --outSAMmapqUnique 37 \
    --outFileNamePrefix ${OUTPUT}/$SEX/${INPUT} \
    --outSAMtype BAM Unsorted \
    --outFilterMultimapNmax 10 \
    --outFilterMismatchNmax 10 \
    --outFilterMismatchNoverLmax 0.1 \
    --quantMode TranscriptomeSAM GeneCounts \
    --outSAMattributes NH HI AS nM \
    --outSAMattrRGline ID:${INPUT} SM:${INPUT} PL:ILLUMINA

rm tmp_$INPUT
#filter reads
#Valeria
#module load StdEnv/2020 samtools/1.15.1

samtools view -h -q 33 "$OUTPUT"/"$SEX"/"$INPUT"Aligned.out.bam > "$OUTPUT"/"$SEX"/"$INPUT"Aligned.uniq.bam

samtools sort -@ $NCPUS "$OUTPUT"/"$SEX"/"$INPUT"Aligned.uniq.bam > "$OUTPUT"/"$SEX"/"$INPUT"Aligned.uniq.sorted.bam
rm "$OUTPUT"/$SEX/"$INPUT"Aligned.out.bam
#
samtools index "$OUTPUT"/"$SEX"/"$INPUT"Aligned.uniq.sorted.bam
