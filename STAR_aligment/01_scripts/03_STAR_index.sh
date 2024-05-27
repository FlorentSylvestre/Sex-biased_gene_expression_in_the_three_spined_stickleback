#!/bin/sh

#uses samtools and star to generate required indexing of refence genome
#Reference genome should be named genome_F and genome_M
#respectively for each sex


########################################
## Global variables
########################################

GENOME="03_reference/genome"
GFF=$2
NCPUS=6

########################################
## Load environment ##
########################################
module load samtools
module load star/2.7.2b
module load java/jdk/1.8.0_102
CSD="/prg/picard-tools/1.119/CreateSequenceDictionary.jar"


########################################
## Execute script ##
########################################

samtools faidx ${GENOME}

#Build reference
echo "Building ${GENOME} index..."

if [ ! -e  ${GENOME%.*}.dict ]
then
    java -jar "$CSD" \
    R="$GENOME"_F.fa \
    O="$GENOME"_F.dict

    java -jar "$CSD" \
    R="$GENOME"_M.fa \
    O="$GENOME"_M.dict
fi

BIT_SIZE=20
echo "STAR indexing \n"
