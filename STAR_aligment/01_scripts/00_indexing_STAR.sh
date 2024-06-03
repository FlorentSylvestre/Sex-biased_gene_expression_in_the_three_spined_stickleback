#!/bin/sh
########################################
## Global variables
########################################

GENOME="03_reference/genome"
GFF=03_reference/Stickleback
NCPUS=6
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
STAR --runMode genomeGenerate \
    --runThreadN ${NCPUS} \
    --genomeDir "03_reference/STARindex_F/" \
    --genomeFastaFiles ${GENOME}_F.fa \
    --sjdbGTFfile ${GFF}_noY.gff \
    --sjdbGTFtagExonParentTranscript Parent \
    --sjdbGTFfeatureExon exon \
    --genomeChrBinNbits ${BIT_SIZE} \
    --genomeSAindexNbases 13

STAR --runMode genomeGenerate \
    --runThreadN ${NCPUS} \
    --genomeDir "03_reference/STARindex_M/" \
    --genomeFastaFiles ${GENOME}_M.fa \
    --sjdbGTFfile ${GFF}_NOPAR.gff \
    --sjdbGTFtagExonParentTranscript Parent \
    --sjdbGTFfeatureExon exon \
    --genomeChrBinNbits ${BIT_SIZE} \
    --genomeSAindexNbases 13
