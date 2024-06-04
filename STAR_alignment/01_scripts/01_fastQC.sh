#!/bin/sh

# This script run fastQC on each files in the listfile
#Nthread: number of CPU tu use
# Output: outputwhere to write data
#listfile: txt file listing all fasta to process

#Global variable:
listfile="$(cat $1|tr '\n' ' ')"
Nthread=$2
Output="05_fastqc"


#loading fastQC module:
module load fastqc/0.11.8

#running fastQC:
fastqc -o $Output -t $Nthread $listfile
