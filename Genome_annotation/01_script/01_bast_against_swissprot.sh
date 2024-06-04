#!/bin/bash

# Global variables
SEQUENCE_FILE=02_sequences/transcriptome.fa
SWISSPROT_RESULT=03_blast_results/analyzed_genes.swissprot
SWISSPROT_HITS=03_blast_results/analyzed_genes.hits
SWISSPROT_DB=99_swissprot_db/swissprot

# Blast all sequences against swissprot (must be installed locally)
# WARNING use `-j N` if you need to limit the number of CPUs used to N <integer>
cat $SEQUENCE_FILE | parallel -k --block 1k --recstart '>' --pipe 'blastx -db '$SWISSPROT_DB' -query - -evalue 1e-3 -outfmt 6 -max_target_seqs 1' > $SWISSPROT_RESULT

# TODO filter blasts on similarity, evalue, length...

# Extract analyzed_genes.hits
awk '{print $1,$2}' $SWISSPROT_RESULT | uniq > $SWISSPROT_HITS
