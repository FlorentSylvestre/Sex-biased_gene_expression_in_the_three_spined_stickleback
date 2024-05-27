#!/bin/sh
gff=$1
genome=$2
out=04_protein/$(basename $genome)
gffread $gff -g $genome -C -V -H -y $out
