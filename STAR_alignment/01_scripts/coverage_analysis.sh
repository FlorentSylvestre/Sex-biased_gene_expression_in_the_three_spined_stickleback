#!/bin/bash


path="04__raw_data"
out="Coverage.txt"

for file in $(ls $path)
    do
        filename=$(echo $file | cut -d"." -f5)
        readscounts=$(echo $(zcat $file | wc -l)/4|bc)
        echo $filename  $readscounts
    done>$out
