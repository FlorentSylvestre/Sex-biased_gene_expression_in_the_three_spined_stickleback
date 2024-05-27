#!/bin/python3

"""
Remove the Pseudo-autosomal region from the Y chromosome GFF annotation 
#Shit coordinates by Size and remove element with negative start
#script <Size> <output>

"""
#modules
import sys

#Variables:
GFFs_path = "03_reference/GCF_016920845.1_GAculeatus_UGA_version5_genomic_chr.gff"
PAR_length = sys.argv[1]
output_path = sys.argv[2]

#main script:

with open(GFFs_path, "r") as GFF:
    with open(output_path, "w") as output:

        for lines in GFF:
            if not lines.startswith("Y"):
                output.write(lines)
                continue

            if lines.startswith("Y"):
                line = lines.split("\t")
                if line[2] == "region":
                    line[4] = str(int(line[4]) - PAR_length)
                    output.write("\t".join(line))
                    continue

                new_start = int(line[3]) - PAR_length
                new_end = int(line[4]) - PAR_length

                if new_start < 1:
                    continue
                if new_start >= 1:
                   line[3] = str(new_start)
                   line[4] = str(new_end)
                   output.write("\t".join(line))
