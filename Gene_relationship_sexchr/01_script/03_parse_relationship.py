#!/bin/envs python3
"""
Script that parse a gff file and Orthofinder output to find ortholous genes between Y and X
output one file with gene orthologues, one with genes with no match on one or the other chromosome
and genes for which transcript are not clustered together
Usage: parse_sex_orthologues.py <gff_file> <path_to_orthofinder_result>
Parth to Orthofinder result exemple: "04_protein/Orthofinder/Results_Oct27"
gff_file correspond to the focal species gff
"""

#import modules
from collections import defaultdict
import sys

#Constants
COMPCHROM = {"chrY", "chrXIX"}

##Parsing input
gff_path = sys.argv[1]
orthologue_path = sys.argv[2].strip("/") + "/Orthogroups/Orthogroups.txt"
out_path = 05_relationship_table/relationship_table.txt

##functions

def parsing_gff(gff_path):
    tmp = defaultdict(str)
    second_pass = defaultdict(list)
    gff_dict_gene = defaultdict(lambda: defaultdict(list))
    gff_dict_tr = defaultdict(lambda: defaultdict(str))

    with open(gff_path, "r") as gff:
        for line in gff:
            if line.startswith("#"):
                continue

            l = line.strip().split("\t")
            seq_type = l[2]
            name = l[8].split(";")[0][3:]
            tmp[name] = seq_type
            if seq_type in ["gene", "pseudogene", "region","cDNA_match", "D_loop"]:
                continue

            parent = l[8].split(";")[1][7:]
            chrom  = l[0]

            if parent not in tmp:
               second_pass[name] = [name, parent, chrom, seq_type]

            if tmp[parent] == "gene":
                gff_dict_gene[parent]["tr"].append(name)
                gff_dict_gene[parent]["chrom"].append(chrom)
                gff_dict_tr[name]["gene"] = parent
                gff_dict_tr[name]["chrom"] = chrom
    return gff_dict_gene, gff_dict_tr


def parsing_orthogroup(orthogroup_file):
    gene_per_ortho = defaultdict(list)
    ortho_info = defaultdict(str)
    count_NA = 0
    with open(orthogroup_file, "r") as ortho:
        ortho.readline() #skipping header
        for line in ortho:
            group, _, stickle = line.strip("\n").split("\t")

            if stickle == "":
                continue

            for seq in stickle.split(", "):
                gene_per_ortho[group].append(seq)
                ortho_info[seq.split(" ")[0]] = group
    return gene_per_ortho, ortho_info

def no_group_transcript(gff_tr, tr_group):
    no_group = dict()
    for tr in gff_tr:
        if tr in tr_group:
            continue

        if gff_tr[tr]["chrom"] not in ["chrY", "chrXIX"]:
            continue

        if gff_tr[tr]["gene"] not in no_group:
            no_group[gff_tr[tr]["gene"]] = gff_tr[tr]["chrom"]
    return no_group


def list_transcript_cluster(gff, transcript_group):
    dispersed_genes = []
    for gene in gff:
        groups = [transcript_group[x] for x in gff[gene]["tr"] if x in transcript_group]

        if len(set(groups)) > 1:
            dispersed_genes.append(gene)

    return dispersed_genes


def check_dispersed(gene, gene_list):
    return gene in gene_list


def check_ortho(orthogroup,gff_tr):
    genes = []
    chrom = []

    for tr in orthogroup:
        if tr.split(" ")[0] in gff_tr:
            genes.append(gff_tr[tr.split(" ")[0]]["gene"])
            chrom.append(gff_tr[tr.split(" ")[0]]["chrom"])


    if set(chrom) != COMPCHROM:
        if set(chrom) == {"chrXIX"}:
            if len(set(genes)) >1:
                return ["Many2none", f"{','.join(set(genes))}\tchrX\tNone\tNone\tMany2none"]
            return ["lostY", f"{','.join(set(genes))}\tchrX\tNone\tNone\tlostY"]

        if set(chrom) == {"chrY"}:
            if len(set(genes)) >1:
                return ["None2many", f"None\tNone\t{','.join(set(genes))}\tchrY\tNone2many"]
            return ["lostX", f"None\tNone\t{','.join(set(genes))}\tchrY\tlostX"]

        if "chrY" in chrom and "chrXIX" in chrom:
            genes_Y = set([y for x,y in  enumerate(genes) if chrom[x] == 'chrY'])
            genes_X = set([y for x,y in  enumerate(genes) if chrom[x] == 'chrXIX'])

            return ["Complex", f"{','.join(set(genes_X))}\tchrX\t{','.join(set(genes_Y))}\tchrY\tComplex"]

        if "chrY" in set(chrom):
            genes_y = set([y for x,y in  enumerate(genes) if chrom[x] == 'chrY'])
            if len(genes_y) >1:
                return ["gain_Y_many", f"None\tNone\t{','.join(genes_y)}\tchrY\tgain_Y_many"]
            return ["gainY", f"None\tNone\t{','.join(genes_y)}\tchrY\tgainY"]

        if "chrXIX" in set(chrom):
            genes_X = set([y for x,y in  enumerate(genes) if chrom[x] == 'chrXIX'])

            if len(genes_X) >1:
                return ["Gain_X_many", f"{','.join(set(genes_X))}\tchrX\tNone\tNone\tGain_X_many"]
            return ["gainX", f"{','.join(genes_X)}\tchrX\tNone\tNone\tgainX"]

        return ["No orthologs", None]

    X = [y for x,y in enumerate(genes) if chrom[x] == "chrXIX"]
    Y = [y for x, y in enumerate(genes) if chrom[x] == "chrY"]
    line = f"{','.join(set(X))}\tchrX\t{','.join(set(Y))}\tchrY"
    line += "\tmany2" if len(set(X)) >=2 else "\tone2"
    line += "many" if len(set(Y)) >=2 else "one"
    return [line.split("\t")[-1], line]

def write_log(expt_table, out):
    log = open(out + ".log", "w")
    for chrom in expt_table:
        for expt in expt_table[chrom]:
            log.write(f"{expt}\t{chrom}\t{expt_table[chrom][expt]}\n")
    log.close()



##main
gff_gene, gff_tr = parsing_gff(gff_path)
group_tr, tr_group = parsing_orthogroup(orthologue_path)
dispersed_genes = list_transcript_cluster(gff_gene, tr_group)
no_group = no_group_transcript(gff_tr, tr_group)


out = open(out_path, "w")
exception_table = defaultdict(lambda: defaultdict(int))

for gene in no_group:
    if gene in dispersed_genes:
        continue

    if no_group[gene] == "chrY":
        out.write(f"None\tNone\t{gene}\tchrY\tNo_group\n")
    else:
        out.write(f"{gene}\tchrX\tNone\tNone\tNo_group\n")
    exception_table[no_group[gene]]["no_group"] += 1

for gene in dispersed_genes:
    exception_table[gff_gene[gene]["chrom"][0]]["Dispersed"] += 1

    if gff_gene[gene]["chrom"][0] == "chrY":
        out.write(f"None\tNone\t{gene}\tchrY\tDispersed\n")
    if gff_gene[gene]["chrom"][0] == "chrXIX":
        out.write(f"{gene}\tchrX\tNone\tNone\tDispersed\n")


for ortho in group_tr:
    dispersed = [check_dispersed(gff_tr[tr.split(" ")[0]]["gene"], dispersed_genes) for tr in group_tr[ortho]]

    if any(dispersed):
        genes = set([gff_tr[tr.split(" ")[0]]["gene"] for x, tr in enumerate(group_tr[ortho]) if not dispersed[x]])
        for g in genes:
            if gff_gene[g]["chrom"][0] == "chrY":
                out.write(f"None\tNone\t{g}\tchrY\tGrouped_with_dispersed\n")
            if gff_gene[g]["chrom"][0] == "chrXIX":
                 out.write(f"{g}\tchrX\tNone\tNone\tGrouped_with_dispersed\n")

            exception_table[gff_gene[g]["chrom"][0]]["Grouped_with_Dispersed"] += 1

        continue

    relationship, line = check_ortho(group_tr[ortho], gff_tr)
    genes_in_group = []
    for tr in group_tr[ortho]:
        genes_in_group.append(gff_tr[tr.split(" ")[0]]["gene"])

    for gene in set(genes_in_group):
        exception_table[gff_gene[gene]["chrom"][0]][relationship] += 1


    if line is None:
        continue

    out.write(line + "\n")
out.close()
write_log(exception_table, out_path)
