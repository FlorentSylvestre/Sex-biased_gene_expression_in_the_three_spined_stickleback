#!/bin/python3
##Script that annotate gene based on gawn annotation table and gff_file
##Usage: python3 annotate_gene.py annotation_table GCF_file output

#Modules:
import sys
from collections import defaultdict
import re

#class def:
class gene_infos:
    def __init__(self, gene_inf):
        self.name = gene_inf[0]
        self.id = gene_inf[1]
        self.chrom = gene_inf[2]
        self.start = gene_inf[3]
        self.end = gene_inf[4]
        self.annot = defaultdict(list)

    def parse_annot(self,tr):
        self.annot["ID"].append(tr["ID"][0])
        if "acc" in tr:
            if tr["acc"] not in self.annot["acc"]:
                self.annot["acc"].append(tr["acc"])
                self.annot["GN"].append(tr["GN"])
        [self.annot["GO"].append(x) for x in tr["GO"] if x not in self.annot["GO"]]

    def get_id(self):
        return self.id

    def __str__(self):
        if "ID" not in self.annot:
            self.annot["ID"].append("-")
        if "GO" not in self.annot:
            self.annot["GO"].append("-")
        if "acc" not in self.annot:
            self.annot["acc"].append("-")
        if "GN" not in self.annot:
            self.annot["GN"].append("-")
        self.annot["GN"] = ";".join(self.annot["GN"])
        return (
            f"{self.name}\t"
            f"{self.chrom}\t"
            f"{self.start}\t"
            f"{self.end}\t"
            f"{';'.join(self.annot['ID'])}\t"
            f"{';'.join(self.annot['acc'])}\t"
            f"{';'.join(set([x.lower() for x in self.annot['GN'].split(';')]))}\t"
            f"{';'.join(self.annot['GO'])};"
        )


#parsing user input:
accession_path = sys.argv[1]
GO_path = sys.argv[2]
GFF3_path = sys.argv[3]
output_path = sys.argv[4]

#FIXED
annotation_path = "05_annotations/genbank_info/"

#parsing GFF3:
transcript = defaultdict(lambda: defaultdict(list))
gene = defaultdict(gene_infos)

with open(GFF3_path, "r") as GFF:
    for line in GFF:
        if line.startswith("#"):
            continue


        l = line.strip().split("\t")
        if l[2] in ["gene", "pseudogene"]:
            tmp_dict = defaultdict(str)
            for x in l[8].split(";"):
                tmp_dict[x.split("=")[0]] = x.split("=")[1]
            gene_name = tmp_dict["Name"]
            gene_id = tmp_dict["ID"].split("gene-")[1]
            gene[gene_name] = gene_infos([gene_name,gene_id,l[0],l[3],l[4]])
            if gene_name == "rns":
                print(gene[gene_name])
                print(gene[gene_name].get_id())
        if l[2] in ["transcript","mRNA","lnc_RNA","guide_RNA","snoRNA", "snRNA","rRNA"]:
           transcript_id = l[8].split(";")[0].split("rna-")[1]
           transcript[transcript_id]["parent"] = l[8].split(";")[1].split("gene-")[1]
           transcript[transcript_id]["ID"].append(transcript_id)

with open(accession_path, "r") as annot_inf:
    for line in annot_inf:
        l = line.strip().split(" ")
        if l[0] not in transcript:
            print("transcript not in gff")
            print(l)
            continue

        transcript[l[0]]["acc"] = l[1]
        swissprot = open(annotation_path + l[0] + ".info", "r")
        swissprot_inf = swissprot.readlines()
        swissprot_inf_GN = [x.strip() for x in swissprot_inf if x.startswith("GN") and not x.startswith("GN   ORFNames")]
        if len(swissprot_inf_GN) >0 :
           # swissprot_inf_GN = ";".join([x.split("=")[1].split(" ")[0].strip(";") for x in swissprot_inf_GN if "=" in x])
            swissprot_inf_GN = [";".join([y.split("=")[1].strip(";") for y in x.split() if y.split("=")[0] in ["Name", "Synonyms"]]) for x in swissprot_inf_GN if "=" in x]
        swissprot_inf_short = [x.strip() for x in swissprot_inf if x.startswith("DE") and "Short" in x]
        if len(swissprot_inf_short) >0 :
            swissprot_inf_short = ";".join([x.split("=")[1].strip(";") for x in swissprot_inf_short])

        transcript[l[0]]["GN"] = "".join(swissprot_inf_GN) + ";" + "".join(swissprot_inf_short)
        transcript[l[0]]["GN"] = transcript[l[0]]["GN"].strip(";")
        transcript[l[0]]["GN"] = re.sub("{.*?}", "", transcript[l[0]]["GN"])

with open(GO_path, "r") as go:
    for line in go:
        l = line.strip().split("\t")
        if len(l) >1 :
            try:
                transcript[l[0]]["GO"] = l[1].split(";")[:-1]
            except:
                print("transcript not")
                print(l)
                continue

for tr in transcript:
    if transcript[tr]["parent"] not in gene:
        gene_name = [x for x in gene if gene[x].get_id() == transcript[tr]["parent"]]
        if len(gene_name) == 0:
            print(f"{tr} transcript parent not in gene list")
            print("could be caused by difference between gene ID and Name in gff")
            continue
        print(gene_name)
        gene[gene_name[0]].parse_annot(transcript[tr])
        continue
    gene[transcript[tr]["parent"]].parse_annot(transcript[tr])

with open(output_path, "w") as outf:
    outf.write("\n".join([str(gene[x]) for x in gene]))
