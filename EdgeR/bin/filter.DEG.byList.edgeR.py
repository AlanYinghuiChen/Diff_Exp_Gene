#!/usr/bin/python
import re,sys,os
import pandas as pd
import numpy as np

diffExp_file=sys.argv[1]
Exp_file=sys.argv[2]
geneName=sys.argv[3]
outFile_all = sys.argv[4]
outFile_filtered_up = sys.argv[5]
outFile_filtered_down = sys.argv[6]
mapStat_file = sys.argv[7]

rpkm=pd.read_csv(Exp_file,sep="\t",index_col=0,header=0)

OUTFILE_ALL = open(outFile_all,"w")
OUTFILE_FILTERED_UP = open(outFile_filtered_up,"w")
OUTFILE_FILTERED_DOWN = open(outFile_filtered_down,"w")
MAPSTAT = open(mapStat_file,"w")

#read gene_ID gene_Name mapping file
ID2NAME=open(geneName,"r")
geneID_name_map = {}
for line in ID2NAME:
	line=line.strip()
	items=re.split("\t",line)
	geneID=items[4]
	geneName = items[3]
	strand = items[5]
	geneType = items[6]
	geneID_name_map[geneID] = {"geneName":geneName,"geneType":geneType,"strand":strand,"mapped":0}
ID2NAME.close()

#write header
DIFF_EXP=open(diffExp_file,"r")
header=DIFF_EXP.readline().strip()
header="\t".join(re.split("\t",header)[1:])
header="gene_ID\tgene_name\tgene_type\tstrand\t"+header+"\t"+"\t".join(list(rpkm.columns))+"\n"
OUTFILE_ALL.write(header)
OUTFILE_FILTERED_UP.write(header)
OUTFILE_FILTERED_DOWN.write(header)

# 1. read the diff_exp file 
# 2. add gene name , gene type and rpkm info.
for line in DIFF_EXP:
	line=line.strip()
	items=re.split("\t",line)
	gene = items[0]
	logFC = items[1]
	pValue = items[4]
	FDR = items[5]
	if gene in geneID_name_map:
		rpkm_list=rpkm.loc[gene,]
		rpkm_list_str =rpkm_list.astype('str')
		geneName = geneID_name_map[gene]["geneName"]
		geneType = geneID_name_map[gene]["geneType"]
		strand  = geneID_name_map[gene]["strand"]
		content = gene+"\t"+ geneName+"\t"+geneType+"\t"+strand+"\t"+"\t".join(items[1:])+"\t"+"\t".join(rpkm_list_str)
		OUTFILE_ALL.write(content+"\n")
		geneID_name_map[gene]["mapped"]+=1
		if logFC == "NA" or FDR == "NA":
			continue
		elif float(logFC) > 0 and float(FDR) <= 0.05:
			OUTFILE_FILTERED_UP.write(content+"\n")
		elif float(logFC) < 0 and float(FDR) <= 0.05:
			OUTFILE_FILTERED_DOWN.write(content+"\n")

for gene in geneID_name_map:
	if geneID_name_map[gene]["mapped"] != 1:
		MAPSTAT.write(gene+"\t"+str(geneID_name_map[gene]["mapped"])+"\n")

DIFF_EXP.close()
OUTFILE_ALL.close()
OUTFILE_FILTERED_UP.close()
OUTFILE_FILTERED_DOWN.close()
MAPSTAT.close()
