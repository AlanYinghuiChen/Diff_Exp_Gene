#!/usr/bin/python
import re,sys,os
import pandas as pd
import numpy as np

sig_type = sys.argv[1]
diffExp_file=sys.argv[2]
Exp_file=sys.argv[3]
geneName=sys.argv[4]
outFile_all = sys.argv[5]
outFile_filtered_up = sys.argv[6]
outFile_filtered_down = sys.argv[7]
mapStat_file = sys.argv[8]

rpkm=pd.read_csv(Exp_file,sep="\t",index_col=0,header=0)

OUTFILE_ALL = open(outFile_all,"w")
OUTFILE_FILTERED_UP = open(outFile_filtered_up,"w")
OUTFILE_FILTERED_DOWN  = open(outFile_filtered_down,"w")
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
diff_exp_header=DIFF_EXP.readline().strip()
header="gene_ID\tgene_name\tgene_type\tstrand\t"+diff_exp_header+"\t"+"\t".join(list(rpkm.columns))+"\n"
OUTFILE_ALL.write(header)
OUTFILE_FILTERED_UP.write(header)
OUTFILE_FILTERED_DOWN.write(header)

# 1. read the diff_exp file 
# 2. add gene name , gene type and rpkm info.
for line in DIFF_EXP:
	line=line.strip()
	items=re.split("\t",line)
	gene = items[0]
	logFC = items[4]
	pValue = items[7]
	FDR = items[8]
	if gene in geneID_name_map:
		rpkm_list=rpkm.loc[gene,]
		rpkm_list_str =rpkm_list.astype('str')
		geneName = geneID_name_map[gene]["geneName"]
		geneType = geneID_name_map[gene]["geneType"]
		strand  = geneID_name_map[gene]["strand"]
		content = gene+"\t"+ geneName+"\t"+geneType+"\t"+strand+"\t"+"\t".join(items[1:])+"\t"+"\t".join(rpkm_list_str)
		OUTFILE_ALL.write(content+"\n")
		geneID_name_map[gene]["mapped"]+=1
		if logFC == "NA" or FDR == "NA" or pValue == "NA":
			continue
		if re.search("FDR",sig_type,re.I):
			sig_filter = float(FDR)
		else:
			sig_filter = float(pValue)

		if float(logFC) > 0 and sig_filter <= 0.01:
			OUTFILE_FILTERED_UP.write(content+"\n")
		if float(logFC) < 0 and sig_filter <= 0.01:
			OUTFILE_FILTERED_DOWN.write(content+"\n")

for gene in geneID_name_map:
	if geneID_name_map[gene]["mapped"] != 1:
		MAPSTAT.write(gene+"\t"+str(geneID_name_map[gene]["mapped"])+"\n")

DIFF_EXP.close()
OUTFILE_ALL.close()
OUTFILE_FILTERED_UP.close()
OUTFILE_FILTERED_DOWN.close()
MAPSTAT.close()
