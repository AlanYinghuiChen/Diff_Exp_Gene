import os,re,sys

input_file = sys.argv[1]
databas_file = sys.argv[2]
geneList_file = sys.argv[3]
out_gene_list = sys.argv[4]
out_unmap_file = sys.argv[5]

geneList = []
GENELIST = open(geneList_file,"r")
for line in GENELIST:
	line = line.strip()
	geneList.append(line)
GENELIST.close()

content_dict = {}
DATABASE = open(databas_file,"r")
header = DATABASE.readline()
for line in DATABASE:
	line = line.strip()
	items = re.split("\t",line)
	geneID = items[0]
	geneSymbol = items[1]
	if geneID not in content_dict:
		content_dict[geneID]= geneSymbol
	else:
		content_dict[geneID]=";".join([content_dict[geneID],geneSymbol])
DATABASE.close()

INPUT = open(input_file,"r")
OUT_GENELIST = open(out_gene_list,"w")
OUT_UNMAP = open(out_unmap_file,"w")
input_header = INPUT.readline()
OUT_GENELIST.write("Gene_ID\tGene_Symbol\t"+input_header)
for line in INPUT:
	line = line.strip()
	items = re.split("\t",line)
	gene_id = items[0]

	if gene_id in content_dict:
		geneSymbol = content_dict[gene_id]
		if geneSymbol in geneList:
			items.insert(1,geneSymbol)
			OUT_GENELIST.write("\t".join(items)+"\n")
	else:
		OUT_UNMAP.write(gene_id+"\n")
OUT_GENELIST.close()
OUT_UNMAP.close()
