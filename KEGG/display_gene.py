import re,os,sys

kegg_file = sys.argv[1]
gene_name_map = sys.argv[2]
diff_gene_file = sys.argv[3]
#out_dir = sys.argv[4]
out_file=sys.argv[4]

kegg_dict = {}
KEGG = open(kegg_file,"r")
KEGG.readline()
for line in KEGG:
	line = line.strip()
	items = re.split("\t",line)
	kegg_id = items[1]
	kegg_name = items[2]
	qvalue = items[7]
	gene_list = re.split("/",items[8])
	kegg_dict[kegg_id] = {"kegg_name":kegg_name,"qvalue":qvalue,"gene_list":gene_list}
KEGG.close()

gene_name_map_dict = {}
GENE_MAP = open(gene_name_map,"r")
GENE_MAP.readline()
for line in GENE_MAP:
	if re.search("^\s$",line):
		continue
	line = line.strip()
	items = re.split("\t",line)
	gene_id= items[1]
	gene_name = items[0]
	gene_name_map_dict[gene_id] = gene_name
GENE_MAP.close()

diff_gene = {}
DIFF_GENE = open(diff_gene_file,"r")
DIFF_GENE.readline()
for line in DIFF_GENE:
	line = line.strip()
	items = re.split("\t",line)
	gene_name = items[1]
	logFC = items[4]
	FDR = items[8]
	diff_gene[gene_name] = {"logFC":logFC,"FDR":FDR}
DIFF_GENE.close()

OUT_FILE = open(out_file,"w")
OUT_FILE.write("KEGG_ID\tKEGG_Pathway\tqvalue(pathway_enrichment)\tgene_name\tlogFC(diff_exp_analysis)\tFDR(diff_exp_analysis)\n")
for kegg_id in kegg_dict:
	qvalue = kegg_dict[kegg_id]["qvalue"]
	kegg_name = kegg_dict[kegg_id]["kegg_name"]
	gene_list = kegg_dict[kegg_id]["gene_list"]
	if float(qvalue) <= 0.05:
		for gene_id in gene_list:
			gene_name = gene_name_map_dict[gene_id]
			logFC = diff_gene[gene_name]["logFC"]
			FDR = diff_gene[gene_name]["FDR"]
			content = [kegg_id,kegg_name,qvalue,gene_name,logFC,FDR]
			OUT_FILE.write("\t".join(content)+"\n")
OUT_FILE.close()
