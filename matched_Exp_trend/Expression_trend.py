import os,re,sys
import pandas as pd
import numpy as np

geneListFile = sys.argv[1]
inputDir = sys.argv[2]
outDir = sys.argv[3]

gene_list_dict = {}
GENE_LIST = open(geneListFile,"r")
for line in GENE_LIST:
	line = line.strip()
	items = re.split("\t",line)
	gene_id = items[4]
	gene_name = items[3]
	gene_list_dict[gene_id] = gene_name
GENE_LIST.close()

expression_data = {}
for root,dirs,files in os.walk(inputDir):
	for file_name in files:
		if re.search(".edgeR.TMM.countPerMillion.xls",file_name):
			sample_name = re.search("(.*).edgeR.TMM.countPerMillion.xls",file_name).group(1)
			file_path = os.path.join(root,file_name)
			for gene_id in gene_list_dict:
				if gene_id not in expression_data:
					expression_data[gene_id] = {}
				expression_data[gene_id][sample_name] = expression_data[gene_id].get(sample_name,{"Control":0,"Cancer":0})
				cpm_df = pd.read_csv(file_path,sep="\t",header=0,index_col=0)
				expression_data[gene_id][sample_name]["Control"] = cpm_df.loc[gene_id,:][0]
				expression_data[gene_id][sample_name]["Cancer"] = cpm_df.loc[gene_id,:][1]


header_list=["Control","Cancer"]
for gene_id in expression_data:
	gene_express_trend_file = os.path.join(outDir,gene_list_dict[gene_id]+".log2cpm.xls")
	up_cancer = os.path.join(outDir,gene_list_dict[gene_id]+"_up_in_cancer.xls")
	down_cancer = os.path.join(outDir,gene_list_dict[gene_id]+"_down_in_cancer.xls")
	gene_express_trend_dict = expression_data[gene_id]
	EPX_FILE = open(gene_express_trend_file,"w")
	UP_CANCER = open(up_cancer,"w")
	DOWN_CANCER = open(down_cancer,"w")
	EPX_FILE.write("Sample_Name\tSample_Type\tlog2CPM\n")
	for sample_name in gene_express_trend_dict:
		gene_name = gene_list_dict[gene_id]
		control_exp = float(gene_express_trend_dict[sample_name]["Control"])
		cancer_exp = float(gene_express_trend_dict[sample_name]["Cancer"])
		if cancer_exp > control_exp:
			UP_CANCER.write(sample_name+"\n")
		elif cancer_exp < control_exp:
			DOWN_CANCER.write(sample_name+"\n")

		EPX_FILE.write(sample_name+"\tControl\t"+str(control_exp)+"\n")
		EPX_FILE.write(sample_name+"\tCancer\t"+str(cancer_exp)+"\n")
	EPX_FILE.close()
	UP_CANCER.close()
	DOWN_CANCER.close()
#	gene_express_trend_df = pd.DataFrame.from_dict(gene_express_trend_dict, orient='index', columns=header_list)
#	gene_express_trend_df.to_csv(gene_express_trend_file,sep="\t",index_label="Sample_Name",header=True,index=True,columns=header_list)

	#write the R script to plot log2CPM
	plot_Rscript = os.path.join(outDir,gene_name+".log2cpm.R")
	PLOT_R = open(plot_Rscript,"w")
	PLOT_R.write("""library("ggplot2")
cpm_df = read.csv('"""+gene_name+""".log2cpm.xls',sep='\\t',header=T)

control <- cpm_df$log2CPM[cpm_df$Sample_Type == "Control"]
cancer <- cpm_df$log2CPM[cpm_df$Sample_Type == "Cancer"]
wilcox.test(control,cancer,paired = TRUE)

pdf('"""+gene_name+""".log2cpm.pdf')

ggplot(data=cpm_df,mapping=aes(x=Sample_Type,y=log2CPM,group=Sample_Name))+
  geom_line()+ geom_point() + labs(title='"""+gene_name+"""_Tissue')+
  theme(plot.title = element_text(hjust = 0.5, size=14,face="bold"), axis.title.x=element_text(size=12,face="bold"), axis.title.y=element_text(size=12,face="bold"),axis.text.x=element_text(size=10,face="bold"),axis.text.y=element_text(size=10,face="bold"))

dev.off()
""")
	PLOT_R.close()
