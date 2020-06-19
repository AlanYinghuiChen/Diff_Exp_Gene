import os,re,sys
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
input = sys.argv[1]
outDir = sys.argv[2]
log = os.path.join(outDir,"gene_sample_filteration.log")
LOG=open(log,"w")

data= pd.read_csv(input,sep="\t",header=0,index_col=0)

data_judge_gene = data.apply(lambda x: x >= 1)
data_judge_gene_sum = data_judge_gene.apply(lambda x: x.sum(),axis=0)
data_gene_num_filtered = data.loc[:,data_judge_gene_sum >= 10000]

#boxplot for gene number
plt.boxplot(list(data_judge_gene_sum))
plt.savefig(os.path.join(outDir,"gene_num.boxplot.png"))
plt.close()

LOG.write("Filtered sample (less than 10000 genes detected, >= 1 spanning read)\n")
LOG.write(str(data_judge_gene_sum[data_judge_gene_sum < 10000]))

total_sample_num = data_gene_num_filtered.shape[1]
num_threshold = int(total_sample_num * 0.1)

data_judge = data_gene_num_filtered.apply(lambda x: x >= 30)
data_judge_sum = data_judge.apply(lambda x: x.sum(),axis=1 )

data_judge_sum_filtered = data_judge_sum.apply(lambda x: x> num_threshold)

data_filtered = data_gene_num_filtered[data_judge_sum_filtered == True]
LOG.write("\ngene_num_after_filteration: %d \n"%data_filtered.shape[0])

data_filtered.to_csv(os.path.join(outDir,"all.intron_spanning.read_count.HTSeq.union.filtered.xls"),sep="\t",header=True,index=True)
