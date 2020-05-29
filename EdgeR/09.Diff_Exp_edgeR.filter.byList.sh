echo start `date`
cd /BioII/lulab_b/chenyinghui/project/exRNA/CRC/cfRNA/09.Diff_Exp_EdgeR

#python ../bin/filter.DEG.byList.edgeR.FDR.py 0 0.01 NC_PKU-VS-CRC.edgeR.diffExp.xls NC_PKU-VS-CRC.edgeR.countPerMillion.xls /BioII/lulab_b/chenyinghui/database/Homo_sapiens/GRCh38/gencode.v32.annotation.gene_info.bed NC_PKU-VS-CRC.edgeR.annot.xls NC_PKU-VS-CRC.edgeR.up.FDR.xls NC_PKU-VS-CRC.edgeR.down.FDR.xls NC_PKU-VS-CRC.edgeR.FDR.mapStat.xls

#sed -n '2,$p' NC_PKU-VS-CRC.edgeR.up.FDR.xls | awk -F "\t" '{print $2}' - > NC_PKU-VS-CRC.edgeR.up.FDR.list
#sed -n '2,$p' NC_PKU-VS-CRC.edgeR.down.FDR.xls | awk -F "\t" '{print $2}' - > NC_PKU-VS-CRC.edgeR.down.FDR.list

cd /BioII/lulab_b/chenyinghui/project/exRNA/CRC/cfRNA/09.Diff_Exp_EdgeR_1.5_0.01

# FDR = 1.5, log2(foldchange) = 0.585, FR=0.01
python ../bin/filter.DEG.byList.edgeR.FDR.py 0.585 0.01 NC_PKU-VS-CRC.edgeR.diffExp.xls NC_PKU-VS-CRC.edgeR.countPerMillion.xls /BioII/lulab_b/chenyinghui/database/Homo_sapiens/GRCh38/gencode.v32.annotation.gene_info.bed NC_PKU-VS-CRC.edgeR.annot.xls NC_PKU-VS-CRC.edgeR.up.FDR.xls NC_PKU-VS-CRC.edgeR.down.FDR.xls NC_PKU-VS-CRC.edgeR.FDR.mapStat.xls

sed -n '2,$p' NC_PKU-VS-CRC.edgeR.up.FDR.xls | awk -F "\t" '{print $2}' - > NC_PKU-VS-CRC.edgeR.up.FDR.list
sed -n '2,$p' NC_PKU-VS-CRC.edgeR.down.FDR.xls | awk -F "\t" '{print $2}' - > NC_PKU-VS-CRC.edgeR.down.FDR.list
echo end `date`
