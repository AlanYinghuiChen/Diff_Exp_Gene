library(edgeR) 
Args <- commandArgs(TRUE)
input <- "/BioII/lulab_b/chenyinghui/project/exRNA/CRC/cfRNA/08.read_count_HTSeq_rmDup_merge/all.intron_spanning.read_count.HTSeq.union.filtered.xls"
outfile <- "/BioII/lulab_b/chenyinghui/project/exRNA/CRC/cfRNA/09.Diff_Exp_EdgeR/NC_PKU-VS-CRC.edgeR.diffExp.xls"
cpmFile <- "/BioII/lulab_b/chenyinghui/project/exRNA/CRC/cfRNA/09.Diff_Exp_EdgeR/NC_PKU-VS-CRC.edgeR.countPerMillion.xls"

raw <- read.table(input,header = T, sep="\t")
countData <-data.frame()
countData <- raw[,2:(ncol(raw))]
row.names(countData) <- raw[,1]

dgListGroups <- c(rep("Control",50),rep("CRC",75))
dgList <- DGEList(counts=countData, genes=rownames(countData),group=factor(dgListGroups))

dgList <- calcNormFactors(dgList, method="TMM")
countsPerMillion <- cpm(dgList, normalized.lib.sizes=TRUE)
write.table(countsPerMillion,file=cpmFile,quote=F, sep="\t", row.names=T, col.names=T)

#countCheck <- countsPerMillion > 1
#keep <- which(rowSums(countCheck) >= 5)
#dgList<-dgList[keep,]
design.mat <- model.matrix(~0 + dgList$sample$group)
colnames(design.mat) <- levels(dgList$sample$group)

d2 <- estimateGLMCommonDisp(dgList, design=design.mat)
d2 <- estimateGLMTrendedDisp(d2, design=design.mat)
d2 <- estimateGLMTagwiseDisp(d2, design=design.mat)

fit <- glmFit(d2, design.mat)
lrt <- glmLRT(fit,contrast=c(-1,1))

edgeR_result <- topTags(lrt,n=nrow(dgList))
write.table(edgeR_result$table,file = outfile,sep = "\t",quote = F,row.names=F,col.names=T)

