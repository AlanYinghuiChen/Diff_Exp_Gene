library(edgeR) 
args <- commandArgs(TRUE)

sample_name <- args[1]
inputDir <- args[2]
outDir <- args[3]

input <- paste(inputDir,"/",sample_name,".featurecounts.txt",sep="")
outfile <- paste(outDir,"/",sample_name,".tissue.edgeR.diffExp.xls",sep="")
cpmFile <- paste(outDir,"/",sample_name,".tissue.edgeR.countPerMillion.xls",sep="")

#input <- "/BioII/lulab_b/chenyinghui/project/exRNA/result_yunfan/CRC/tissue/02.featurecount_merge/all.featurecounts.txt"
#outfile <- "/BioII/lulab_b/chenyinghui/project/exRNA/result_yunfan/CRC/tissue/03.Diff_Exp_edgeR/NC-VS-CRC.tissue.edgeR.diffExp.xls"
#cpmFile <- "/BioII/lulab_b/chenyinghui/project/exRNA/result_yunfan/CRC/tissue/03.Diff_Exp_edgeR/NC-VS-CRC.tissue.edgeR.countPerMillion.xls"

raw <- read.table(input,header = T, sep="\t")
countData <-data.frame()
countData <- raw[,2:(ncol(raw))]
row.names(countData) <- raw[,1]

dgListGroups <- c(rep("Control",1),rep("CRC",1))
dgList <- DGEList(counts=countData, genes=rownames(countData),group=factor(dgListGroups))

countsPerMillion <- cpm(dgList, normalized.lib.sizes=TRUE, log=TRUE)
write.table(countsPerMillion,file=cpmFile,quote=F, sep="\t", row.names=T, col.names=T)

#countCheck <- countsPerMillion > 1
#keep <- which(rowSums(countCheck) >= 5)
#dgList<-dgList[keep,]

dgList <- calcNormFactors(dgList, method="TMM")

design.mat <- model.matrix(~0 + dgList$sample$group)
colnames(design.mat) <- levels(dgList$sample$group)

d2 <- estimateGLMCommonDisp(dgList, method="deviance", robust=TRUE, subset=NULL)
#如果对照或实验组只有1个样本，那么下面2种分布函数不可用
#d2 <- estimateGLMTrendedDisp(d2, design=design.mat)
#d2 <- estimateGLMTagwiseDisp(d2, design=design.mat)

fit <- glmFit(d2, design.mat)
lrt <- glmLRT(fit,contrast=c(-1,1))

edgeR_result <- topTags(lrt,n=nrow(dgList))
write.table(edgeR_result$table,file = outfile,sep = "\t",quote = F,row.names=F,col.names=T)

