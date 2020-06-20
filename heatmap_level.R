library('pheatmap')
library('gplots')
Args <- commandArgs(TRUE)
input1 <- Args[1]
outfile1 <- Args[2]
#outfile2 <- Args[3]

countData <- read.table(input1, header = T, sep="\t")
countData_ok <- as.data.frame(countData[,2:ncol(countData)])
rownames(countData_ok) <- countData[,1]
sample_list <- colnames(countData_ok)

bk <- c(seq(-3,-0.1,by=0.01),seq(0,3,by=0.01))
mycolor <- c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2))
# pheatmap
#mycolor = colorRampPalette(c("blue","white","red"))(256)
pdf(outfile1,width=3,height=5)
#mycolors <- colorRampPalette(c("skyblue","white","firebrick1"))(100)
pheatmap(countData_ok, cluster_row=F, cluster_cols=F, legend=TRUE,  show_rownames=T, show_colnames=T,color = mycolor,legend_breaks=seq(-3,3,1),breaks=bk, annotation_legend=T)
dev.off()

# heatmap2
#countData_ok <- as.matrix(countData_ok)
#sample_color <- c(rep("skyblue",6),rep("hotpink",8))

#pdf(outfile2,pointsize=10)
#cr <- colorRampPalette(c("blue","white","red"))
#heatmap.2(countData_ok, scale="row", ColSideColors = sample_color, key=T, keysize=0.7, col=cr(1000), na.rm = TRUE, density.info="none",trace="none", Colv = F, Rowv = T, margins=c(9,9))
#dev.off()

