library('pheatmap')
library('gplots')
Args <- commandArgs(TRUE)
input1 <- Args[1]
outfile1 <- Args[2]
outfile2 <- Args[3]

countData <- read.table(input1, header = T, sep="\t")
countData_ok <- as.data.frame(countData[,3:ncol(countData)])
rownames(countData_ok) <- countData[,2]
sample_list <- colnames(countData_ok)

#type <- c(rep("Control",54),rep("CRC",15))
type <- c(rep("CRC",75),rep("Control",50))
sample <- data.frame(type)
rownames(sample) <- c(sample_list)
colnames(sample) <- "condition"

# pheatmap
sample_color <- list(sample_type=c(Control="skyblue",CRC="hotpink"))

pdf(outfile1,pointsize=10)
mycolors <- colorRampPalette(c("navy","white","firebrick3"))(1000)
pheatmap(countData_ok, clustering_distance_rows="correlation", clustering_distance_cols="correlation", clustering_method="ward.D", cluster_row=T, cluster_cols=T, scale="row", legend=TRUE, display_numbers=F, show_rownames=T, show_colnames=F, annotation_col=sample, annotation_colors=sample_color, col=mycolors)
dev.off()

# heatmap2
countData_ok <- as.matrix(countData_ok)
sample_color <- c(rep("hotpink",75),rep("skyblue",50))

pdf(outfile2,pointsize=10)
cr <- colorRampPalette(c("blue","white","red"))
heatmap.2(countData_ok, scale="row", ColSideColors = sample_color, key=T, keysize=0.7, col=cr(1000), na.rm = TRUE, density.info="none",trace="none", Colv = T, Rowv = T, margins=c(9,9))
dev.off()

