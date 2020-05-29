library("ggplot2")
cpm_df = read.csv('IRF7.log2cpm.xls',sep='\t',header=T)

control <- cpm_df$log2CPM[cpm_df$Sample_Type == "Control"]
cancer <- cpm_df$log2CPM[cpm_df$Sample_Type == "Cancer"]
wilcox.test(control,cancer,paired = TRUE)

pdf('IRF7.log2cpm.pdf')

ggplot(data=cpm_df,mapping=aes(x=Sample_Type,y=log2CPM,group=Sample_Name))+
  geom_line()+ geom_point() + labs(title='IRF7_Tissue')+
  theme(plot.title = element_text(hjust = 0.5, size=14,face="bold"), axis.title.x=element_text(size=12,face="bold"), axis.title.y=element_text(size=12,face="bold"),axis.text.x=element_text(size=10,face="bold"),axis.text.y=element_text(size=10,face="bold"))

dev.off()
