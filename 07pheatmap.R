
rt=read.table("risk.txt",sep="\t",header=T,row.names=1,check.names=F)
rt=rt[order(rt$riskScore),]
rt1=rt[c(3:(ncol(rt)-2))]
rt1=t(rt1)

#rt1=log2(rt1+0.001)
library(pheatmap)
annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)

pdf(file="heatmap.pdf",width = 6,height = 6)
#tiff(file="heatmap.tiff",width = 35,height = 10,units ="cm",compression="lzw",bg="white",res=600)
pheatmap(rt1, 
         annotation=annotation, 
         cluster_cols = FALSE,
         fontsize_row=11,
         fontsize_col=5,
         scale="row",
         color = colorRampPalette(c("limegreen", "khaki1", "red3"))(50) )
dev.off()

