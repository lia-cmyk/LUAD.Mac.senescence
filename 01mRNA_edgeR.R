######Video source: https://shop119322454.taobao.com
#source("http://bioconductor.org/biocLite.R")   #source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#install.packages("gplots")

foldChange=1
padj=0.05
                   #???ù???Ŀ¼
library("edgeR")
rt=read.table("mRNA_symbol.txt",sep="\t",header=T,check.names=F)  #?ĳ??Լ????ļ???
rt=as.matrix(rt)
mRNA=rt[,1]
mRNA <- gsub("(.*?)\\|.*","\\1",mRNA)
rownames(rt)=mRNA
exp=rt[,2:ncol(rt)]

dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
# data=data[rowMeans(data)>1,]

group=c(rep("normal",59),rep("tumor",535))                         #???հ?֢????????Ʒ??Ŀ?޸?
design <- model.matrix(~group)
y <- DGEList(counts=data,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y,pair = c("normal","tumor"))
topTags(et)
ordered_tags <- topTags(et, n=100000)

allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff
newData=y$pseudo.counts

write.table(diff,file="edgerOut.xls",sep="\t",quote=F)
diffSig = diff[(diff$FDR < padj & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]
write.table(diffSig, file="diffSig.xls",sep="\t",quote=F)
diffUp = diff[(diff$FDR < padj & (diff$logFC>foldChange)),]
write.table(diffUp, file="up.xls",sep="\t",quote=F)
diffDown = diff[(diff$FDR < padj & (diff$logFC<(-foldChange))),]
write.table(diffDown, file="down.xls",sep="\t",quote=F)

normalizeExp=rbind(id=colnames(newData),newData)
write.table(normalizeExp,file="normalizeExp.txt",sep="\t",quote=F,col.names=F)   #???????л???У?????ı???ֵ??normalizeExp.txt??
diffExp=rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file="diffmRNAExp.txt",sep="\t",quote=F,col.names=F)         #????????????У?????ı???ֵ??diffmRNAExp.txt??

heatmapData <- newData[rownames(diffSig),]

#volcano
pdf(file="vol.pdf")
#tiff(file="vol.tiff",width =12,height =12,units ="cm",compression="lzw",bg="white",res=400)
xMax=300	#????PValue=0??????
#xMax=max(-log10(allDiff$FDR))+1
yMax=12
plot(-log10(allDiff$FDR), allDiff$logFC, xlab="-log10(FDR)",ylab="logFC",
     main="mRNA_Volcano", xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.4, cex.main=1.5, cex.lab=1.3)
diffSub=allDiff[allDiff$FDR<padj & allDiff$logFC>foldChange,]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="red",cex=0.4)
diffSub=allDiff[allDiff$FDR<padj & allDiff$logFC<(-foldChange),]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="green",cex=0.4)
abline(h=0,lty=2,lwd=3)
dev.off()

#heatmaps
hmExp=log10(heatmapData+0.001)
library('gplots')
hmMat=as.matrix(hmExp)
#pdf(file="heatmap.pdf",width=60,height=90)
tiff(file="heatmap.tiff",width =60,height =90,units ="cm",compression="lzw",bg="white",res=400)
par(oma=c(10,3,3,7))
heatmap.2(hmMat,col='greenred',trace="none")
dev.off()
