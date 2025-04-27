
library("limma")
normalCount=59     
tumorCount=535        
          
rt=read.table("mRNA_symbol.txt",sep="\t",header=T,check.names=F)          
rt=as.matrix(rt)
mRNA=rt[,1]
mRNA <- gsub("(.*?)\\|.*","\\1",mRNA)
rownames(rt)=mRNA
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

group=c(rep("normal",normalCount),rep("tumor",tumorCount))
design <- model.matrix(~factor(group))
colnames(design)=levels(factor(group))
rownames(design)=colnames(data)
v <-voom(data, design = design, plot = F, save.plot = F)
out=v$E
out=rbind(ID=colnames(out),out)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)

