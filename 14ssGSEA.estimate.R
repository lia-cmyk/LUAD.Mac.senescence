
library(limma)
library(estimate)


inputFile="raw_uniq.symbol.txt"                                                  #?????ļ?????


rt=read.table(inputFile,sep="\t",header=T,check.names=F,row.names=1)


group=sapply(strsplit(colnames(rt),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
rt=rt[,group==0]
rt=rt[rowMeans(rt)>0,]
out=rbind(ID=colnames(rt),rt)
#???????????ľ????ļ?
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)

#????estimate??
filterCommonGenes(input.f="uniq.symbol.txt", 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct")

#????ÿ????Ʒ?Ĵ???
scores=read.table("estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
out=rbind(ID=colnames(scores),scores)
write.table(out,file="scores.txt",sep="\t",quote=F,col.names=F)
