
inputFile="uniq.symbol.txt"                                         #?????ļ?
gmtFile="immune.gmt"                                           #GMT?ļ?

#???ð?
library(GSVA)
library(limma)
library(GSEABase)

rt <- as.matrix(read.table(inputFile,sep="\t",header=T,check.names=F,row.names=1))

group <- sapply(strsplit(colnames(rt),"\\-"),"[",4)
group <- sapply(strsplit(group,""),"[",1)
group <- gsub("2","1",group)
rt <- rt[,group==0]
mat <- rt[rowMeans(rt)>1,]

geneSet <- getGmt(gmtFile, 
               geneIdType=SymbolIdentifier())

#ssgsea????
ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE,min.sz=1.1)
#????ssGSEA score????????
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
#??ssGSEA score???н???
ssgseaOut=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
write.table(ssgseaOut,file="ssgseaOut.txt",sep="\t",quote=F,col.names=F)
