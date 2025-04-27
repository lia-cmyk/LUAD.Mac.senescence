HLA = read.table("genes.txt",sep="\t",header=F,check.names=F)
rt = read.table("normalizeExp.txt",sep="\t",header=T,row.names = 1,check.names=F)   #??ȡ?ļ?
group <- sapply(strsplit(colnames(rt),"\\-"),"[",4)
group <- sapply(strsplit(group,""),"[",1)
group <- gsub("2","1",group)
rt <- rt[,group==0]

rt_HLA = rt[rownames(rt)%in%HLA$V1, ]
rt_HLA <- cbind(id=rownames(rt_HLA),rt_HLA)
write.table(rt_HLA, file="HLAexp.txt", sep="\t", quote=F, col.names=T, row.names=F)
