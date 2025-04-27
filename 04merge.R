
rt <- read.table("tumor_clinical.txt",sep="\t",header=T,check.names=F,row.names=1)
rt1 <- rt[,c(1:2)]
rt2 <- rt[,3:ncol(rt)]

genes <- read.table("0.05genes.txt")
rt2 <- rt2[,which(colnames(rt2)%in%genes$V1)]
rt <- cbind(id=rownames(rt),rt1,rt2)
write.table(rt,"multiInput.txt",sep="\t",quote=F,row.names=F)