
library("limma")
 

rt <- read.table("mRNA_symbol.txt",sep="\t",header=T,check.names=F)          
rt <- as.matrix(rt)
mRNA <- rt[,1]
mRNA <- gsub("(.*?)\\|.*","\\1",mRNA)
rownames(rt) <- mRNA
exp <- rt[,2:ncol(rt)]
dimnames <- list(rownames(exp),colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data <- avereps(data)
data <- data[rowMeans(data)>0,]

group <- sapply(strsplit(colnames(data), "\\-"),"[",4)
group <- sapply(strsplit(group,""),"[",1)
group <- gsub("2","1",group)
data <- data[,group==0]
data <- t(data[rowMeans(data)>0,])

rownames(data) <- gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data <- avereps(data)
data <- as.data.frame(t(data))
normalize <- t(apply(data, 1, function(x)x-(mean(x))))
normalize <- cbind(id=rownames(normalize), normalize)
write.table(normalize, "TIDE.input.txt", sep="\t", quote=F, row.names=F)


TIDE.Score <- read.table("TIDE.output.csv", sep=",", header=T, check.names=F)
TIDE.Score <- TIDE.Score[,c("Patient","TIDE")]
colnames(TIDE.Score)[1] <- "id"

rt <- read.table("risk.txt",header = T,sep = "\t",check.names = F)
clsdata <- rt[,c("id","risk")]

data <- merge(clsdata, TIDE.Score, how="inner", by="id")
data$risk <- factor(data$risk,levels=c("low","high"))

pdf("TIDE.Score.pdf")
ggplot(data, aes(x=risk, y=TIDE)) + geom_violin(aes(fill=risk), trim=FALSE) + geom_boxplot(width=0.2,position=position_dodge(0.9)) + scale_fill_manual(values=c("#5555FF", "#FF7744")) + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange", color = "red") + stat_compare_means(method="t.test", label="p.signif")  + theme_bw()
dev.off()