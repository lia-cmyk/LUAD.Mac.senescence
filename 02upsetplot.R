library(UpSetR)

x <- read.table("Mac.senescence.txt",sep="\t")
y <- read.table("DEmRNA.txt",sep="\t")

data <- intersect(x$V1, y$V1)
write.table(data,"inner.txt",sep="\t",row.names=F,col.names=F,quote=F)

a <- list(Mac.senescence=x$V1, DEmRNA=y$V1)


pdf("upset.pdf", onefile=F)
# upset(fromList(a), sets.bar.color=c("red", "blue", "black", "yellow", "pink"))
upset(fromList(a), sets.bar.color=c("red","blue"))
dev.off()