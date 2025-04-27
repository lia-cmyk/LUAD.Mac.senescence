
rt=read.table("risk.txt",header=T,sep="\t",check.names=F,row.names=1)
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>10]=10

pdf(file="riskScore.pdf",width = 15, height = 6)
par(mar=c(5,5,4,4))
#tiff(file="riskScore.tiff"width = 20, height = 12,units ="cm",compression="lzw",bg="white",res=600)
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("green",lowLength),
     rep("red",highLength)), cex.lab=1.8)
abline(h=median(rt$riskScore),v=lowLength,lty=2)
dev.off()
