
rt=read.table("risk.txt",header=T,sep="\t",check.names=F,row.names=1)
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
color=as.vector(rt$fustat)
color[color==1]="red"
color[color==0]="green"

pdf(file="survStat.pdf",width = 15, height = 6)
par(mar=c(5,5,4,4))
#tiff(file="survStat.tiff",width = 20, height = 12,units ="cm",compression="lzw",bg="white",res=600)
plot(rt$futime,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color, cex.lab=1.8)
legend("topright", 
       c("Dead", "Alive"),
       pch=19,
       col=c("red","green"))
abline(v=lowLength,lty=2)
dev.off()

