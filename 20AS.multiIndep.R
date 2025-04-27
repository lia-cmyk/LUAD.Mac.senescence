
library(tibble)
library(survival)
library(forestplot)

clrs <- fpColors(box="red",line="darkblue", summary="royalblue")  

rt=read.table("clinical_cox.txt",header=T,sep="\t",check.names=F)
risk <- read.table("risk.txt", header=T, sep="\t", check.names=F)
risk <- risk[,c(1,ncol(risk)-1)]
rt <- merge(rt, risk, by='id', all=T, how='inner')
rt <- column_to_rownames(rt, 'id')


multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
             HR=multiCoxSum$conf.int[,"exp(coef)"],
             HR.95L=multiCoxSum$conf.int[,"lower .95"],
             HR.95H=multiCoxSum$conf.int[,"upper .95"],
             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

#????ɭ??ͼ
rt=read.table("multiCox.xls",header=T,sep="\t",row.names=1,check.names=F)
# rownames(rt) <- gsub("`", "", rownames(rt))
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          #????ͼƬ????

pdf("forest.pdf",width=8,height=4,onefile=F)

forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.3,
           xlab="Hazard ratio"
           )
dev.off()
