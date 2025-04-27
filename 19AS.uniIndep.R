
library(survival)
library(forestplot)
                           #???ù???Ŀ¼

clrs <- fpColors(box="green",line="darkblue", summary="royalblue")             #????ɭ??ͼ??ɫ

rt=read.table("clinical_cox.txt",header=T,sep="\t",check.names=F)
risk <- read.table("risk.txt", header=T, sep="\t", check.names=F)
risk <- risk[,c(1,ncol(risk)-1)]
rt <- merge(rt, risk, by='id', how='inner')

outTab=data.frame()
for(i in colnames(rt[,4:ncol(rt)])){
	 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	 coxSummary = summary(cox)
	 coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	 outTab=rbind(outTab,
	              cbind(id=i,
	              HR=coxSummary$conf.int[,"exp(coef)"],
	              HR.95L=coxSummary$conf.int[,"lower .95"],
	              HR.95H=coxSummary$conf.int[,"upper .95"],
	              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
	              )
}
write.table(outTab,file="uniCox.xls",sep="\t",row.names=F,quote=F)

#????ɭ??ͼ
rt=read.table("uniCox.xls",header=T,sep="\t",check.names=F,row.names=1)
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

