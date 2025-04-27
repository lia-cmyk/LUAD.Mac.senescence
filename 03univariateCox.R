outTab=data.frame()

library(survival)
rt=read.table("tumor_clinical.txt",header=T,sep="\t",row.names=1,check.names=F)
rt1=log2(rt[,3:ncol(rt)]+0.01)
rt=cbind(rt[,1:2],rt1)

for(i in colnames(rt[,3:ncol(rt)])){
 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
 coxSummary = summary(cox)
 outTab=rbind(outTab,cbind(gene=i,HR=coxSummary$coefficients[,"exp(coef)"],
 z=coxSummary$coefficients[,"z"],
 pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
}

write.table(outTab,file="univariateCox.xls",sep="\t",row.names=F,quote=F)


