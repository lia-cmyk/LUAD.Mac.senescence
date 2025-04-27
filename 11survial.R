
library(survival)
library(survminer)
rt=read.table("risk.txt",header=T,sep="\t")
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
#pValue=round(pValue,3)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
summary(fit)    #???????????

pdf(file="survival_TCGA.pdf", width=5, height=5, onefile=F)


ggsurvplot(fit, data =rt,
           pval =T, # 在图上添加log rank检验的p值
           conf.int =T,# 添加置信区间
           risk.table=T, # 在图下方添加风险表
           legend.labs=c("high risk","low risk"), 
           legend.title="strata",#表头标签
           title="Survival curve",#改一下整体名称
           ylab="Survival rate",
           xlab = " Time (Years)",#修改X轴Y轴名称
           risk.table.col = "strata", # 风险表加颜色
           linetype = "strata", # 生存曲线的线型
           surv.median.line = "hv", # 标注出中位生存时间
           ggtheme = theme_bw(), #背景布局
           palette = c("red", "blue")) # 图形颜色风格
dev.off()
