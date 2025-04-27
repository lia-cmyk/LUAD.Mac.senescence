library(ggpubr)
library(ggforce)

rt <- read.table("risk.txt",header = T,sep = "\t",check.names = F)
clsdata <- rt[,c("id","risk")]
TCIA <- read.table("TCIA-ClinicalData.txt", header=T, check.names=F, sep="\t")

data <- merge(clsdata, TCIA, how="inner", by="id")
data$risk <- factor(data$risk,levels=c("low","high"))

# compare_means(expr ~ group, data = Gene_exp)
pdf("ips_ctla4_pos_pd1_pos.pdf")
ggplot(data, aes(x=risk, y=ips_ctla4_pos_pd1_pos)) + geom_violin(aes(fill=risk), trim=FALSE) + geom_boxplot(width=0.2,position=position_dodge(0.9)) + scale_fill_manual(values=c("#5555FF", "#FF7744")) + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange", color = "red") + stat_compare_means(method="t.test", label="p.signif")  + theme_bw()
dev.off()
