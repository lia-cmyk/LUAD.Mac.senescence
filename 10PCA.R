
library("FactoMineR")#画主成分分析图需要加载这两个包
library("factoextra") 


# The variable group_list (index = 54676) is removed
dat <- read.table("risk.txt", sep="\t", header=T, check.names=F, row.names=1)
dat <- dat[,-c(1,2,(ncol(dat)-1))]
dat <- dat[order(dat[,ncol(dat)], decreasing=T),]
dat$risk <- ifelse(dat$risk=="low", "low_risk", "high_risk")

# before PCA analysis
pdf("pca.pdf")
dat.pca <- PCA(dat[,-ncol(dat)],graph = FALSE)#现在dat最后一
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = dat$risk, # color by groups
            addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
dev.off()

