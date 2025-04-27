
library(survival)
library(timeROC)


title <- "TCGA_ROC"
rt <- read.table("risk.txt",header=T,sep="\t",check.names=F,row.names=1)
ROC_rt <- timeROC(T=rt$futime,delta=rt$fustat,
                    marker=rt$riskScore,cause=1,
                    weighting='marginal',
                    times=c(1,3,5),ROC=TRUE)

pdf(file=paste0("ROC_TCGA",".pdf"), width=5, height=5)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
plot(ROC_rt,time=1,title=FALSE,col='red',lwd=2)

plot(ROC_rt,time=3,title=FALSE,lwd=2,col='green',add=TRUE)

plot(ROC_rt,time=5,title=FALSE,lwd=2,col='blue',add=TRUE)
title(main=title)
legend('bottomright',
        c(paste0('AUC at 1 years: ',round(ROC_rt$AUC[1],2)),
        paste0('AUC at 3 years: ',round(ROC_rt$AUC[2],2)),
        paste0('AUC at 5 years: ',round(ROC_rt$AUC[3],2))),
        col=c('red','green','blue'),lwd=2,bty='n')
dev.off()
