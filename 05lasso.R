
library("glmnet")
library("survival")

rt=read.table("multiInput.txt",header=T,sep="\t",row.names=1,check.names=FALSE)
rt1=log2(rt[,3:ncol(rt)]+0.01)
rt=cbind(rt[,1:2],rt1)
rt$futime <- ifelse(rt$futime==0,1,rt$futime)
rt[,"futime"]=rt[,"futime"]/365

# gene=read.table("gene.txt",header=F)
# rt=rt[,c("futime","fustat",as.vector(gene[,1]))]

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))

fit <- glmnet(x, y, family = "cox", maxit =5000)
pdf("lambda.pdf")
#tiff("lambda.tiff",width =18,height =18,units ="cm",compression="lzw",bg="white",res=400)
plot(fit, xvar = "lambda", label = TRUE)
abline(v=c(fit$lambda.min,fit$lambda.1se),lty="dashed")
dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit =5000)
pdf("cvfit.pdf")
#tiff("cvfit.tiff",width =18,height =18,units ="cm",compression="lzw",bg="white",res=400)
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef.txt",sep="\t",quote=F,row.names=F)
