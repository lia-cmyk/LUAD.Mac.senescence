library(tibble)
library(rms)
library(foreign)
library(survival)

rt <- read.table("clinical_nomo.txt", header=T,sep="\t",check.names = F,stringsAsFactors = F)
types <- read.table("risk.txt", header=T, sep="\t", check.names=F, stringsAsFactors=F)
types <- types[,c(1, ncol(types)-1)]
rt <- merge(rt, types, how="inner", by="id")
clinical <- column_to_rownames(rt, "id")

clinical$age<-factor(clinical$age,labels=c("<=65",">65"))
clinical$gender<-factor(clinical$gender,labels=c("male","female"))
clinical$stage<-factor(clinical$stage,labels=c("stage i + stage ii","stage iii + stage iv"))
clinical$T<-factor(clinical$T,labels=c("T1+T2","T3+T4"))
clinical$M<-factor(clinical$M,labels=c("M0","M1"))
clinical$N<-factor(clinical$N,labels=c("N0","N1+N2+N3"))

#??????????Cox?ع?ģ??



#1-year
cox1 <- cph(Surv(futime,fustat) ~ age + gender + stage +  T + M + N + riskScore,surv=T,x=T, y=T,time.inc = 1*30*12*1,data=clinical) 

cal <- calibrate(cox1, cmethod="KM", method="boot", u=1*30*12*1, m=105, B=319)

#??У׼ͼ
pdf("calibrate1.pdf",12,8)
par(mar = c(10,5,3,2),cex = 1.0)
plot(cal,lwd=3,lty=2,errbar.col="black",xlim = c(0,1),ylim = c(0,1),xlab ="Nomogram-Predicted Probability of 1-Year Survival",ylab="Actual 1-Year Survival",col="blue")
lines(cal[,c('mean.predicted','KM')],type = 'b',lwd = 3,col ="black" ,pch = 16)
box(lwd = 1)
abline(0,1,lty = 3,lwd = 3,col = "black")
dev.off()

#3-year
cox2 <- cph(Surv(futime,fustat) ~ age + gender + stage +  T + M + N + riskScore,surv=T,x=T, y=T,time.inc = 1*30*12*3,data=clinical) 

cal <- calibrate(cox2, cmethod="KM", method="boot", u=1*30*12*3, m=105, B=319)

#??У׼ͼ
pdf("calibrate3.pdf",12,8)
par(mar = c(10,5,3,2),cex = 1.0)
plot(cal,lwd=3,lty=2,errbar.col="black",xlim = c(0,1),ylim = c(0,1),xlab ="Nomogram-Predicted Probability of 3-Year Survival",ylab="Actual 3-Year Survival",col="blue")
lines(cal[,c('mean.predicted','KM')],type = 'b',lwd = 3,col ="black" ,pch = 16)
box(lwd = 1)
abline(0,1,lty = 3,lwd = 3,col = "black")
dev.off()

#5-year
cox3 <- cph(Surv(futime,fustat) ~ age + gender + stage +  T + M + N + riskScore,surv=T,x=T, y=T,time.inc = 1*30*12*5,data=clinical) 

cal <- calibrate(cox3, cmethod="KM", method="boot", u=1*30*12*5, m=105, B=319)

#??У׼ͼ
pdf("calibrate5.pdf",12,8)
par(mar = c(10,5,3,2),cex = 1.0)
plot(cal,lwd=3,lty=2,errbar.col="black",xlim = c(0,1),ylim = c(0,1),xlab ="Nomogram-Predicted Probability of 5-Year Survival",ylab="Actual 5-Year Survival",col="blue")
lines(cal[,c('mean.predicted','KM')],type = 'b',lwd = 3,col ="black" ,pch = 16)
box(lwd = 1)
abline(0,1,lty = 3,lwd = 3,col = "black")
dev.off()
