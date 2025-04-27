#install.packages("rms")


library(survival)
library(rms)
library(tibble)

inputfile="clinical_nomo.txt"

rt <- read.table(inputfile,header=T,sep="\t",check.names = F,stringsAsFactors = F) 
types <- read.table("risk.txt", header=T, sep="\t", check.names=F, stringsAsFactors=F)
types <- types[,c(1, ncol(types)-1)]
rt <- merge(rt, types, how="inner", by="id")
rt <- column_to_rownames(rt, "id")



rt[,"futime"] <- rt[,"futime"]/365

ddist <- datadist(rt)
options(datadist='ddist')



cox <- cph(Surv(futime,fustat) ~age+gender+T+N+M+stage+riskScore,surv=TRUE,x=TRUE, y=TRUE,data=rt)
surv <- Survival(cox)
sur_1_year<-function(x)surv(1,lp=x)#1??????
sur_3_year<-function(x)surv(3,lp=x)#3??????
sur_5_year<-function(x)surv(5,lp=x)#5??????

nom_sur<- nomogram(cox,fun=list(sur_1_year,sur_3_year,sur_5_year),lp= F,
                   funlabel=c('1-Year Survival','3-Year Survival','5-Year survival'),
                   maxscale=100,
                   fun.at=c('0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1'))

pdf("nom.pdf",16,8)
plot(nom_sur)
dev.off()

nom_sur
