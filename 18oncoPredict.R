rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)

dir='./DataFiles/Training Data/'
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 


# testExpr<- GDSC2_Expr[,sample(1:ncol(GDSC2_Expr),10)]
# testExpr[1:4,1:4]
testExpr <- read.table("uniq.symbol.txt",sep="\t",header=T,check.names=F,row.names=1)
testExpr <- testExpr[,60:ncol(testExpr)] %>% as.matrix()
# colnames(testExpr)=paste0('test',colnames(testExpr))
dim(testExpr)


calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )