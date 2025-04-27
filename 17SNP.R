
library(maftools)

risk <- read.table("risk.txt",sep="\t",header=T,check.names=F,row.names=1)
maf<- read.maf(maf = 'input.maf')
# maf@data$Tumor_Sample_Barcode 
mut <- maf@data
mut$Tumor_Sample_Barcode <- substring(mut$Tumor_Sample_Barcode,1,12)

mut.High <- mut[(mut$Tumor_Sample_Barcode %in% rownames(risk)[risk$risk=="high"]),]
mut.Low <- mut[(mut$Tumor_Sample_Barcode %in% rownames(risk)[risk$risk=="low"]),]
maf.High <-read.maf(maf = mut.High,isTCGA = T)
maf.Low <- read.maf(maf=mut.Low,isTCGA=T)

# 颜色
col = RColorBrewer::brewer.pal(n=10, name='Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Ins','In_Frame_Ins', 'Splice_Site', 'In_Frame_Del','Nonstop_Mutation','Translation_Start_Site','Multi_Hit')

#人种
racecolors = RColorBrewer::brewer.pal(n=4, name='Spectral')
names(racecolors) = c("ASIAN", "WHITE", "BLACK_OR_AFRICAN_AMERICAN",  "AMERICAN_INDIAN_OR_ALASKA_NATIVE")


pdf(file="summary_high.pdf",width=7,height=6)
plotmafSummary(maf=maf.High, rmOutlier=TRUE, addStat='median', dashboard=TRUE, titvRaw=FALSE)
dev.off()

pdf(file="summary_low.pdf",width=7,height=6)
plotmafSummary(maf=maf.Low, rmOutlier=TRUE, addStat='median', dashboard=TRUE, titvRaw=FALSE)
dev.off()

pdf(file="waterfall_high.pdf",width=7,height=6)
oncoplot(maf=maf.High, top=20, fontSize=0.8 ,showTumorSampleBarcodes=F)
dev.off()

pdf(file="waterfall_low.pdf",width=7,height=6)
oncoplot(maf=maf.Low, top=20, fontSize=0.8, showTumorSampleBarcodes=F)
dev.off()


pdf(file="interaction_high.pdf",width=7,height=6)
somaticInteractions(maf=maf.High, top=20, pvalue=c(0.05, 0.001), showCounts=FALSE, fontSize=0.6)
dev.off()

pdf(file="interaction_low.pdf",width=7,height=6)
somaticInteractions(maf=maf.Low, top=20, pvalue=c(0.05, 0.001), showCounts=FALSE, fontSize=0.6)
dev.off()