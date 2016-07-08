setwd('/Users/Jonathan/Documents/sequencing/RNA-seq/LoaP1/')
library('DESeq2')

count_data = read.csv('GeneCounts.csv', header=TRUE, skip=0, row.names=1)
count_matrix = count_data[,6:16]
sample_info = read.table('sample_info.txt', header=TRUE)

loap_des = DESeqDataSetFromMatrix(count_matrix, sample_info, ~ SampleType)
loap_xyl_des = DESeqDataSetFromMatrix(count_matrix, sample_info, ~ Xylose + LoaP)

loap_des$Strain = relevel(loap_des$SampleType, "WT")
loap_xyl_des$LoaP = relevel(loap_xyl_des$LoaP, "Y")
loap_xyl_des$Xylose = relevel(loap_xyl_des$Xylose, "N")

loap_des = DESeq(loap_des)
loap_xyl_des = DESeq(loap_xyl_des)

loap_blind_vst = varianceStabilizingTransformation(loap_des)
loap_vst = varianceStabilizingTransformation(loap_des, blind=FALSE)

dL_v_WT = results(loap_des, contrast=c('SampleType', 'dL', 'WT'))
Cm_v_Cp = results(loap_des, contrast=c('SampleType', 'CompNoX', 'CompX'))
dL_v_Cp = results(loap_des, contrast=c('SampleType', 'dL', 'CompX'))
Cm_v_WT = results(loap_des, contrast=c('SampleType', 'CompNoX', 'WT'))
Cp_v_WT = results(loap_des, contrast=c('SampleType', 'CompX', 'WT'))


noLoaP_v_LoaP = results(loap_xyl_des)

png(filename="loap_deseq/loaP_PCA.png")
plotPCA(loap_blind_vst, intgroup=c("SampleType")) + ggtitle("loaP PCA") + theme(plot.title = element_text(lineheight=.8, face="bold"))
dev.off()

png(filename="loap_deseq/loap_MA.png")
plotMA(dL_v_WT, alpha=0.05, main="dL vs WT")
dev.off()

png(filename="loap_deseq/loap_comp_MA.png")
plotMA(Cm_v_Cp, alpha=0.05, main="CompNoX vs CompX")
dev.off()

png(filename="loap_deseq/loap_all_MA.png")
plotMA(noLoaP_v_LoaP, alpha=0.05, main="No LoaP v LoaP")
dev.off()



write.csv(as.data.frame(dL_v_WT), file="loap_deseq/dL v WT.csv")
write.csv(as.data.frame(Cm_v_Cp), file="loap_deseq/CompNoX v CompX.csv")
write.csv(as.data.frame(noLoaP_v_LoaP), file="loap_deseq/noLoaP vs LoaP.csv")
