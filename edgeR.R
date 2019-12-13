if (! require(edgeR)) {
   source("https://bioconductor.org/biocLite.R")
   biocLite("edgeR")
   library(edgeR)
}

data = read.table("/home/silvia/rna_seq/samples/Lonardo_counts.txt", header=T, row.names=1, com='')
col_ordering = c(4,5,6,1,2,3)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
conditions = factor(c(rep("NODAL", 3), rep("CTR", 3)))

exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study = calcNormFactors(exp_study)
exp_study = estimateCommonDisp(exp_study)
exp_study = estimateTagwiseDisp(exp_study)
et = exactTest(exp_study, pair=c("NODAL", "CTR"))
tTags = topTags(et,n=NULL)
result_table = tTags$table
result_table = data.frame(sampleA="NODAL", sampleB="CTR", result_table)
result_table$logFC = -1 * result_table$logFC
write.table(result_table, file='Lonardo_counts.txt.NODAL_vs_CTR.edgeR.DE_results', sep='	', quote=F, row.names=T)
write.table(rnaseqMatrix, file='Lonardo_counts.txt.NODAL_vs_CTR.edgeR.count_matrix', sep='	', quote=F, row.names=T)
source("/home/silvia/rna_seq/R/rnaseq_plot_funcs.R")
pdf("Lonardo_counts.txt.NODAL_vs_CTR.edgeR.DE_results.MA_n_Volcano.pdf")
plot_MA_and_Volcano(rownames(result_table), result_table$logCPM, result_table$logFC, result_table$FDR)
dev.off()
