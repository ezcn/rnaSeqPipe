####~~~~~ take gene names from gtf file ~~~~~~#####
#zcat Homo_sapiens.GRCh38.95.gtf.gz | awk ' {if ($3=="gene") print  $10,  $14 } ' | tr "\"" " "  > gene_name.tsv

library(tidyverse)
library(gplots)
library(ggrepel)
library(reshape2)


#1) Plot reads type

myd<-read.table("countsLonardo.txt.summary", sep="\t", header=T)
mydl <- gather(myd, key=readtype , value=number, -Status)
mydl$readtype=factor(mydl$readtype, levels=c("CTR.1", "CTR.2", "CTR.3", "NODAL.1", "NODAL.2", "NODAL.3", "POS.1", "POS.2", "POS.3", "NEG.1", "NEG.2", "NEG.3"))

pread<-ggplot(mydl, aes(readtype,  number, fill=Status) ) +geom_bar(stat="identity" )+theme_bw()+ggtitle("countset1 - -p -B -C -t exon -g gene_id ")+xlab("") +ylab("number of reads") + theme(axis.text.x = element_text(angle = 90))

ggsave("/home/silvia/rna_seq/results/reads_type.png", plot= pread, device="png", width = 40, height = 23, units = "cm", dpi = 300)



#2) Add gene names to final files 

myd<-read.table("NEW_Lonardo_counts.txt.CTR_vs_NODAL.edgeR.DE_results", header=T, sep="\t")
myn<-read.table("/home/silvia/Downloads/gene_name.tsv", sep="\t", header=T)
tot<-merge(myd,myn, by="geneid")
write.table(tot, "Lonardo_CTR_Vs_NODAL_results", quote=F, sep="\t" , row.names=F, col.names=T)


myd<-read.table("NEW_Lonardo_counts.txt.POS_vs_NEG.edgeR.DE_results", header=T, sep="\t")
myn<-read.table("/home/silvia/Downloads/gene_name.tsv", sep="\t", header=T)
tot<-merge(myd,myn, by="geneid")
write.table(tot, "Lonardo_POS_Vs_NEG.DE_results", quote=F, sep="\t" , row.names=F, col.names=T)


#3) Add read counts to final files

myc<-read.table("/home/silvia/rna_seq/samples/Lonardo_counts.txt", sep="\t", header=T)
myc2<- myc %>% select (geneid, CTR.1, CTR.2, CTR.3, NODAL.1, NODAL.2, NODAL.3 )
final<-merge(tot,myc2, by="geneid")

write.table(final, "/home/silvia/rna_seq/results/Lonardo_NODAL_Vs_CTR_withcounts.DE_results", quote=F, sep="\t" , row.names=F, col.names=T)

myc<-read.table("/home/silvia/rna_seq/samples/Lonardo_counts.txt", sep="\t", header=T)
myc2<- myc %>% select (geneid, POS.1, POS.2, POS.3, NEG.1, NEG.2, NEG.3 )
final<-merge(tot,myc2, by="geneid")

write.table(final, "/home/silvia/rna_seq/results/Lonardo_POS_Vs_NEG_withcounts.DE_results", quote=F, sep="\t" , row.names=F, col.names=T)

#4) Check for FDR 

ggplot(final, aes(FDR)) + geom_density()

#5) Heat map on counts

myd<-read.table("/mpbastudies3/lonardo/Lonardo_counts.txt", sep="\t", header=T)
mydf<- myd %>% select( CTR.1, CTR.2, CTR.3, NODAL.1, NODAL.2, NODAL.3, POS.1, POS.2, POS.3, NEG.1, NEG.2, NEG.3)
myd1<- mydf[ rowSums(mydf > 0) >= 12, ]
rnames <- myd1[,1]
mat_data <- data.matrix(myd1[,2:ncol(myd1)])
rownames(mat_data) <- rnames
png("/mpbastudies3/lonardo/counts3.png", width = 40, height = 23, units = "cm", res = 300)
heatmap(mat_data, Rowv=NA)
dev.off()

#6) DE plots

myd<-read.table("/home/silvia/rna_seq/results/Lonardo_POS_Vs_NEG_withcounts.DE_results", sep="\t", header=T)

logFC.thres=1.5
FDR.thres=0.05
myd$col <- as.factor(ifelse(abs(myd$logFC)>logFC.thres & myd$FDR<FDR.thres, "disregulated", "unchanged"))

ggplot(myd, aes(logFC, -1*log(FDR, 10) , colour=col) ) +geom_point( ) +theme_bw() +geom_vline(xintercept=c(-logFC.thres, logFC.thres), color="grey" ) +geom_hline(yintercept=-1*log(FDR.thres,10), color="grey") +scale_color_manual(values=c("#3C415E" ,"#DFE2E2") ) +ggtitle("POSITIVE versus NEGATIVE -  FDR=0.05 logFC=1.5" ) +geom_text_repel(data=subset(myd, col=="disregulated"),  aes(label=name) )

ggsave("POSNEG.png", device="png", width = 30, height = 30, units = "cm", dpi = 300)
