library(tibble)
library(dplyr) 
library(ggplot2) 
library(ggrepel) 
library(tidyr) 
library(BBmisc) 
#library(rtracklayer) only needed for readGFF 


args = commandArgs(trailingOnly=TRUE)
expDir=args[1]
logFC.thres=as.numeric(args[2])
FDR.thres=as.numeric(args[3])
textrepel.thres=as.numeric(args[4]) 

#expDir="Micro.onlyprimary"
#expDir="TUCREnsemblExon"
#expDir="Ensembl"
#logFC.thres=1.5
#FDR.thres=0.001
#textrepel.thres= 20 #-1*log(FDR.thres,10)

samples_data = read.table("defalco.contrast.pairs", header=F)
sample_types = as.character(unique(samples_data[,1]))
sample_sets = as.character(samples_data[,2])
names(sample_sets)=sample_types


#samples_dictionary=read.table("defalco.samples.tvs", header=F) 
#sample_class=samples_dictionary[,1]
#sample_name=samples_dictionary[,2]
#names(sample_class)=sample_name 


#mygff=readGFF("infiles/Homo_sapiens.GRCh38.95.gtf" ) #library(rtracklayer) 
#mygffmin=cbind.data.frame(mygff$gene_id, mygff$gene_name, mygff$gene_biotype)
mygff=read.table("Homo_sapiens.GRCh38.95.gtf.table", header=T )

outdir=paste("/home/enza/oogaProtocol/DeFalco/201910_analisi2/", expDir, "/postEdgeR/", sep="") 
dir.create(file.path(outdir), showWarnings = FALSE)


mycol=c("#fad284", "#f38181", "#705772")


all <- data.frame(gene_id=character(),
                  sampleA=character(), 
                  sampleA=character(),
                  logFC=numeric(),
                  logCPM=numeric(),
                  PValue=numeric(),
                 FDR=numeric(),
                 experiment=character(),
                 set=character(), 
                 stringsAsFactors=TRUE)

for (st in sample_types){ 
   
    #~~~~  1. load Differential Expressin results 
    myde=read.table(paste(expDir, "/", expDir, ".counts.format.mat.", st,".edgeR.DE_results", sep=""), header=T , row.names=1 ) 
    mydenames=rownames_to_column(myde, var = "gene_id")
    mydenames$experiment=st
    mydenames$set=sample_sets[[st]]
    mydenames$posneg=ifelse(mydenames$logFC>0 , "y", "n")
   
   
    #~~~~  2. List of DE genes 
    DEgenes<- mydenames %>% filter (FDR< FDR.thres & abs(logFC)>logFC.thres ) %>% select("gene_id")  
 
   
    colorbar=c( "#2a1a5e", "#f45905")
    #~~~~  3. barplots per experiment 
    mydbar.tmp=merge (DEgenes, mydenames,  by="gene_id", all.x=F )
    mydbar = merge (mydbar.tmp, mygff,  by="gene_id", all.x=F )
    
    
    #ggplot(mydbar, aes(x = reorder(gene_name, logFC) , logFC) )+geom_bar(stat="identity", aes(fill=experiment), position = position_dodge(preserve = 'single'))+coord_flip()  +theme_bw() + scale_fill_manual(values=mycol) +ylab (paste("Fold change (FRD <",FDR.thres,")")) +xlab("")+ theme(legend.title = element_blank()) +ggtitle( paste(expDir, st, " FDR <", FDR.thres, sep=" " )) 
    bp<-ggplot(mydbar.tmp, aes(x = reorder(gene_id, logFC) , logFC) )+geom_bar(stat="identity", aes(fill=posneg), position = position_dodge(preserve = 'single'))+coord_flip()  +theme_bw() +ylab (paste("Fold change (FRD <",FDR.thres,")")) +xlab("")+ theme(legend.title = element_blank()) +ggtitle( paste("microRNA Primary Transcripts", st, "\nFDR <", FDR.thres, "abs(logFC) >", logFC.thres, sep=" " )) +
    theme(legend.position = "none", legend.title=element_blank(),legend.box.background=element_blank(),
      axis.title.x = element_text(face=1, size=20),
      axis.text.x  = element_text( size=20),
      axis.title.y = element_text(face=1, size=20),
      axis.text.y  = element_text(vjust=0.5, size=20), 
      plot.title = element_text(size = 25)) + scale_fill_manual(values= colorbar)  
      
    png(paste(outdir,"bar.", st, ".png", sep=""), res=300 ,width=30, height=20  ,units="cm")
    print(bp)
    dev.off()          
    
    
    #~~~~  4. heatmap per experiment 
    #~~~~  4.1. load raw counts to normalize 
    myrawcounts=read.table(paste(expDir, "/", expDir, ".counts.format.mat.", st, ".edgeR.count_matrix", sep=""), header=T , row.names=1) 
    myrawcountsnames=rownames_to_column(myrawcounts, var = "gene_id")  
    mycounts.norm <- gather(as.data.frame(cbind.data.frame (gene_id=myrawcountsnames$gene_id ,  normalize(myrawcountsnames[2:ncol(myrawcountsnames)], method = "standardize")  )), "sample", "counts", -gene_id)
    #mycounts.norm <- gather(as.data.frame(cbind.data.frame (gene_id=myrawcountsnames$gene_id ,  normalize(myrawcountsnames[2:ncol(myrawcountsnames)], method = "range", range = c(-1, 1))  )), "sample", "counts", -gene_id)
    #mycounts.norm <- gather(as.data.frame(cbind.data.frame (gene_id=myrawcountsnames$gene_id ,  scale(myrawcountsnames[2:ncol(myrawcountsnames)])  )), "sample", "counts", -gene_id)
    #~~~~  4.2. heatmap  
    mydheat.tmp= merge (DEgenes, mycounts.norm,  by="gene_id", all.x=F )
    mydheat = merge (mydheat.tmp, mygff,  by="gene_id", all.x=F )
    hp<- ggplot(mydheat.tmp, aes(x = reorder(gene_id, counts), sample) ) +geom_tile(aes(fill=counts) ) +theme_bw() + scale_fill_gradient2(low = "#2a1a5e", high = "#f45905")  +ggtitle( paste( st, " FDR <", FDR.thres, "abs(logFC) >", logFC.thres,  sep=" " ) ) +xlab("") + ylab("")+
theme(legend.box.background=element_blank(),
      axis.title.x = element_text(face=1, size=20),
      axis.text.x  = element_text( size=18, angle =0),
      axis.title.y = element_text(face=1, size=20),
      axis.text.y  = element_text(vjust=0.5, size=20),
      plot.title = element_text(size = 25))+ labs(fill = "standardized\ncounts")
    
    
    png(paste(outdir, "heat.", st , ".png", sep=""), res=300 ,width=40, height=15  ,units="cm")
    print( hp)
     dev.off() 
    
       
    #~~~~  5. fill the dataframe for plots by set (out of the cycle) 
    all=rbind(all, mydenames) 
    
    
    #~~~~  6.write reports for single experiments 
    mycounts=read.table(paste(expDir, "/", expDir, ".counts.format.mat.", st, ".edgeR.count_matrix", sep=""), header=T , row.names=1, col.names=c("condArep1", "condArep2", "condArep3", "condBrep1", "condBrep2", "condBrep3")   )
    mycountsnames=rownames_to_column(mycounts, var = "gene_id")
    report.tmp=merge(mydenames, mycountsnames, by="gene_id")
    report=merge(report.tmp, mygff, by="gene_id", all.x=F)
    write.table(report, file = paste(outdir,"report.",st , ".tsv", sep="") , append = FALSE, quote = FALSE,  sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE)
}



#all %>% mutate(ucname=paste("uc.", strsplit(gene_id, "[_]|[.]")[[1]][3], sep="")) %>% head()  ##rename uc 
all$col <- as.factor(ifelse(abs(all$logFC)>logFC.thres & all$FDR<FDR.thres, "disregulated", "unchanged"))
all$updown=ifelse(all$col=="disregulated", ifelse(all$logFC>0, "upregulated", "downregulated" ) , "unchanged")
mytally=all %>% group_by(set, experiment, updown) %>% tally()

write.table(mytally, file = paste(outdir, "tally.logFC", logFC.thres,".FDR" ,FDR.thres, ".tsv", sep="") , append = FALSE, quote = FALSE,  sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE)


all$experiment <- factor(all$experiment, levels=c("NOX_vs_Parental", "ACID_vs_Parental", "ACID_vs_NOX", "NOX48HACID_vs_NOX48H", "NOX24HHYP_vs_NOX48H" ))

### VOLCANO 
#png(paste(outdir,"set1.volc.png", sep=""), res=300 ,width=45, height=15  ,units="cm")
#ggplot(subset(all, set=="set1"), aes(logFC, -1*log(FDR, 10) , colour=col) ) +geom_point( )+facet_grid (.  ~ experiment, drop=TRUE  ) +theme_bw() +geom_vline(xintercept=c(-logFC.thres, logFC.thres), color="grey" ) +geom_hline(yintercept=-1*log(FDR.thres,10), color="grey") +scale_color_manual(values=c("#3c415e" ,"#dfe2e2") ) +geom_text_repel(data=subset(all, -1*log(FDR,10) > textrepel.thres &  abs(logFC) >logFC.thres  & set=="set1"  ), aes(label=gene_id) )+ theme(legend.position="none", strip.text= element_text(size = 15) ) +ggtitle( paste(expDir, st, " FDR <", FDR.thres, sep=" " ))
#dev.off()

png(paste(outdir,"set1.volc.png", sep=""), res=300 ,width=45, height=15  ,units="cm")

ggplot(subset(all, set=="set1"), aes(logFC, -1*log(FDR, 10) , colour=col) ) + 
geom_point( )+ 
facet_grid (.  ~ experiment, drop=TRUE  ) +
theme_bw() +
geom_vline(xintercept=c(-logFC.thres, logFC.thres), color="grey" ) + 
geom_hline(yintercept=-1*log(FDR.thres,10), color="grey") + 
scale_color_manual(values=c("#3c415e" ,"#dfe2e2") ) +
geom_text_repel(data=subset(all, -1*log(FDR,10) > textrepel.thres &  abs(logFC)>logFC.thres  & set=="set1"  ), aes(label=gene_id) ) + 
ylab("-log10 FDR") +
ggtitle( paste("microRNA Primary Transcripts", st , sep=" " )) + 
theme(legend.position="none", 
	strip.text= element_text(size = 20),
	axis.title.x = element_text(face=1, size=20), 
	axis.text.x  = element_text( size=20, angle =0),
	axis.title.y = element_text(face=1, size=20),
	axis.text.y  = element_text(vjust=0.5, size=20), 
	plot.title = element_text(size = 25))
dev.off()


### BAR 

png(paste(outdir,"set1.bar.png", sep=""), res=300 ,width=15, height=20  ,units="cm")
ggplot(subset(all, FDR<FDR.thres &  abs(logFC)>logFC.thres & set=="set1" ) , aes(x = reorder(gene_id, logFC) , logFC) )+geom_bar(stat="identity", aes(fill=experiment), position = position_dodge(preserve = 'single'))+coord_flip()  +theme_bw() + scale_fill_manual(values=mycol) +ylab (paste("Fold change (FRD <",FDR.thres,")")) +xlab("")+ theme(legend.title = element_blank()) +ggtitle( paste(expDir, st, " FDR <", FDR.thres, sep=" " ))
dev.off() 


png(paste(outdir,"set2.bar.png", sep=""), res=300 ,width=15, height=20  ,units="cm")
ggplot(subset(all, FDR<FDR.thres & abs(logFC)>logFC.thres & set=="set2") , aes(x = reorder(gene_id, logFC) , logFC) )+geom_bar(stat="identity", aes(fill=experiment), position = position_dodge(preserve = 'single'))+coord_flip()  +theme_bw() + scale_fill_manual(values=mycol) +ylab (paste("Fold change (FRD <",FDR.thres,")")) +xlab("")+ theme(legend.title = element_blank()) +ggtitle( paste(expDir, st, " FDR <", FDR.thres, sep=" " ))
dev.off() 
