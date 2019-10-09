# rnaSeqPipe

lab pipeline for RNA seq analyses 



1. alignment using BWA 

2. counts using featurecounts 

3. DE using EdgeR 

perl scr/run_DE_analysis.pl  --matrix infiles/UCbase.nostrand.micro.counts.format.mat  --method edgeR --dispersion 0.1 --samples_file defalco.samples.tvs --contrasts defalco.contrast.tsv  --output UCbase.nostrand.micro

4. postedgeR

library(tibble)
library(dplyr) 
library(ggplot2) 
library(ggrepel) 
library(tidyr) 

samples_data = read.table("defalco.contrast.pairs", header=F)
sample_types = as.character(unique(samples_data[,1]))
sample_sets = as.character(samples_data[,2])
names(sample_sets)=sample_types

#expDir="Micro.onlyprimary"
#expDir="TUCREnsemblExon"
expDir="Ensembl"
outdir=paste("/home/enza/oogaProtocol/DeFalco/201910_analisi2/", expDir, "/postEdgeR/", sep="") 
dir.create(file.path(outdir), showWarnings = FALSE)


logFC.thres=1.5
FDR.thres=0.0001

mycol=c("#fad284", "#f38181", "#705772")


all <- data.frame(ensid=character(),
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
    
    ##  raw counts 
   
    mycounts=read.table(paste(expDir, "/", expDir, ".counts.format.mat.", st, ".edgeR.count_matrix", sep=""), header=T , row.names=1, col.names=c("condArep1", "condArep2", "condArep3", "condBrep1", "condBrep2", "condBrep3")   )
    mycountsnames=rownames_to_column(mycounts, var = "ensid")

    ##  raw counts to normalize 
    myrawcounts=read.table(paste(expDir, "/", expDir, ".counts.format.mat.", st, ".edgeR.count_matrix", sep=""), header=T , row.names=1) 
    myrawcountsnames=rownames_to_column(myrawcounts, var = "ensid")
    mycounts.norm <- gather(as.data.frame(cbind.data.frame (ensid=myrawcountsnames$ensid ,  scale(myrawcountsnames[2:ncol(myrawcountsnames)])  )), "sample", "counts", -ensid)

    ##  DE results 
    myde=read.table(paste(expDir, "/", expDir, ".counts.format.mat.", st,".edgeR.DE_results", sep=""), header=T , row.names=1 ) 
    mydenames=rownames_to_column(myde, var = "ensid")
    mydenames$experiment=st
    mydenames$set=sample_sets[[st]]
    
    ##  List of DE genes 
    DEgenes<- mydenames %>% filter (FDR< FDR.thres) %>% select("ensid")  
    mydheat= merge (DEgenes, mycounts.norm, by="ensid" )
    png(paste(outdir, "heat.", st , ".png", sep=""), res=300 ,width=40, height=15  ,units="cm") 
    print(ggplot(mydheat, aes(x = reorder(ensid, counts), sample) ) +geom_tile(aes(fill=counts) ) +theme_bw() + scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) + theme(axis.text.x = element_text(angle = 90))) +ggtitle( paste(expDir, st, " FDR <", FDR.thres, sep=" " ))
    dev.off() 
    
       ## BAR 
    png(paste(outdir,"bar.", st, ".png", sep=""), res=300 ,width=15, height=20  ,units="cm")
    print(ggplot(subset(mydenames, FDR<FDR.thres) , aes(x = reorder(ensid, logFC) , logFC) )+geom_bar(stat="identity", aes(fill=experiment), position = position_dodge(preserve = 'single'))+coord_flip()  +theme_bw() + scale_fill_manual(values=mycol) +ylab (paste("Fold change (FRD <",FDR.thres,")")) +xlab("")+ theme(legend.title = element_blank()) +ggtitle( paste(expDir, st, " FDR <", FDR.thres, sep=" " )) )
    dev.off()   
       
    ## fill the dataframe
    all=rbind(all, mydenames) 
    
    ## write reports for single experiments 
    report=merge(mydenames, mycountsnames, by="ensid")
    write.table(report, file = paste(outdir,"report.",st , ".tsv", sep="") , append = FALSE, quote = FALSE,  sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE)
}



#all %>% mutate(ucname=paste("uc.", strsplit(ensid, "[_]|[.]")[[1]][3], sep="")) %>% head()  ##rename uc 
all$col <- as.factor(ifelse(abs(all$logFC)>logFC.thres & all$FDR<FDR.thres, "disregulated", "unchanged"))
all$updown=ifelse(all$col=="disregulated", ifelse(all$logFC>0, "upregulated", "downregulated" ) , "unchanged")
mytally=all %>% group_by(set, experiment, updown) %>% tally()

write.table(mytally, file = paste(outdir, "tally.logFC", logFC.thres,".FDR" ,FDR.thres, ".tsv", sep="") , append = FALSE, quote = FALSE,  sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE)


all$experiment <- factor(all$experiment, levels=c("NOX_vs_Parental", "ACID_vs_Parental", "ACID_vs_NOX", "NOX48HACID_vs_NOX48H", "NOX24HHYP_vs_NOX48H" ))

### VOLCANO 
png(paste(outdir,"set1.volc.png", sep=""), res=300 ,width=45, height=15  ,units="cm")
ggplot(subset(all, set=="set1"), aes(logFC, -1*log(FDR, 10) , colour=col) ) +geom_point( )+facet_grid (.  ~ experiment, drop=TRUE  ) +theme_bw() +geom_vline(xintercept=c(-logFC.thres, logFC.thres), color="grey" ) +geom_hline(yintercept=-1*log(FDR.thres,10), color="grey") +scale_color_manual(values=c("#3c415e" ,"#dfe2e2") ) +geom_text_repel(data=subset(all, -1*log(FDR,10) > -1*log(FDR.thres,10) & set=="set1"  ), aes(label=ensid) )+ theme(legend.position="none", strip.text= element_text(size = 15) ) +ggtitle( paste(expDir, st, " FDR <", FDR.thres, sep=" " ))
dev.off()

png(paste(outdir,"set2.volc.png", sep=""), res=300 ,width=45, height=15  ,units="cm")
ggplot(subset(all, set=="set2"), aes(logFC, -1*log(FDR, 10) , colour=col) ) +geom_point( )+facet_grid (.  ~ experiment, drop=TRUE  ) +theme_bw() +geom_vline(xintercept=c(-logFC.thres, logFC.thres), color="grey" ) +geom_hline(yintercept=-1*log(FDR.thres,10), color="grey") +scale_color_manual(values=c("#3c415e" ,"#dfe2e2") ) +geom_text_repel(data=subset(all, -1*log(FDR,10) > -1*log(FDR.thres,10) & set=="set2"  ), aes(label=ensid) )+ theme(legend.position="none", strip.text= element_text(size = 15) ) +ggtitle( paste(expDir, st, " FDR <", FDR.thres, sep=" " ))
dev.off()


### BAR 

png(paste(outdir,"set1.bar.png", sep=""), res=300 ,width=15, height=20  ,units="cm")
ggplot(subset(all, FDR<FDR.thres & set=="set1") , aes(x = reorder(ensid, logFC) , logFC) )+geom_bar(stat="identity", aes(fill=experiment), position = position_dodge(preserve = 'single'))+coord_flip()  +theme_bw() + scale_fill_manual(values=mycol) +ylab (paste("Fold change (FRD <",FDR.thres,")")) +xlab("")+ theme(legend.title = element_blank()) +ggtitle( paste(expDir, st, " FDR <", FDR.thres, sep=" " ))
dev.off() 


png(paste(outdir,"set2.bar.png", sep=""), res=300 ,width=15, height=20  ,units="cm")
ggplot(subset(all, FDR<FDR.thres & set=="set2") , aes(x = reorder(ensid, logFC) , logFC) )+geom_bar(stat="identity", aes(fill=experiment), position = position_dodge(preserve = 'single'))+coord_flip()  +theme_bw() + scale_fill_manual(values=mycol) +ylab (paste("Fold change (FRD <",FDR.thres,")")) +xlab("")+ theme(legend.title = element_blank()) +ggtitle( paste(expDir, st, " FDR <", FDR.thres, sep=" " ))
dev.off() 
