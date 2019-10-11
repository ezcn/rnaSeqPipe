library(cluster)
library(Biobase)
library(qvalue)
NO_REUSE = F


args = commandArgs(trailingOnly=TRUE)
mycounts=args[1]
mysamples=args[2]
print(mycounts) 

# try to reuse earlier-loaded data if possible
#if (file.exists(paste(mycounts, ".RData", sep="" )) && ! NO_REUSE) {
#	    print('RESTORING DATA FROM EARLIER ANALYSIS')
#    load("mycounts")
#} else {
#	    print('Reading matrix file.')
 #   primary_data = read.table("countsIGB.txt.tableTRI.txt", header=T, com='', sep="\t", row.names=1, check.names=F)
#        primary_data = as.matrix(primary_data)
#}

#source("/usr/local/src/trinityrnaseq-Trinity-v2.3.2/Analysis/DifferentialExpression/R/heatmap.3.R")
#source("/usr/local/src/trinityrnaseq-Trinity-v2.3.2/Analysis/DifferentialExpression/R/misc_rnaseq_funcs.R")
#source("/usr/local/src/trinityrnaseq-Trinity-v2.3.2/Analysis/DifferentialExpression/R/pairs3.R")

source("/home/enza/oogaProtocol/DeFalco/201910_analisi2/trinitysources/heatmap.3.R")
source("/home/enza/oogaProtocol/DeFalco/201910_analisi2/trinitysources/misc_rnaseq_funcs.R")
source("/home/enza/oogaProtocol/DeFalco/201910_analisi2/trinitysources/pairs3.R")

###### counts 
primary_data = read.table(mycounts, header=T, com='', sep="\t", row.names=1, check.names=F)
data = primary_data

##### samples 
samples_data = read.table(mysamples, header=F, check.names=F, fill=T)
samples_data = samples_data[samples_data[,2] != '',]
sample_types = as.character(unique(samples_data[,1]))
rep_names = as.character(samples_data[,2])

data = data[, colnames(data) %in% samples_data[,2], drop=F ]
nsamples = length(sample_types)
sample_colors = rainbow(nsamples)
names(sample_colors) = sample_types
sample_type_list = list()
for (i in 1:nsamples) {
	    samples_want = samples_data[samples_data[,1]==sample_types[i], 2]
    sample_type_list[[sample_types[i]]] = as.vector(samples_want)
}
sample_factoring = colnames(data)
for (i in 1:nsamples) {
	    sample_type = sample_types[i]
    replicates_want = sample_type_list[[sample_type]]
        sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
# reorder according to sample type.
tmp_sample_reordering = order(sample_factoring)
data = data[,tmp_sample_reordering,drop=F]
sample_factoring = sample_factoring[tmp_sample_reordering]
initial_matrix = data # store before doing various data transformations
data = log2(data+1)
sample_factoring = colnames(data)
for (i in 1:nsamples) {
	    sample_type = sample_types[i]
    replicates_want = sample_type_list[[sample_type]]
        sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
sampleAnnotations = matrix(ncol=ncol(data),nrow=nsamples)
for (i in 1:nsamples) {
	  sampleAnnotations[i,] = colnames(data) %in% sample_type_list[[sample_types[i]]]
}
sampleAnnotations = apply(sampleAnnotations, 1:2, function(x) as.logical(x))
sampleAnnotations = sample_matrix_to_color_assignments(sampleAnnotations, col=sample_colors)
rownames(sampleAnnotations) = as.vector(sample_types)
colnames(sampleAnnotations) = colnames(data)
data = as.matrix(data) # convert to matrix
write.table(data, file=paste(mycounts, ".log2.dat", sep=""), quote=F, sep="\t");
pdf(paste(mycounts,"log2.pc.pdf", sep=""))
data = as.matrix(data)
# Z-scale and center the genes across all the samples for PCA
prin_comp_data = initial_matrix
prin_comp_data = log2(prin_comp_data+1)
prin_comp_data = scale(prin_comp_data)
prin_comp_data = t(scale(t(prin_comp_data), center=TRUE, scale=F)) # just center trans expr level, retain original effect size.
pca = prcomp(t(prin_comp_data), center = FALSE, scale. = FALSE)
pc_pct_variance = (pca$sdev^2)/sum(pca$sdev^2)
def.par <- par(no.readonly = TRUE) # save default, for resetting...
gridlayout = matrix(c(1:4),nrow=2,ncol=2, byrow=TRUE);
layout(gridlayout, widths=c(1,1));
write.table(pca$rotation, file=paste(mycounts, ".log2.PCA.loadings", sep=""), quote=F, sep="\t")
write.table(pca$x, file=paste(mycounts, ".log2.PCA.scores", sep=""), quote=F, sep="\t")
for (i in 1:(max(3,2)-1)) {
	    xrange = range(pca$x[,i])
    yrange = range(pca$x[,i+1])
        samples_want = rownames(pca$x) %in% sample_type_list[[sample_types[1]]]
        pc_i_pct_var = sprintf("(%.2f%%)", pc_pct_variance[i]*100)
	    pc_i_1_pct_var = sprintf("(%.2f%%)", pc_pct_variance[i+1]*100)
	    plot(pca$x[samples_want,i], pca$x[samples_want,i+1], xlab=paste('PC',i, pc_i_pct_var), ylab=paste('PC',i+1, pc_i_1_pct_var), xlim=xrange, ylim=yrange, col=sample_colors[1])
	        for (j in 2:nsamples) {
			        samples_want = rownames(pca$x) %in% sample_type_list[[sample_types[j]]]
	            points(pca$x[samples_want,i], pca$x[samples_want,i+1], col=sample_colors[j], pch=j)
		        }
	        plot.new()
		    legend('topleft', as.vector(sample_types), col=sample_colors, pch=1:nsamples, ncol=2)
}

par(def.par)
pcscore_mat_vals = pca$rotation[,1:3]
pcscore_mat = matrix_to_color_assignments(pcscore_mat_vals, col=colorpanel(256,'purple','black','yellow'), by='row')
colnames(pcscore_mat) = paste('PC', 1:ncol(pcscore_mat))
dev.off()
gene_cor = NULL
