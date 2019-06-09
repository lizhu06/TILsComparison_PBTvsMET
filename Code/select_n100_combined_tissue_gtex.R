rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")

## Gtex tpm
load("Data/gtex_gene_annot_uni_genes.RData")
load("Data/gtex_tpm_uni_genes.RData")
load("Data/gtex_tissues.RData")

tissue_short <- sapply(1:length(tissues), function(x) strsplit(tissues[x], split=" - ")[[1]][1])
names(tissue_short) <- names(tissues)
gene_annot <- gtex_gene_annot_uni_genes
log2tpm <- data.matrix(log2(gtex_tpm_uni_genes+1))

rm(gtex_tpm_uni_genes)

all(colnames(log2tpm) == names(tissue_short))

id_keep <- tissue_keep <- NULL
uni_tissue <- unique(tissue_short)
length(uni_tissue) #31

set.seed(201703)
for(i in 1:length(uni_tissue)){
	if(sum(tissue_short == uni_tissue[i]) < 100) {
		id_keep <- c(id_keep, which(tissue_short == uni_tissue[i]))
		tissue_keep <- c(tissue_keep, rep(uni_tissue[i], sum(tissue_short == uni_tissue[i])))
		print(paste("keep ", sum(tissue_short == uni_tissue[i]), " ", uni_tissue[i], sep=""))
	}else{
		id_keep <- c(id_keep, sample(which(tissue_short == uni_tissue[i]), size=100, replace=FALSE))
		tissue_keep <- c(tissue_keep, rep(uni_tissue[i], 100))
		print(paste("keep 100 ", uni_tissue[i], sep=""))
	}
}

length(id_keep) # 2771
log2tpm_keep <- log2tpm[, id_keep]
names(tissue_keep) <- colnames(log2tpm_keep)
length(tissue_keep) # 2771

id_keep <- colnames(log2tpm_keep)
save(id_keep, file="Data/Gtex_n100_combined_tissues_id_keep.RData")
save(tissue_keep, file="Data/Gtex_n100_combined_tissues_tissues_keep.RData")
save(log2tpm_keep, tissue_keep, file="Data/Gtex_n100_combined_tissues_log2tpm_keep.RData")



######## NOT RUN BELOW ########
## PCA plot 
library(ggplot2)
library(grid)
library(gridExtra)

## filter genes based on SD
sd <- apply(log2tpm_keep, 1, sd)
log2tpm_keep_fil <- log2tpm_keep[sd > 0, ]
dim(log2tpm_keep_fil) # 53642  2771

# create data frame with scores
pca_res = prcomp(t(log2tpm_keep_fil), scale= TRUE)
df_out = as.data.frame(pca_res$x)

df_out$labels <- factor(tissue_keep)

# plot of observations
percentage <- round(pca_res$sdev / sum(pca_res$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )

pdf(file="Results/PCA_Gtex_n100_combined_tissues.pdf", width=15, height=5)
p<-ggplot(df_out, aes(x=PC1, y=PC2, color=labels, shape=labels))+ 
	scale_shape_manual(values=seq(1,length(df_out$labels))+60)+
	geom_point()+xlab(percentage[1]) + ylab(percentage[2])
p
dev.off()
