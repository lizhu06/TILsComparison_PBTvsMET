rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")
# R 3.4.0
library(limma) #3.32.6
library(GSVA)  # 1.24.2

## load  scores
davolin_ave <- get(load("Data/ave_score_davoli_tpm_Chung_pooled.RData"))
davolin_gsva <- get(load("Data/gsva_score_davoli_tpm_Chung_pooled.RData"))
tamborero_ave <- get(load("Data/ave_score_tamborero_tpm_Chung_pooled.RData"))
tamborero_gsva <- get(load("Data/gsva_score_tamborero_tpm_Chung_pooled.RData"))
ciber_abs <- read.csv("Results/cibersort_res_chuang_absolute.csv")
ciber_abs <- t(ciber_abs)
colnames(ciber_abs) <- ciber_abs[1,]
ciber_abs <- (ciber_abs[-1,])
ciber_abs <- ciber_abs[seq(1,22),]

ciber_rel <- read.csv("Results/cibersort_res_chuang_relative.csv")
ciber_rel <- t(ciber_rel)
colnames(ciber_rel) <- ciber_rel[1,]
ciber_rel <- (ciber_rel[-1,])
ciber_rel <- ciber_rel[seq(1,22),]

## manually load in cell type proportion
count <- c(9, 38, 8, 3, 23, 2, 
	9, 6, 4, 7, 2, 27, 
	3, 10, 1, 22,2)
bc1 <- c(0, 0, 0)
bc2 <- c(0, 0, 0)
bc3 <- c(9, 9, 0)
bc3ln <- c(38, 6, 0)
bc04 <- c(8, 0, 0)
bc5 <- c(0, 0, 0)
bc06 <- c(0, 4, 3)
bc07 <- c(3, 7, 10)
bc07ln <- c(23, 2, 0)
bc08 <- c(0, 0, 1)
bc09 <- c(2, 27, 22)
bc10 <- c(0, 0, 0)
count <- cbind(bc1, bc2, bc3, bc3ln, bc04, bc5, bc06, bc07ln, bc08, bc09, bc10)
colnames(count) <- c("BC01","BC02","BC03","BC03LN","BC04","BC05","BC06","BC07LN",
	"BC08","BC09","BC10")

rownames(count) <- c("Bcell","Tcell","macrophage")
count_cor <- cor(t(count), method="spearman")
round(count_cor, digits=2)

## plot heatmap of count
library("gplots")
source("/home05/liz86/Steffi/Functions/heatmap.3_v2.R")
mydist_e = function(c) {dist(c,method="euclidian")}
myclust_w = function(c) {hclust(c,method="ward.D")}
mydist_c = function(c) as.dist((1-cor(t(c),method="spearman"))/2)

pdf(file="Results/Heatmap_count_cor_chuang_scrna.pdf",
	width=6, height=6)
heatmap.2(count_cor, hclustfun=myclust_w, distfun=mydist_e, 
	na.rm = TRUE, scale="none", dendrogram="both", 
	margins=c(9,9), Rowv=FALSE, Colv=FALSE,
	symbreaks=FALSE, key=TRUE, keysize=0.6, symkey=FALSE, 
	density.info="none", trace="none", labCol=c("B_cell","T_cell","Macrophage"), 
	labRow=c("B_cell","T_cell","Macrophage"), col=greenred(75), 
	KeyValueName="spearman",side.height.fraction=0.1,
	lwid=c(0.2, 0.5),lhei=c(0.2, 0.5), cexRow=1, cexCol=1)
dev.off()


### Davolin AVE
corr_davolin_ave <- matrix(NA, nrow(davolin_ave), 3)
for(i in 1:nrow(davolin_ave)){
	for(j in 1:3){
		corr_davolin_ave[i, j] <- cor(davolin_ave[i,-12], count[j,], method="spearman")
	}
}
rownames(corr_davolin_ave) <- rownames(davolin_ave)
colnames(corr_davolin_ave) <- c("B_cell_count","T_cell_count","Macrophage_count")
round(corr_davolin_ave, digits=2)
## plot heatmap
library("gplots")
source("/home05/liz86/Steffi/Functions/heatmap.3_v2.R")
mydist_e = function(c) {dist(c,method="euclidian")}
myclust_w = function(c) {hclust(c,method="ward.D")}
mydist_c = function(c) as.dist((1-cor(t(c),method="spearman"))/2)

pdf(file="Results/Heatmap_davolin_ave_cor_chuang_scrna.pdf",
	width=5.5, height=4.5)
heatmap.3(corr_davolin_ave, hclustfun=myclust_w, distfun=mydist_e, 
	na.rm = TRUE, scale="none", dendrogram="both", 
	margins=c(9,9), Rowv=FALSE, Colv=FALSE,
	symbreaks=FALSE, key=TRUE, keysize=0.6, symkey=FALSE, 
	density.info="none", trace="none", labCol=colnames(corr_davolin_ave), 
	labRow=rownames(corr_davolin_ave), col=greenred(75), 
	KeyValueName="spearman",side.height.fraction=0.1,
	lwid=c(0.25, 0.6),lhei=c(0.2, 0.6), cexRow=1, cexCol=1)
dev.off()


### Davolin GSVA
corr_davolin_gsva <- matrix(NA, nrow(davolin_gsva), 3)
for(i in 1:nrow(davolin_gsva)){
	for(j in 1:3){
		corr_davolin_gsva[i, j] <- cor(davolin_gsva[i,-12], count[j,], method="spearman")
	}
}
rownames(corr_davolin_gsva) <- rownames(davolin_gsva)
colnames(corr_davolin_gsva) <- c("B_cell_count","T_cell_count","Macrophage_count")
round(corr_davolin_gsva, digits=2)

## plot heatmap
library("gplots")
source("/home05/liz86/Steffi/Functions/heatmap.3_v2.R")
mydist_e = function(c) {dist(c,method="euclidian")}
myclust_w = function(c) {hclust(c,method="ward.D")}
mydist_c = function(c) as.dist((1-cor(t(c),method="spearman"))/2)

pdf(file="Results/Heatmap_davolin_gsva_cor_chuang_scrna.pdf",
	width=5.5, height=4.5)
heatmap.3(corr_davolin_gsva, hclustfun=myclust_w, distfun=mydist_e, 
	na.rm = TRUE, scale="none", dendrogram="both", 
	margins=c(9,9), Rowv=FALSE, Colv=FALSE,
	symbreaks=FALSE, key=TRUE, keysize=0.6, symkey=FALSE, 
	density.info="none", trace="none", labCol=colnames(corr_davolin_gsva), 
	labRow=rownames(corr_davolin_gsva), col=greenred(75), 
	KeyValueName="spearman",side.height.fraction=0.1,
	lwid=c(0.25, 0.6),lhei=c(0.2, 0.5), cexRow=1, cexCol=1)
dev.off()

### Tamborero AVE
corr_tamborero_ave <- matrix(NA, nrow(tamborero_ave), 3)
for(i in 1:nrow(tamborero_ave)){
	for(j in 1:3){
		corr_tamborero_ave[i, j] <- cor(tamborero_ave[i,-12], count[j,], method="spearman")
	}
}
rownames(corr_tamborero_ave) <- rownames(tamborero_ave)
colnames(corr_tamborero_ave) <- c("B_cell_count","T_cell_count","Macrophage_count")
round(corr_tamborero_ave, digits=2)

## plot heatmap
library("gplots")
source("/home05/liz86/Steffi/Functions/heatmap.3_v2.R")
mydist_e = function(c) {dist(c,method="euclidian")}
myclust_w = function(c) {hclust(c,method="ward.D")}
mydist_c = function(c) as.dist((1-cor(t(c),method="spearman"))/2)

pdf(file="Results/Heatmap_tamborero_ave_cor_chuang_scrna.pdf",
	width=5.5, height=4.5)
heatmap.3(corr_tamborero_ave, hclustfun=myclust_w, distfun=mydist_e, 
	na.rm = TRUE, scale="none", dendrogram="both", 
	margins=c(9,9), Rowv=FALSE, Colv=FALSE,
	symbreaks=FALSE, key=TRUE, keysize=0.4, symkey=FALSE, 
	density.info="none", trace="none", labCol=colnames(corr_tamborero_ave), 
	labRow=rownames(corr_tamborero_ave), col=greenred(75), 
	KeyValueName="spearman",side.height.fraction=0.1,
	lwid=c(0.25, 0.6),lhei=c(0.15, 0.45), cexRow=1, cexCol=1)
dev.off()

### Tamborero GSVA
corr_tamborero_gsva <- matrix(NA, nrow(tamborero_gsva), 3)
for(i in 1:nrow(tamborero_gsva)){
	for(j in 1:3){
		corr_tamborero_gsva[i, j] <- cor(tamborero_gsva[i,-12], count[j,], method="spearman")
	}
}
rownames(corr_tamborero_gsva) <- rownames(tamborero_gsva)
colnames(corr_tamborero_gsva) <- c("B_cell_count","T_cell_count","Macrophage_count")
round(corr_tamborero_gsva, digits=2)

## plot heatmap
library("gplots")
source("/home05/liz86/Steffi/Functions/heatmap.3_v2.R")
mydist_e = function(c) {dist(c,method="euclidian")}
myclust_w = function(c) {hclust(c,method="ward.D")}
mydist_c = function(c) as.dist((1-cor(t(c),method="spearman"))/2)

pdf(file="Results/Heatmap_tamborero_gsva_cor_chuang_scrna.pdf",
	width=5.5, height=4.5)
heatmap.3(corr_tamborero_gsva, hclustfun=myclust_w, distfun=mydist_e, 
	na.rm = TRUE, scale="none", dendrogram="both", 
	margins=c(9,9), Rowv=FALSE, Colv=FALSE,
	symbreaks=FALSE, key=TRUE, keysize=0.4, symkey=FALSE, 
	density.info="none", trace="none", labCol=colnames(corr_tamborero_gsva), 
	labRow=rownames(corr_tamborero_gsva), col=greenred(75), 
	KeyValueName="spearman",side.height.fraction=0.1,
	lwid=c(0.25, 0.6),lhei=c(0.15, 0.45), cexRow=1, cexCol=1)
dev.off()

### ciber_abs
corr_ciber_abs <- matrix(NA, nrow(ciber_abs), 3)
for(i in 1:nrow(ciber_abs)){
	for(j in 1:3){
		corr_ciber_abs[i, j] <- cor(as.numeric(ciber_abs[i,-12]), 
			count[j,], method="spearman")
	}
}
rownames(corr_ciber_abs) <- rownames(ciber_abs)
colnames(corr_ciber_abs) <- c("B_cell_count","T_cell_count","Macrophage_count")
round(corr_ciber_abs, digits=2)

## plot heatmap
library("gplots")
source("/home05/liz86/Steffi/Functions/heatmap.3_v2.R")
mydist_e = function(c) {dist(c,method="euclidian")}
myclust_w = function(c) {hclust(c,method="ward.D")}
mydist_c = function(c) as.dist((1-cor(t(c),method="spearman"))/2)

pdf(file="Results/Heatmap_ciber_abs_cor_chuang_scrna.pdf",
	width=5.5, height=6)
heatmap.3(corr_ciber_abs, hclustfun=myclust_w, distfun=mydist_e, 
	na.rm = TRUE, scale="none", dendrogram="both", 
	margins=c(9,12), Rowv=FALSE, Colv=FALSE,
	symbreaks=FALSE, key=TRUE, keysize=0.4, symkey=FALSE, 
	density.info="none", trace="none", labCol=colnames(corr_ciber_abs), 
	labRow=rownames(corr_ciber_abs), col=greenred(75), 
	KeyValueName="spearman",side.height.fraction=0.1,
	lwid=c(0.2, 0.5),lhei=c(0.15, 0.6), cexRow=1, cexCol=1)
dev.off()

### ciber_res
count_rel <- sweep(count, 2, colSums(count), "/")
rownames(count_rel) <- paste(rownames(count_rel), "rel",sep="-")
count_rel[is.na(count_rel)] <- 0
corr_ciber_rel <- matrix(NA, nrow(ciber_rel), 3)
for(i in 1:nrow(ciber_rel)){
	for(j in 1:3){
		corr_ciber_rel[i, j] <- cor(as.numeric(ciber_rel[i,-12]), 
			count_rel[j,], method="spearman")
	}
}
colnames(corr_ciber_rel) <- rownames(count_rel)
rownames(corr_ciber_rel) <- rownames(ciber_rel)
round(corr_ciber_rel, digits=2)

## plot heatmap
library("gplots")
source("/home05/liz86/Steffi/Functions/heatmap.3_v2.R")
mydist_e = function(c) {dist(c,method="euclidian")}
myclust_w = function(c) {hclust(c,method="ward.D")}
mydist_c = function(c) as.dist((1-cor(t(c),method="spearman"))/2)

pdf(file="Results/Heatmap_ciber_rel_cor_chuang_scrna.pdf",
	width=5.5, height=6)
heatmap.3(corr_ciber_rel, hclustfun=myclust_w, distfun=mydist_e, 
	na.rm = TRUE, scale="none", dendrogram="both", 
	margins=c(9,12), Rowv=FALSE, Colv=FALSE,
	symbreaks=FALSE, key=TRUE, keysize=0.4, symkey=FALSE, 
	density.info="none", trace="none", labCol=colnames(corr_ciber_rel), 
	labRow=rownames(corr_ciber_rel), col=greenred(75), 
	KeyValueName="spearman",side.height.fraction=0.1,
	lwid=c(0.2, 0.5),lhei=c(0.15, 0.6), cexRow=1, cexCol=1)
dev.off()






