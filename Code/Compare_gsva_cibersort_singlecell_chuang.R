rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")
# R 3.4.0
library(limma) #3.32.6
library(GSVA)  # 1.24.2

## load  scores
davoli_gsva <- get(load("Results_v2/gsva_score_davoli_log2tpm_Chung_pooled.RData"))
tamborero_gsva <- get(load("Results_v2/gsva_score_tamborero_tpm_Chung_pooled.RData"))
ciber_abs <- get(load("Results_v2/cibersort_absolute_singleCell_chuang.RData"))
ciber_rel <- get(load("Results_v2/cibersort_relative_singleCell_chuang.RData"))
timer <- get(load("Results_v2/TIMER_res_singleCell_Chuang.RData"))

## calculate proportion
load("Data_v2/cell_type_chung.RData")
id <- sapply(1:nrow(cell_type), function(x) strsplit(rownames(cell_type)[x],
  split="_")[[1]][1])
type <- sapply(1:nrow(cell_type), function(x) strsplit(rownames(cell_type)[x],
  split="_")[[1]][2])
keep_index <- which(!(type %in% c("Pooled", "Re")))
cell_type_keep <- cell_type[keep_index, ]
id_keep <- id[!(type %in% c("Pooled", "Re"))]
type_keep <- type[!(type %in% c("Pooled", "Re"))]

uni_id_keep <- sort(unique(id_keep))
count <- matrix(NA, length(uni_id_keep), 3)
for(i in 1:length(uni_id_keep)){
  count[i,1] <- sum(cell_type_keep[id_keep==uni_id_keep[i],4]=="Bcell")
  count[i,2] <- sum(cell_type_keep[id_keep==uni_id_keep[i],4]=="Tcell")
  count[i,3] <- sum(cell_type_keep[id_keep==uni_id_keep[i],4]=="Myeloid")
}
rownames(count) <- uni_id_keep
colnames(count) <- c("Bcell","Tcell","Macrophages")

prop_totalImmune <- sweep(count, 1, apply(count,1, sum), "/")
prop_totalImmune[is.na(prop_totalImmune)]<-0

# remove BC07, which doesn't have RNAseq data
count <- count[which(rownames(count) != "BC07"), ]
prop_totalImmune <- prop_totalImmune[which(rownames(prop_totalImmune) 
  != "BC07"), ]

## organize GSVA and deconvolution results
all <- list(davoli_gsva, tamborero_gsva, ciber_rel, timer)
names(all) <- c("Davoli", "Tamborero", "CIBERSORT", "TIMER")
for(i in 1:length(all)){
  id <- sapply(1:ncol(all[[i]]), function(x) strsplit(colnames(all[[i]])[x],
    split="_")[[1]][1])
  print(all(id==rownames(count)))
  colnames(all[[i]]) <- id
}

cor_list <- list()
for(i in 1:length(all)){
  cor_list[[i]] <- matrix(NA, nrow(all[[i]]), 3)
  for(s in 1:nrow(all[[i]])){
    for(t in 1:3){
      if(i != 3){
        cor_list[[i]][s,t] <- cor(all[[i]][s,], count[,t], method="spearman")
      } else{
        cor_list[[i]][s,t] <- cor(all[[i]][s,], prop_totalImmune[,t], method="spearman")
      }
      
    }
  }
  rownames(cor_list[[i]]) <- rownames(all[[i]])
  colnames(cor_list[[i]]) <- colnames(count)
}
names(cor_list) <- names(all)

## plot heatmap
height_vec <- c(4.5, 5, 5, 4)
library("gplots")
source("/home05/liz86/Steffi/Functions/heatmap.3_v2.R")
mydist_e = function(c) {dist(c,method="euclidian")}
myclust_w = function(c) {hclust(c,method="ward.D")}
mydist_c = function(c) as.dist((1-cor(t(c),method="spearman"))/2)

i <- 1
pdf(file=paste0("Results_v2/Heatmap_corr_", names(cor_list)[i], 
  "_singlecell_chung.pdf"), width=4.5, height=4)
heatmap.3(cor_list[[i]], hclustfun=myclust_w, distfun=mydist_e, 
  na.rm = TRUE, scale="none", dendrogram="both", 
  margins=c(9,9), Rowv=FALSE, Colv=FALSE,
  symbreaks=FALSE, key=TRUE, keysize=0.6, symkey=FALSE, 
  density.info="none", trace="none", 
  labCol=colnames(cor_list[[i]]), 
  labRow=rownames(cor_list[[i]]), col=greenred(75), 
  KeyValueName="Spearman cor",side.height.fraction=0.1,
  lwid=c(0.25, 0.6),lhei=c(0.2, 0.5), cexRow=1, cexCol=1)
dev.off()

i <- 2
pdf(file=paste0("Results_v2/Heatmap_corr_", names(cor_list)[i], 
  "_singlecell_chung.pdf"), width=4.5, height=4.5)
heatmap.3(cor_list[[i]], hclustfun=myclust_w, distfun=mydist_e, 
  na.rm = TRUE, scale="none", dendrogram="both", 
  margins=c(9,9), Rowv=FALSE, Colv=FALSE,
  symbreaks=FALSE, key=TRUE, keysize=0.6, symkey=FALSE, 
  density.info="none", trace="none", 
  labCol=colnames(cor_list[[i]]), 
  labRow=rownames(cor_list[[i]]), col=greenred(75), 
  KeyValueName="Spearman cor",side.height.fraction=0.1,
  lwid=c(0.25, 0.6),lhei=c(0.2, 0.6), cexRow=1, cexCol=1)
dev.off()

i <- 3
pdf(file=paste0("Results_v2/Heatmap_corr_", names(cor_list)[i], 
  "_singlecell_chung.pdf"), width=4.5, height=5)
heatmap.3(cor_list[[i]], hclustfun=myclust_w, distfun=mydist_e, 
  na.rm = TRUE, scale="none", dendrogram="both", 
  margins=c(9,9), Rowv=FALSE, Colv=FALSE,
  symbreaks=FALSE, key=TRUE, keysize=0.6, symkey=FALSE, 
  density.info="none", trace="none", 
  labCol=colnames(cor_list[[i]]), 
  labRow=rownames(cor_list[[i]]), col=greenred(75), 
  KeyValueName="Spearman cor",side.height.fraction=0.1,
  lwid=c(0.25, 0.6),lhei=c(0.2, 0.65), cexRow=0.75, cexCol=1)
dev.off()

i <- 4
pdf(file=paste0("Results_v2/Heatmap_corr_", names(cor_list)[i], 
  "_singlecell_chung.pdf"), width=4.5, height=4)
heatmap.3(cor_list[[i]], hclustfun=myclust_w, distfun=mydist_e, 
  na.rm = TRUE, scale="none", dendrogram="both", 
  margins=c(9,9), Rowv=FALSE, Colv=FALSE,
  symbreaks=FALSE, key=TRUE, keysize=0.6, symkey=FALSE, 
  density.info="none", trace="none", 
  labCol=colnames(cor_list[[i]]), 
  labRow=rownames(cor_list[[i]]), col=greenred(75), 
  KeyValueName="Spearman cor",side.height.fraction=0.1,
  lwid=c(0.25, 0.6),lhei=c(0.2, 0.5), cexRow=1, cexCol=1)
dev.off()


