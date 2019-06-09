rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")
# R 3.4.0
library(limma) #3.32.6
library(GSVA)  # 1.26.0

#### load data
load("Data_v2/scRNA_Chung_tpm_noFiltering.RData")
tpm <- tpm_noFiltering
rm(tpm_noFiltering)
log2tpm <- log2(tpm+1)

## pick type
cell_names <- sapply(1:ncol(log2tpm), function(x) 
	strsplit(colnames(log2tpm)[x], split="_")[[1]][1])
type <- sapply(1:ncol(log2tpm), function(x) 
	strsplit(colnames(log2tpm)[x], split="_")[[1]][2])
type[type=="Tumor" | type=="Pooled"]
cor(log2tpm[,"BC01_Tumor"], log2tpm[,"BC01_Pooled"]) #0.842

## only selected pooled samples
tpm <- tpm[, type=="Pooled"]
dim(tpm)  #55860    12
log2tpm <- log2tpm[, type=="Pooled"]
dim(log2tpm) #55860    12

################################
### compute davoli signatures
################################
load("/net/wong05/home/liz86/Steffi/Kevin_IDC_ILC_DE/RawData/Immune_cell_signatures_Davoli2016.RData")

## check if any genes are not in the dataset
sapply(1:length(immune_cell_signatures_davoli), function(x) 
	immune_cell_signatures_davoli[[x]][!(immune_cell_signatures_davoli[[x]] 
		%in% rownames(log2tpm))])

fit <- gsva(log2tpm, immune_cell_signatures_davoli, rnaseq=FALSE)
gsva_score_davoli <- fit
save(gsva_score_davoli, file="Results_v2/gsva_score_davoli_log2tpm_Chung_pooled.RData")

################################
# Tamborero signature
###############################
load("/net/wong05/home/liz86/Steffi/Kevin_IDC_ILC_DE/RawData/Immune_cell_signatures_Tamborero2017.RData")
sapply(1:length(immune_cell_signatures_tamborero2017), function(x) 
	length(immune_cell_signatures_tamborero2017[[x]]))
## check if any genes are not in the dataset
sapply(1:length(immune_cell_signatures_tamborero2017), function(x) 
	immune_cell_signatures_tamborero2017[[x]][!
	(immune_cell_signatures_tamborero2017[[x]] %in% 
		rownames(log2tpm))])

# To alias
for(i in 1:length(immune_cell_signatures_tamborero2017)){
	out_gene <- immune_cell_signatures_tamborero2017[[i]][!
	(immune_cell_signatures_tamborero2017[[i]] %in% 
		rownames(log2tpm))]
	if(length(out_gene)>0){
		symbol2 <- alias2SymbolTable(out_gene, species="Hs")
		symbol2[is.na(symbol2)] <- "NA2"
		symbol2_in <- symbol2 %in% rownames(log2tpm)
		out_gene_hit <- out_gene[symbol2_in]
		out_gene_hit_symbol <- symbol2[symbol2_in]
		immune_cell_signatures_tamborero2017[[i]][match(out_gene_hit,
			immune_cell_signatures_tamborero2017[[i]])] <- out_gene_hit_symbol
	}
}

immune_cell_signatures_tamborero2017[[5]] <- gsub("6-Mar", "MARCH6",
	immune_cell_signatures_tamborero2017[[5]])

sapply(1:length(immune_cell_signatures_tamborero2017), function(x) 
	immune_cell_signatures_tamborero2017[[x]][!(immune_cell_signatures_tamborero2017[[x]] %in% 
		rownames(log2tpm))])

## GSVA (use log2 cpm)
tamborero_genes_common <- lapply(1:length(immune_cell_signatures_tamborero2017), 
	function(x) intersect(immune_cell_signatures_tamborero2017[[x]], 
		rownames(log2tpm)))
names(tamborero_genes_common) <- names(immune_cell_signatures_tamborero2017)

fit <- gsva(log2tpm, tamborero_genes_common,rnaseq=FALSE)
gsva_score_tamborero <- fit
save(gsva_score_tamborero, file="Results_v2/gsva_score_tamborero_tpm_Chung_pooled.RData")

#################################
# Prepare for TIMER
#################################
tpm_timer <- cbind(rownames(tpm), tpm)
write.csv(tpm_timer, file="Data_v2/tpm_timer_singlecell_chuang.csv", 
	quote=FALSE, row.names=FALSE) # file format doesn't work for cibersort

# load results
timer_raw <- read.csv("Results_v2/TIMER_res_singleCell_Chuang.csv")
timer <- timer_raw[, -1]
rownames(timer) <- timer_raw[,1]
timer <- t(timer)
save(timer,file="Results_v2/TIMER_res_singleCell_Chuang.RData")

#################################
# Prepare for CIBERSORT
#################################
tpm_for_cibersort <- cbind(rownames(tpm), tpm)
write.table(tpm_for_cibersort, file="Data_v2/tpm_for_cibersort_singlecell_chuang.txt",
	quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

# load results
ciber_raw <- read.csv("Results_v2/CIBERSORT_relative_singleCell_chuang.csv")
ciber_relative <- ciber_raw[,seq(2, 23)]
rownames(ciber_relative) <- ciber_raw[,1]
ciber_relative <- t(ciber_relative)
save(ciber_relative, file="Results_v2/cibersort_relative_singleCell_chuang.RData")

# load results
ciber_raw <- read.csv("Results_v2/CIBERSORT_absolute_singleCell_chuang.csv")
ciber_absolute <- ciber_raw[,seq(2, 23)]
rownames(ciber_absolute) <- ciber_raw[,1]
ciber_absolute <- t(ciber_absolute)
save(ciber_absolute, file="Results_v2/cibersort_absolute_singleCell_chuang.RData")





