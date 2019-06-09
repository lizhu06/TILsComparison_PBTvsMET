rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")
# R 3.4.0
library(limma) #3.32.6
library(GSVA)  # 1.24.2

#### load data
data <- read.table("RawData/GSE75688_scRAN/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt", 
	header=TRUE)
data2 <- data.matrix(data[, -c(1,2,3)])
gene_annot <- data[,c(1,2,3)]

## sum up isoform
uni_gene <- unique(gene_annot[, "gene_name"])
length(uni_gene) # 55860

dup_genes <- gene_annot[which(duplicated(gene_annot[,"gene_name"])),"gene_name"]
uni_dup_genes <- unique(dup_genes)
not_dup_gene_names <- gene_annot[which(!(gene_annot[,"gene_name"] %in% uni_dup_genes)), 
	"gene_name"]
length(not_dup_gene_names) #55571
length(uni_dup_genes) # 289

sum_iso <- function(gene_name){
	sub_data <- data2[which(gene_annot[,"gene_name"] == gene_name), ]
	return(apply(sub_data, 2, sum))
}

dup_gene_sum_iso <- sapply(1:length(uni_dup_genes), 
	function(x) sum_iso(uni_dup_genes[x]))
data_sum_iso <- rbind(data2[match(not_dup_gene_names, 
	gene_annot[,"gene_name"]), ], 
	t(dup_gene_sum_iso))
rownames(data_sum_iso) <- c(not_dup_gene_names, uni_dup_genes)
dim(data_sum_iso) #55860   563 

# keep a copy of data set without filtering
tpm_noFiltering <- data_sum_iso
save(tpm_noFiltering, file="Data_v2/scRNA_Chung_tpm_noFiltering.RData")
log2tpm_noFiltering <- log2(tpm_noFiltering+1)
save(log2tpm_noFiltering, file="Data_v2/scRNA_Chung_log2tpm_noFiltered.RData")

# TPM < 1 to TPM=0
data_sum_iso[which(data_sum_iso<1)] <- 0 

# remove genes expressed in <10% of all tumor groups
get_id <- function(sample_name){
	return(substr(sample_name, start=1, stop=nchar(sample_name)-3))	
}

id <- c(colnames(data2)[1:14], sapply(15:ncol(data2), function(x) 
	get_id(colnames(data2)[x])))

tb <- table(id)
tumor_id <- names(tb)[tb>1]
zero_per_group <- function(x){
	zeor_rate_per_group <- sapply(1:length(tumor_id), function(i) 
		mean(x[which(id==tumor_id[i])]==0))
	return(zeor_rate_per_group)
}
zero_rate <- apply(data_sum_iso, 1, zero_per_group)

data_filter <- data_sum_iso[apply(t(zero_rate)>=0.9, 1, sum)<14,]
dim(data_filter) #17802   563

tpm_filtered <- data_filter
save(tpm_filtered, file="Data_v2/scRNA_Chung_tpm_filtered.RData")

log2tpm_filtered <- log2(tpm_filtered+1)
save(log2tpm_filtered, file="Data_v2/scRNA_Chung_log2tpm_filtered.RData")

#### load phenotype data
load("Data_v2/scRNA_Chung_log2tpm_filtered.RData")
cell_type <- read.table("RawData_v2/single_cell_final_sample_information.txt", header=TRUE)
rownames(cell_type) <- cell_type[,"sample"]
cell_type <- cell_type[, -1]
rownames(cell_type)[!(rownames(cell_type) %in% colnames(log2tpm_filtered))]
save(cell_type, file="Data_v2/cell_type_chung.RData")











