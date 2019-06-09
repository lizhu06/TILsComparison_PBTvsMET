rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")
# R 3.4.0
library(edgeR) # 3.20.4
library(GSVA)  # 1.26.0

#### load data
load("RawData_v2/txi.salmon082.PanMets.GeneTPM.Rda")
load("Data_v2/log2tpm_unique.RData")
log2tpm <- log2tpm_unique
tpm <- tpm[rownames(log2tpm_unique), ] # has more samples!
load("Data_v2/sample_annot.RData")
load("Data_v2/sample_pair_info.RData")
#load("Data/estimate_purity_score.RData")
load("Data_v2/gene_annot_biomart_unique.RData")
rownames(log2tpm) <- rownames(tpm) <- gene_annot_biomart_unique[,"external_gene_name_v2"]
dim(log2tpm) #55557   105

dup_pri_id <- sample_pair_info[which(duplicated(sample_pair_info[,"primary_id"])),"primary_id"]
dup_mets_id <- lapply(1:length(dup_pri_id), function(x) sample_pair_info[which(
	sample_pair_info[,"primary_id"]==dup_pri_id[x]), "mets_id"])

log2tpm_cor <- sapply(1:length(dup_mets_id), function(x) 
	cor(log2tpm[,dup_mets_id[[x]]], method="spearman")[1,2])
#0.9066897 0.7821558 0.6636428

colnames(log2tpm)[c(1,21)]
cor(log2tpm[,c(1,21)], method="spearman")


################################
### compute davoli signatures
################################
load("Results_v2/gsva_score_davoli_log2tpm.RData")
gsva_davoli_cor <- sapply(1:length(dup_mets_id), function(x) 
	cor(gsva_score_davoli_log2tpm[,dup_mets_id[[x]]], method="spearman")[1,2])
gsva_davoli_cor
#0.7212121 0.6000000 0.6000000

all_cor <- cor(gsva_score_davoli_log2tpm, 
	method="spearman")
all_cor[lower.tri(all_cor, diag = TRUE)] <- NA
summary(c(all_cor), na.rm=TRUE)

################################
# Tamborero signature
###############################
load("Results_v2/gsva_score_tamborero_log2tpm.RData")
gsva_tamborero_cor <- sapply(1:length(dup_mets_id), function(x) 
	cor(gsva_score_tamborero_log2tpm[,dup_mets_id[[x]]], 
		method="spearman")[1,2])
gsva_tamborero_cor
#0.91470588  0.03235294 -0.05000000
all_cor <- cor(gsva_score_tamborero_log2tpm, 
	method="spearman")
all_cor[lower.tri(all_cor, diag = TRUE)] <- NA
summary(c(all_cor), na.rm=TRUE)

#################################
# Prepare for TIMER
#################################
load("Results_v2/TIMER_res_pri_mets.RData")
timer_cor <- sapply(1:length(dup_mets_id), function(x) 
	cor(timer[,dup_mets_id[[x]]], method="spearman")[1,2])
timer_cor
#1.0000000 1.0000000 0.8285714

#################################
# Prepare for CIBERSORT
#################################
load("Results_v2/cibersort_relative.RData")
cibersort_cor <- sapply(1:length(dup_mets_id), function(x) 
	cor(ciber_relative[,dup_mets_id[[x]]], method="spearman")[1,2])
cibersort_cor




