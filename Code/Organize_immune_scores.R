rm(list=ls())
options(stringsAsFactors = FALSE)
library(survival)
library(rms)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")

#### load data
load("Data_v2/log2tpm.RData")
load("Data_v2/sample_annot.RData")
load("Data_v2/sample_pair_info.RData")
#load("Data/gene_annot_biomart.RData")

## load immune scores
load("Results_v2/gsva_score_davoli_log2tpm.RData")
load("Results_v2/gsva_score_tamborero_log2tpm.RData")
rownames(gsva_score_tamborero_log2tpm) <- gsub(" ", "_", 
	rownames(gsva_score_tamborero_log2tpm))
rownames(gsva_score_tamborero_log2tpm) <- gsub("\\+", "_", 
	rownames(gsva_score_tamborero_log2tpm))
rownames(gsva_score_tamborero_log2tpm) <- gsub(
  " ", ".", rownames(gsva_score_tamborero_log2tpm))
load("Results_v2/cibersort_relative.RData")
load("Results_v2/cibersort_absolute.RData")
load("Results_v2/TIMER_res_pri_mets.RData")
rownames(timer) <- gsub(
  "_", ".", rownames(timer))
load("Results_v2/log2tpm_estimate_purity_score.RData")

id <- colnames(gsva_score_davoli_log2tpm)
exc_id <- c("7M_RCS","7P_RCS")
id <- id[!(id %in% exc_id)]
all(colnames(gsva_score_davoli_log2tpm) == colnames(gsva_score_tamborero_log2tpm))
all(colnames(gsva_score_davoli_log2tpm) %in% colnames(ciber_absolute))
all(colnames(estimate_purity_score) == colnames(gsva_score_davoli_log2tpm))

all_res <- list(t(estimate_purity_score[, id]), 
	t(gsva_score_davoli_log2tpm[,id]), t(gsva_score_tamborero_log2tpm[,id]),
	t(ciber_relative[, id]), 
	t(ciber_absolute[, id]), 
	t(timer[, id]))
names(all_res) <- c("estimate", "davoli", "tambo",
	"ciber_rel", "ciber_abs", "timer")

sapply(all_res, dim)
save(all_res, file="Results_v2/immune_res.RData")

# exclude outliers, file from Nolan
final_id <- read.csv("Data/20170223_FINAL_Pitt_BrainMet-clinicalData.csv")
final_id <- final_id[1:21,]
exc_id <- c("7M_RCS","7P_RCS")

sample_pair_info <- sample_pair_info[which(
	sample_pair_info[,"mets_id"] != "7M_RCS"), ]
sample_annot <- sample_annot[which(!(
	sample_annot[,"ID"] %in% c("7M_RCS","7P_RCS"))), ]

# take average of mets matched to the same primary
ave_dup <- function(score_matrix, sample_pair_info){
	dup_id <- sample_pair_info[duplicated(
		sample_pair_info[,"primary_id"]), "primary_id"]
	#sample_pair_info[which(sample_pair_info[,"primary_id"] %in% dup_id), ]
	score_matrix2 <- score_matrix
	for(i in 1:length(dup_id)){
		mets_id <- sample_pair_info[which(sample_pair_info[,
			"primary_id"] == dup_id[i]), "mets_id"]
		ave_data <- apply(score_matrix[mets_id, ], 2, mean)
		score_matrix2[mets_id[1], ] <- 
			score_matrix2[mets_id[2], ] <- ave_data
	}
	return(score_matrix2)
}

all_res_aveDup <- lapply(1:length(all_res), 
	function(x) ave_dup(all_res[[x]], sample_pair_info))
names(all_res_aveDup) <- names(all_res)
immune_res_aveDup <- all_res_aveDup
save(immune_res_aveDup, file="Results_v2/immune_res_aveDup.RData")

