rm(list=ls())
options(stringsAsFactors = FALSE)
library(survival)
library(rms)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")
library(ggplot2)
library(gridExtra)

#### load data
load("RawData_v2/TILpaper/brain_mets_matched_tilcount.RData")
load("RawData_v2/TILpaper/primary_matched_tilcount.RData")
load("RawData_v2/TILpaper/clinical.RData")

pbt_immune <- primary_matched_tilcount[,"Final.TIL.Analysis"]
brm_immune <- brain_mets_matched_tilcount[,"Final.TIL.Analysis"]
delta_immune <- brm_immune - pbt_immune
names(pbt_immune) <- names(brm_immune) <- names(delta_immune) <- brain_mets_matched_tilcount[,"ID"]

age_cut <- rep(NA, nrow(BRM_TIL_clin))
age_cut[BRM_TIL_clin[,"age"] >= 50] <- "age>=median"
age_cut[BRM_TIL_clin[,"age"] < 50] <- "age<median"
BRM_TIL_clin[,"age"] <- age_cut
all_clin <- BRM_TIL_clin
rownames(all_clin) <- all_clin[,"ID"]
## variables to test
var_for_test <- c("age", "race", "HR", "HER2")
all(var_for_test %in% colnames(all_clin))


primary_vec <- c(TRUE, FALSE, FALSE)
mets_vec <- c(FALSE, TRUE, FALSE)
delta_vec <- c(FALSE, FALSE, TRUE)
table_name_vec <- c("primary", "mets", "delta")
note <- NULL
res_matrix_all <- list()
res_list <- list()

for(j in 1:3){
	res_list[[j]] <- list()
	for(x in 1:length(var_for_test)){
		clin_var_name <- var_for_test[x]
		id <- all_clin[!is.na(all_clin[,clin_var_name]), "ID"]
		var <- all_clin[id, clin_var_name]

		primary <- primary_vec[j]
		mets <- mets_vec[j]
		delta <- delta_vec[j]

		if(primary==TRUE){
			immune <- pbt_immune[id]
		}else if(mets==TRUE){
			immune <- brm_immune[id]
		}else{
			immune <- delta_immune[id]
		} 

		# calculate mean at each categories
		uni_var_cat <- unique(var)
		res <- matrix(NA, 1, length(uni_var_cat)*3+1)
		for(i in 1:length(uni_var_cat)){
			immune_cat <- immune[var==uni_var_cat[i]]
			res[1, (i-1)*3+1] <- median(immune_cat)
			#res[1, ((i-1)*3+2):((i-1)*3+3)] <- quantile(immune_cat, prob=c(0.25, 0.75))
		}
		if(length(uni_var_cat)>2){
			res[1, ncol(res)] <- kruskal.test(immune, factor(var))$p.value
		}else{
			res[1, ncol(res)] <- wilcox.test(immune[var==uni_var_cat[1]], 
				immune[var==uni_var_cat[2]])$p.value
		}
		colnames(res) <- c(sapply(1:length(uni_var_cat), function(x) paste(uni_var_cat[x], 
			c("median", "25perc", "75perc"), sep="_")), "pval")	

		res_list[[j]][[x]] <- res
	}
	names(res_list[[j]]) <- var_for_test
}
names(res_list) <- c("primary", "mets", "delta")









