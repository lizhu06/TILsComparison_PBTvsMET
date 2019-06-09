rm(list=ls())
options(stringsAsFactors = FALSE)
library(survival)
library(rms)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")
library(ggplot2)
library(gridExtra)

#### load data
load("Data_v2/sample_annot.RData")
load("Data_v2/sample_pair_info.RData")
load("Results_v2/immune_res_aveDup.RData")
load("Data_v2/all_clinical_info.RData")
all_clin[, "Met.Location"] <- gsub("bone", "Bone", all_clin[, "Met.Location"])
all_clin[, "Met.Location"] <- gsub("ovary", "Ovary", all_clin[, "Met.Location"])
age_cut_median <- all_clin[,"Dx.Age"] > median(all_clin[, "Dx.Age"])
age_cut_median2 <- rep(NA, length(age_cut_median))
age_cut_median2[age_cut_median==1] <- "age>=median"
age_cut_median2[age_cut_median==0] <- "age<median"
all_clin[,"Dx.Age"] <- age_cut_median2

age_cut_median <- all_clin[,"mets.age"] > median(all_clin[, "mets.age"], 
	na.rm=TRUE)
age_cut_median2 <- rep(NA, length(age_cut_median))
age_cut_median2[age_cut_median==1] <- "mets.age>=median"
age_cut_median2[age_cut_median==0] <- "mets.age<median"
all_clin[,"mets.age"] <- age_cut_median2

tissue <- sample_annot[match(rownames(immune_res_aveDup$ciber_rel), 
	sample_annot[,"ID"]), "site"]

## variables to test
var_for_test <- c("Dx.Age", "mets.age", 
	"Met.Location", "Race", 
	"Primary.Histology.ILCIDC", "Menopausal.Status", 
	"Pathological.Stage.Num", 
	"ER.Prim", "ER.MET", "PR.Prim", "PR.MET", "HER2.Prim", "HER2.MET",
	"HER2.Prim.PosNeg", "HER2.MET.PosNeg", 
	"HR.Prim", "HR.MET", "HR.HER2.Prim", "HR.HER2.MET",
	"PreEndocrine", "PreHER2", "PreChemo")
all(var_for_test %in% colnames(all_clin))

#all_clin[] <- lapply(all_clin, factor)

#for(i in 1:length(var_for_test)){
#	print(var_for_test[i])
#	print(levels(all_clin[,i]))
#}

## function to test clinical var with immune
test_with_immune_per_var <- function(clin_var_name, all_clin, 
	immune_var_name, score_matrix, 
	primary=TRUE, mets=FALSE, delta=FALSE){

	if(FALSE){
		clin_var_name <- "Met.Location"
		score_matrix <- immune_res_aveDup$estimate
		score_matrix_name <- "estimate"
		if(score_matrix_name == "estimate"){
			score_matrix <- score_matrix[, -c(1, 3,4), drop=FALSE]
		}
		primary <- TRUE
		continuous <- FALSE
		immune_var_name <- colnames(score_matrix)[1]
	}

	id <- rownames(all_clin)[!is.na(all_clin[,clin_var_name])]
	var <- all_clin[id, clin_var_name]

	if(primary==TRUE){
		immune <- score_matrix[id, immune_var_name] 
	}else if(mets==TRUE){
		id <- sample_pair_info[match(id, 
			sample_pair_info[,"primary_id"]), "mets_id"]
		immune <- score_matrix[id, immune_var_name]
	}else{
		primary_id <- id
		mets_id <- sample_pair_info[match(id, 
				sample_pair_info[,"primary_id"]), "mets_id"]
		immune <- (score_matrix[mets_id, immune_var_name] - 
			score_matrix[primary_id, immune_var_name])
	} 

	# calculate mean at each categories
	uni_var_cat <- unique(var)
	res <- matrix(NA, 1, length(uni_var_cat)*3+1)
	for(i in 1:length(uni_var_cat)){
		immune_cat <- immune[var==uni_var_cat[i]]
		res[1, (i-1)*3+1] <- median(immune_cat)
		res[1, ((i-1)*3+2):((i-1)*3+3)] <- quantile(immune_cat, prob=c(0.25, 0.75))
	}
	if(length(uni_var_cat)>2){
		res[1, ncol(res)] <- kruskal.test(immune, factor(var))$p.value
	}else{
		res[1, ncol(res)] <- wilcox.test(immune[var==uni_var_cat[1]], 
			immune[var==uni_var_cat[2]])$p.value
	}
	colnames(res) <- c(sapply(1:length(uni_var_cat), function(x) paste(uni_var_cat[x], 
		c("median", "25perc", "75perc"), sep="_")), "pval")	
	return(res)
}

#res <- test_with_immune_per_var("Dx.Age", all_clin, 
#	"ImmuneScore", immune_res_aveDup$estimate, 
#	primary=TRUE, mets=FALSE, delta=FALSE, posNeg=FALSE)

#res <- test_with_immune_per_var("Race", all_clin, 
#	"ImmuneScore", immune_res_aveDup$estimate, 
#	primary=TRUE, mets=FALSE, delta=FALSE, posNeg=FALSE)

test_with_immune <- function(var_for_test, 
	score_matrix, immune_var_name){

	if(FALSE){
		score_matrix <- immune_res_aveDup$estimate
		score_matrix_name <- "estimate"
		immune_var_name <- "ImmuneScore"
	}

	#if(score_matrix_name == "estimate"){
	#	score_matrix <- score_matrix[, -c(1, 3, 4), drop=FALSE]
	#}

	primary_vec <- c(TRUE, FALSE, FALSE)
	mets_vec <- c(FALSE, TRUE, FALSE)
	delta_vec <- c(FALSE, FALSE, TRUE)
	table_name_vec <- c("primary", "mets", "delta")
	note <- NULL
	res_matrix_all <- list()
	res <- list()

	for(j in 1:3){
		res[[j]] <- list()
		for(x in 1:length(var_for_test)){
			res[[j]][[x]] <- test_with_immune_per_var(clin_var_name=var_for_test[x],
				all_clin, immune_var_name, 
				score_matrix, 
				primary=primary_vec[j], mets=mets_vec[j], 
				delta=delta_vec[j])
		}
		names(res[[j]]) <- var_for_test
	}
	names(res) <- c("primary", "mets", "delta")
	return(res)
}


foo <- test_with_immune(var_for_test, immune_res_aveDup[[1]], "ImmuneScore")







