rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")
library(tidyr)
library(reshape2)
library(ggplot2)
library(multcomp)

## load  scores
load("Results_v2/immune_res_aveDup.RData")

#all_res_norm <- lapply(1:length(all_res), function(x) sweep(all_res[[x]],
#  1, 1-estimate_purity_score["TumorPurity", id], "/"))
all_res_norm <- immune_res_aveDup


names(all_res_norm) <- c("estimate", "davoli", "tambo",
  "ciber_rel", "timer")
sapply(all_res_norm, dim)
sapply(1:length(all_res_norm), function(x) rownames(all_res_norm[[x]]))

## load sample data
load("Data_v2/sample_pair_info.RData")
load("Data_v2/sample_annot.RData")
sample_pair_info <- sample_pair_info[
  which(sample_pair_info[,"mets_id"] != "7M_RCS"), ]
sample_annot <- sample_annot[which(
  !(sample_annot[,"ID"] %in% c("7M_RCS","7P_RCS"))), ]

mets_id <- sample_pair_info[,"mets_id"]
pri_id <- sample_pair_info[,"primary_id"]
pri_id[duplicated(pri_id)]
#delta_purity <- estimate_purity_score["TumorPurity", mets_id] - 
#  estimate_purity_score["TumorPurity", pri_id]
mets_site <- as.factor(sample_pair_info[,"mets_site"])
mets_site <- factor(mets_site, levels=c("Brain", "Ovary", "Bone", "GI"))
uni_mets_site <- levels(mets_site)

table(sample_pair_info$mets_site)
table(sample_pair_info$HR_prim)
table(sample_pair_info$HER2_prim)
table(sample_pair_info$HR_prim, sample_pair_info$HER2_prim)

sample_pair_info$HR_HER2 <- rep(NA, nrow(sample_pair_info))

sample_pair_info$HR_HER2[which(sample_pair_info$HR_prim=="Pos" & 
  sample_pair_info$HER2_prim=="Pos")] <- "HR+/HER2+"

sample_pair_info$HR_HER2[which(sample_pair_info$HR_prim=="Pos" & 
  sample_pair_info$HER2_prim=="Neg")] <- "HR+/HER2-"

sample_pair_info$HR_HER2[which(sample_pair_info$HR_prim=="Neg" & 
  sample_pair_info$HER2_prim=="Pos")] <- "HR-/HER2+"

sample_pair_info$HR_HER2[which(sample_pair_info$HR_prim=="Neg" & 
  sample_pair_info$HER2_prim=="Neg")] <- "TNBC"

sample_pair_info$HR_HER2[which(sample_pair_info$HR_prim=="Pos" & 
  sample_pair_info$HER2_prim=="Equi")] <- "HR+/HER2="

HR_HER2 <- sample_pair_info$HR_HER2
HR_HER2 <- factor(HR_HER2, levels=c("HR+/HER2+", "HR+/HER2-", "HR-/HER2+", 
  "TNBC", "HR+/HER2="))
table(HR_HER2)
#HRpos_HER2pos  HRpos_HER2neg  HRneg_HER2pos      tripleNeg HRpos_HER2equi
#            8             25              3             10              4
uni_HR_HER2 <- levels(HR_HER2)



table(sample_pair_info[which(sample_pair_info$mets_site=="Brain"), "HR_HER2"])
#HRneg_HER2pos HRpos_HER2neg HRpos_HER2pos     tripleNeg
#            3             5             5             8
uni_HR_HER2_brain <- c("HR+/HER2+", "HR+/HER2-", "HR-/HER2+", 
  "TNBC")

## function to plot and test
source("R_20180418/function_test_immune_plot_delta.R")

###### Combine tissues ######
foo <- paired_test_one_box(all_res_norm$estimate, 
  "Results/test_raw_immune/test_immune_estimate", 
  width=2, height=3)

foo <- paired_test_one_box(all_res_norm$davoli, 
  "Results/test_raw_immune/test_immune_davoli", 
  width=5, height=3)

foo <- paired_test_one_box(all_res_norm$tambo, 
  "Results/test_raw_immune/test_immune_tambo", 
  width=8, height=3)

foo <- paired_test_one_box(all_res_norm$timer, 
  "Results/test_raw_immune/test_immune_timer", 
  width=2.5, height=2.5)

foo <- paired_test_one_box(all_res_norm$ciber_rel, 
  "Results/test_raw_immune/test_immune_ciber_relative", 
  width=8, height=3.5)

ciber_rel_add_ratio <- data.frame(all_res_norm$ciber_rel)
ciber_rel_add_ratio$M2vsM1 <- all_res_norm$ciber_rel[,"Macrophages.M2"] /(all_res_norm$ciber_rel[, "Macrophages.M1"]+0.0000000001)

foo <- paired_test_one_box(ciber_rel_add_ratio, 
  "Results/test_raw_immune/test_immune_ciber_relative_add_ratio", 
  width=8, height=3.5)

###### tissue specific ######
foo <- paired_test_tissue_sep_one_box(all_res_norm$estimate, 
  "Results/test_raw_immune/test_immune_tissue_sep_estimate", 
  width=3, height=3)

foo <- paired_test_tissue_sep_one_box(all_res_norm$davoli, 
  "Results/test_raw_immune/test_immune_tissue_sep_davoli", 
  width=5.5, height=3)

foo <- paired_test_tissue_sep_one_box(all_res_norm$tambo, 
  "Results/test_raw_immune/test_immune_tissue_sep_tambo", 
  width=6.5, height=3)

foo <- paired_test_tissue_sep_one_box(all_res_norm$timer, 
  "Results/test_raw_immune/test_immune_tissue_sep_timer", 
  width=3.5, height=3)

foo <- paired_test_tissue_sep_one_box(all_res_norm$ciber_rel, 
  "Results/test_raw_immune/test_immune_tissue_sep_ciber_relative", 
  width=10, height=3.5)

#### PR and HER2  ####
foo <- paired_test_HR_HER2(all_res_norm$estimate, 
  "Results/test_raw_immune/test_immune_HR_HER2_estimate", 
  width=3, height=3, onlyBrain=FALSE)

foo <- paired_test_HR_HER2(all_res_norm$davoli, 
  "Results/test_raw_immune/test_immune_HR_HER2_davoli", 
  width=8, height=3, onlyBrain=FALSE)

foo <- paired_test_HR_HER2(all_res_norm$tambo, 
  "Results/test_raw_immune/test_immune_HR_HER2_tambo", 
  width=8, height=3, onlyBrain=FALSE)

foo <- paired_test_HR_HER2(all_res_norm$timer, 
  "Results/test_raw_immune/test_immune_HR_HER2_timer", 
  width=4, height=3, onlyBrain=FALSE)

foo <- paired_test_HR_HER2(all_res_norm$ciber_rel, 
  "Results/test_raw_immune/test_immune_HR_HER2_ciber_rel", 
  width=12, height=3, onlyBrain=FALSE)

#### PR and HER2 only Brain
foo <- paired_test_HR_HER2(all_res_norm$estimate, 
  "Results/test_raw_immune/test_immune_HR_HER2_estimate", 
  width=3, height=3, onlyBrain=TRUE)

foo <- paired_test_HR_HER2_ESTIMATE(all_res_norm$estimate, 
  "Results/test_raw_immune/test_immune_HR_HER2_estimate", 
  width=5, height=3, onlyBrain=TRUE)

foo <- paired_test_HR_HER2(all_res_norm$davoli, 
  "Results/test_raw_immune/test_immune_HR_HER2_davoli", 
  width=8, height=3, onlyBrain=TRUE)

foo <- paired_test_HR_HER2(all_res_norm$tambo, 
  "Results/test_raw_immune/test_immune_HR_HER2_tambo", 
  width=8, height=3, onlyBrain=TRUE)

foo <- paired_test_HR_HER2(all_res_norm$timer, 
  "Results/test_raw_immune/test_immune_HR_HER2_timer", 
  width=4, height=3, onlyBrain=TRUE)

foo <- paired_test_HR_HER2(all_res_norm$ciber_rel, 
  "Results/test_raw_immune/test_immune_HR_HER2_ciber_rel", 
  width=10, height=3, onlyBrain=TRUE)



