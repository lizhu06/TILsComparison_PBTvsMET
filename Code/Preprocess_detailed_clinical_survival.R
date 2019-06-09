rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")
library(plyr)
library(survival)
library(rms)

load("Data_v2/sample_annot.RData")  # sample annot
load("Data_v2/sample_pair_info.RData") # pair information
exc_id <- c("7M_RCS","7P_RCS")
sample_pair_info <- sample_pair_info[which(sample_pair_info[,"mets_id"] != "7M_RCS"), ]
sample_annot <- sample_annot[which(!(sample_annot[,"ID"] %in% c("7M_RCS","7P_RCS"))), ]
dim(sample_annot) # 106 9
dim(sample_pair_info) # 53 4

## load survival (has more clinical information) 
bone <- get(load("Data/bone_mets_clinical.RData")) #11 pairs
dim(bone) # 11 31
brain <- get(load("Data/brain_mets_clinical.RData")) #21
dim(brain) # 21 55
gi_ov <- get(load("Data/OV_GI_clinical.RData")) #missed 3
dim(gi_ov[which(gi_ov[,"Study.Met.Location"]=="ovary"),]) # 12 39
dim(gi_ov[which(gi_ov[,"Study.Met.Location"]=="GI"),]) # 5 39

ov_primary_id <- sample_pair_info[which(sample_pair_info[,
	"mets_site"]=="Ovary"), "primary_id"] 
ov_primary_id[!(ov_primary_id %in% gi_ov[,"primary_id"])] #0032-50T no clinical info (*******)

## clinical each data set
bone$Post.Endocrine <- bone$Post.HER2 <- rep(NA, nrow(bone))
bone_v2 <- cbind(rep("bone", nrow(bone)), bone[,c("Case","Dx.Age","Ethnicity",
		"Histological.Subtype", "Menopausal.Status", 
	"Pathological.Stage",
	"ER.Prim", "PR.Prim", "HER2.Prim", 
	"ER.Met", "PR.Met", "HER2.Met", 
	"Endocrine.Tx", "HER2.Tx", "PreBoM.Chemotherapy",
	"Post.Endocrine", "Post.HER2", "PostBoM.Chemotherapy", 
	"DFS", "BMFS", "SAR", "OS", "Survival.Status")])
bone_v2[which(!(bone_v2[,"Histological.Subtype"] %in% c("ILC", "IDC"))), 
	"Histological.Subtype"] <- "Mixed"

brain_v2 <- brain[,c("Study.Met.Location", "Case", "Dx.Age", "Race", 
	"Primary.Histology", "Menopausal.Status.At.Dx", 
		"Pathological.Stage", 
		"ER.Prim", "PR.Prim", "HER2.Prim", 
		"ER.BrM", "PR.BrM", "HER2.BrM", 
		"Adj.Endocrine.Tx", "Adj.HER2.Tx", "Adj.Chemotherapy",
		"PostBrM.Endocrine.Tx", "PostBrM.HER2.Tx", "PostBrM.Chemotherapy",
		"DFS", "BMFS", "SPBM", "OS","Vital.Status")]
brain_v2[which(brain_v2[,"Primary.Histology"]=="Mucinous Ductal Carcinoma"),
	"Primary.Histology"] <- "IDC"
brain_v2[which(brain_v2[,"Primary.Histology"]=="Invasive carcinoma w/ ductal and lobular features"),
	"Primary.Histology"] <- "Mixed"
brain_v2[which(brain_v2[,"Menopausal.Status.At.Dx"]=="N/A"), 
"Menopausal.Status.At.Dx"] <- NA

gi_ov$Post.Endocrine <- gi_ov$Post.HER2 <- 
	gi_ov$Post.Chemo <- rep(NA, nrow(gi_ov))
gi_ov_v2 <- gi_ov[,c("Study.Met.Location", "primary_id", 
	"Dx.Age", "Race", "Primary.Histology",
	"Menopausal.Status.At.Dx..Pre.or.Post.",
		"Pathological.Stage", 
		"ER.Prim", "PR.Prim", "HER2.Prim", 
		"ER.MET", "PR.MET", "HER2.MET", 
		"Adj.Endocrine.Tx", "Adj.HER2.Tx", "Adj.Chemotherapy",
		"Post.Endocrine", "Post.HER2", "Post.Chemo",
		"DFS", "MFS", "SPM", "OS", "Vital.Status")]

# combine clinical data
colnames(bone_v2) <- colnames(brain_v2) <- colnames(gi_ov_v2) <-
	c("Met.Location", "primary_id", "Dx.Age", "Race", 
		"Primary.Histology", "Menopausal.Status",
		"Pathological.Stage", 
			"ER.Prim", "PR.Prim", "HER2.Prim", 
			"ER.MET", "PR.MET", "HER2.MET", 
			"PreEndocrine", "PreHER2", "PreChemo",
			"PostEndocrine", "PostHER2", "PostChemo",
			"DFS", "MFS", "SPM", "OS", "vital")
all_clin <- rbind(bone_v2, brain_v2, gi_ov_v2)
all_clin[which(is.na(all_clin[,"MFS"])), ]
all_clin[which(is.na(all_clin[,"SPM"])), ]

# if ER/PR/HER2 missing in mets, set the same as primary
all_clin[which(is.na(all_clin[,"ER.MET"])),"ER.MET"] <- 
	all_clin[which(is.na(all_clin[,"ER.MET"])),"ER.Prim"]
all_clin[which(is.na(all_clin[,"PR.MET"])),"PR.MET"] <- 
	all_clin[which(is.na(all_clin[,"PR.MET"])),"PR.Prim"]
all_clin[which(is.na(all_clin[,"HER2.MET"])),"HER2.MET"] <- 
	all_clin[which(is.na(all_clin[,"HER2.MET"])),"HER2.Prim"]

# add other variables
all_clin$mets.age <- all_clin$Dx.Age + all_clin$MFS
all_clin$ER.MET[all_clin[,"ER.MET"]=="wkPos"] <- "Pos"
all_clin$HER2.Prim.PosNeg <- all_clin$HER2.Prim
all_clin$HER2.Prim.PosNeg[all_clin$HER2.Prim.PosNeg=="Equi"] <- "Neg"
all_clin$HER2.MET.PosNeg <- all_clin$HER2.MET
all_clin$HER2.MET.PosNeg[all_clin$HER2.MET.PosNeg=="Equi"] <- "Neg"
all_clin$Primary.Histology.ILCIDC <- all_clin$Primary.Histology
all_clin$Primary.Histology.ILCIDC[which(!
	(all_clin$Primary.Histology.ILCIDC %in% c("ILC", "IDC")))] <- NA
#all_clin$Pathological.Stage.Num <- all_clin$Pathological.Stage
#all_clin$Pathological.Stage.Num[which(all_clin$Pathological.Stage
#	%in% c("3A", "IIIA", "IIIC"))] <- "III"
#all_clin$Pathological.Stage.Num[which(all_clin$Pathological.Stage
#	%in% c("IA", "IB", "T1cN1a", "T1cN1b3"))] <- "I"
#all_clin$Pathological.Stage.Num[which(all_clin$Pathological.Stage
#	%in% c("IIA", "IIB", "T2N0", "T2N1", "T2N1b3", "T2N3"))] <- "II"
all_clin$Pathological.Stage.Num <- all_clin$Pathological.Stage
all_clin$Pathological.Stage.Num[which(all_clin$Pathological.Stage
	%in% c("3A", "IIIA", "IIIC"))] <- "III"
all_clin$Pathological.Stage.Num[which(all_clin$Pathological.Stage
	%in% c("I", "IA", "IB", "T1cN1a", "T1cN1b3", "IIA", "IIB", "T2N0", "T2N1", "T2N1b3", "T2N3"))] <- "I/II"
all_clin$Pathological.Stage.Num <- as.factor(all_clin$Pathological.Stage.Num)
all_clin$Pathological.Stage.Num <- factor(all_clin$Pathological.Stage.Num,levels=c("I/II", "III", "IV"))
all_clin$ER.changed <- (all_clin$ER.Prim != all_clin$ER.MET)*1
all_clin$PR.changed <- (all_clin$PR.Prim != all_clin$PR.MET)*1
all_clin[which(all_clin$PR.changed==1), ]
all_clin$HER2.changed <- (all_clin$HER2.Prim != all_clin$HER2.MET)*1
#all_clin[which(all_clin$HER2.changed==1),]
rownames(all_clin) <- all_clin[,"primary_id"]
#all_clin[,c("Pathological.Stage", "Pathological.Stage.Num")]

## add HR subgroup of primary and mets
all_clin$HR.Prim <- rep("Pos", nrow(all_clin))
all_clin$HR.Prim[which(all_clin$ER.Prim=="Neg" & all_clin$PR.Prim=="Neg")] <- "Neg"
table(all_clin$HR.Prim)

all_clin$HR.MET <- rep("Pos", nrow(all_clin))
all_clin$HR.MET[which(all_clin$ER.MET=="Neg" & all_clin$PR.MET=="Neg")] <- "Neg"
table(all_clin$HR.MET)

## add HR and HER2 subgroup (primary)
all_clin$HR.HER2.Prim <- rep(NA, nrow(all_clin))
all_clin$HR.HER2.Prim[which(all_clin$HR.Prim=="Pos" & 
  all_clin$HER2.Prim=="Pos")] <- "HR+/HER2+"
all_clin$HR.HER2.Prim[which(all_clin$HR.Prim=="Pos" & 
  all_clin$HER2.Prim=="Neg")] <- "HR+/HER2-"
all_clin$HR.HER2.Prim[which(all_clin$HR.Prim=="Neg" & 
  all_clin$HER2.Prim=="Pos")] <- "HR-/HER2+"
all_clin$HR.HER2.Prim[which(all_clin$HR.Prim=="Neg" & 
  all_clin$HER2.Prim=="Neg")] <- "TNBC"
all_clin$HR.HER2.Prim[which(all_clin$HR.Prim=="Pos" & 
  all_clin$HER2.Prim=="Equi")] <- "HR+/HER2="
table(all_clin$HR.HER2.Prim)

## add HR and HER2 subgroup (MET)
all_clin$HR.HER2.MET <- rep(NA, nrow(all_clin))
all_clin$HR.HER2.MET[which(all_clin$HR.MET=="Pos" & 
  all_clin$HER2.MET=="Pos")] <- "HR+/HER2+"
all_clin$HR.HER2.MET[which(all_clin$HR.MET=="Pos" & 
  all_clin$HER2.MET=="Neg")] <- "HR+/HER2-"
all_clin$HR.HER2.MET[which(all_clin$HR.MET=="Neg" & 
  all_clin$HER2.MET=="Pos")] <- "HR-/HER2+"
all_clin$HR.HER2.MET[which(all_clin$HR.MET=="Neg" & 
  all_clin$HER2.MET=="Neg")] <- "TNBC"
all_clin$HR.HER2.MET[which(all_clin$HR.MET=="Pos" & 
  all_clin$HER2.MET=="Equi")] <- "HR+/HER2="
table(all_clin$HR.HER2.MET)

save(all_clin, file="Data_v2/all_clinical_info.RData")

# code dummy variable
all_clin_dummy <- data.frame(rep(NA, nrow(all_clin)))
all_clin_dummy$age.old <- ifelse(all_clin$Dx.Age > median(all_clin$Dx.Age), 1, 0)
all_clin_dummy$Met.Location.bone <- ifelse(all_clin$Met.Location=="bone", 1, 0)
all_clin_dummy$Met.Location.GI <- ifelse(all_clin$Met.Location=="GI", 1, 0)
all_clin_dummy$Met.Location.ovary <- ifelse(all_clin$Met.Location=="ovary", 1, 0)
all_clin_dummy$RaceBlack <- ifelse(all_clin$Race=="Black", 1, 0)
all_clin_dummy$ILC <- ifelse(all_clin$Primary.Histology.ILCIDC=="ILC", 1, 0)
all_clin_dummy$postMenopausal <- ifelse(all_clin$Menopausal.Status=="Post", 1, 0)
all_clin_dummy$stageIII <- ifelse(all_clin$Pathological.Stage.Num=="III", 1, 0)
all_clin_dummy$stageIV <- ifelse(all_clin$Pathological.Stage.Num=="IV", 1, 0)
all_clin_dummy$ERpos.Prim <- ifelse(all_clin$ER.Prim == "Pos", 1, 0)
all_clin_dummy$PRpos.Prim <- ifelse(all_clin$PR.Prim == "Pos", 1, 0)
all_clin_dummy$HER2pos.Prim <- ifelse(all_clin$HER2.Prim.PosNeg == "Pos", 1, 0)
all_clin_dummy$ERpos.MET <- ifelse(all_clin$ER.MET == "Pos", 1, 0)
all_clin_dummy$PRpos.MET <- ifelse(all_clin$PR.MET == "Pos", 1, 0)
all_clin_dummy$HER2pos.MET <- ifelse(all_clin$HER2.MET.PosNeg == "Pos", 1, 0)
all_clin_dummy$preEndo <- ifelse(all_clin$PreEndocrine == "Yes", 1, 0)
all_clin_dummy$preHer2 <- ifelse(all_clin$PreHER2 == "Yes", 1, 0)
all_clin_dummy$preChemo <- ifelse(all_clin$PreChemo == "Yes", 1, 0)
all_clin_dummy$postEndo <- ifelse(all_clin$PostEndocrine == "Yes", 1, 0)
all_clin_dummy$postHer2 <- ifelse(all_clin$PostHER2 == "Yes", 1, 0)
all_clin_dummy$postChemo <- ifelse(all_clin$PostChemo == "Yes", 1, 0)
all_clin_dummy$mets.age.old <- ifelse(all_clin$mets.age > median(all_clin$mets.age, na.rm=TRUE), 1, 0)
all_clin_dummy$ER.changed <- all_clin$ER.changed
all_clin_dummy$PR.changed <- all_clin$PR.changed
all_clin_dummy$HER2.changed <- all_clin$HER2.changed
all_clin_dummy <- data.matrix(all_clin_dummy[, -1])
rownames(all_clin_dummy) <- rownames(all_clin)
save(all_clin_dummy, file="Data_v2/all_clin_dummy.RData")

## create survival data
MFS <- Surv(time=all_clin[,"MFS"], event=rep(1, nrow(all_clin)))
MFS[which(all_clin[,"Met.Location"]=="Brain"), ]
SPM <- Surv(time=all_clin[,"SPM"], event=(all_clin[,"vital"]=="Dead")*1)
pri_id <- all_clin[,"primary_id"]
mets_id <- sample_pair_info[match(pri_id, sample_pair_info[,"primary_id"]), 
	"mets_id"]
rownames(SPM) <- rownames(MFS) <- pri_id
MFS <- MFS[!is.na(MFS[,1]), ]
SPM <- SPM[!is.na(SPM[,1]), ]
dim(MFS) # 47
dim(SPM) # 46
save(MFS, file="Data_v2/MFS.RData")
save(SPM, file="Data_v2/SPM.RData")





