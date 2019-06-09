rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")
library(lubridate)
library(survival)

clin <- read.csv("RawData_v2/TILpaper/3.23.17-DUKE&UNC - TIL clinical Data - preliminary analysis -   2.2.17_v2.csv")
dim(clin) # 145  31
length(unique(clin[,"subject_id"])) # 125

load("RawData_v2/TILpaper/common_id.RData")
sum(clin[,"subject_id"] %in% common_id) # 49

clin_matched <- clin[match(common_id, clin[,"subject_id"]), ]
dim(clin_matched) # 49 31

table(clin_matched[,"race"])
race <- rep(NA, nrow(clin_matched))
race[clin_matched[,"race"] == 0] <- "White"
race[clin_matched[,"race"] == 1] <- "African-American"
race[clin_matched[,"race"] == 555] <- "Other"
race[clin_matched[,"race"] == 999] <- "Unknown"
tb <- table(race)
prop.table(tb)

table(clin_matched[, "insitution"])
insitution <- clin_matched[,"insitution"]

table(clin_matched[, "bcdiagnosis"])
histology <- rep(NA, nrow(clin_matched))
histology[which(clin_matched[,"bcdiagnosis"] == 1)] <- "Ductal"
histology[which(clin_matched[,"bcdiagnosis"] == 3)] <- "Lobular"
table(histology)
prop.table(table(histology))

# also load til table (for sample date variables)
tilcount <- read.csv("RawData_v2/TILpaper/3.23.17-TILbiopsydata-preliminaryanalysis-2.2.17_v2.csv")
pbt_til <- tilcount[which(tilcount[,"Sample.Location"] == "breast primary"), ]
pbt_til_matched <- pbt_til[match(common_id, pbt_til[,"ID"]), ]
dim(pbt_til_matched) #49 12

# ER
er <- clin_matched[,"er_status"]
er[er==1] <- "Pos"
er[er==2] <- "Neg"
er[!is.na(clin_matched[, "path_ER"])]
er[!is.na(clin_matched[, "path_ER"])] <- clin_matched[!is.na(
  clin_matched[,"path_ER"]), "path_ER"]
er[er==1] <- "Pos"
er[er==0] <- "Neg"
table(er) # 26 negative, 23 positive

# PR
pr <- clin_matched[,"pr_status"]
pr[pr==1] <- "Pos"
pr[pr==2] <- "Neg"
pr[!is.na(clin_matched[, "path_PR"])]
pr[!is.na(clin_matched[, "path_PR"])] <- clin_matched[!is.na(
  clin_matched[,"path_PR"]), "path_PR"]
pr[pr==1] <- "Pos"
pr[pr==0] <- "Neg"
pr[is.na(pr)] <- "Unknown"
table(pr) # 28 negative, 20 positive, 1 unknown
table(as.character(er), as.character(pr)) 

# HR
HR <- rep("Neg", length(pr))
HR[er=="Pos" | pr=="Pos"] <- "Pos"
table(HR) # 24 neg, 25 pos
prop.table(table(HR))

# HER2
HER2 <- rep(NA, nrow(clin_matched))
HER2_status <- clin_matched[,"her2_status"]
HER2_fish <- clin_matched[,"her2fish"]
HER2[(HER2_status == 3) | (HER2_status == 2 & HER2_fish == 1)] <- "Pos"
HER2[which((HER2_status %in% c(0, 1))|(HER2_status == 2 & HER2_fish == 2) )] <- "Neg"
HER2[clin_matched[, "path_HER2"] != ""]
HER2[which(clin_matched[, "path_HER2"]=="positive")] <- "Pos"
HER2[which(clin_matched[, "path_HER2"]=="negative")] <- "Neg"
table(HER2)
prop.table(table(HER2))

# BCsubtype
BCsubtype <- rep(NA, nrow(clin_matched))
BCsubtype[which(HR=="Pos" & HER2=="Pos")] <- "HR+/HER2+"
BCsubtype[which(HR=="Pos" & HER2=="Neg")] <- "HR+/HER2-"
BCsubtype[which(HR=="Neg" & HER2=="Pos")] <- "HR-/HER2+"
BCsubtype[which(HR=="Neg" & HER2=="Neg")] <- "HR-/HER2-"
table(BCsubtype)
prop.table(table(BCsubtype))

# age
dob <- as.Date(clin_matched[, "pt_dob_enrl"], format='%m/%d/%Y')
#dob <- as.Date(dob)
diag_date <- as.Date(clin_matched[, "biopsy_dt_diag"], 
  format='%m/%d/%Y')
mets_date <- as.Date(clin_matched[, "bmdiagnosisdate"], 
  format='%m/%d/%Y')
pbt_sample_date <- as.Date(pbt_til_matched[, "Sample.Date"], 
  format='%m/%d/%Y')
last_date <- as.Date(clin_matched[, "last_date"], 
  format='%m/%d/%Y')
death_ind <- clin_matched[,"death_status"]

#age <- difftime(diag_date, dob, units='weeks')/52.25
age <- time_length(difftime(diag_date, dob), "years")
#age <- time_length(difftime(mets_date, dob), "years") # wrong?
age <- as.numeric(age)
age[!is.na(clin_matched[,"age"])]
age[!is.na(clin_matched[,"age"])] <- clin_matched[!is.na(clin_matched[,"age"]),"age"]
mean(age, na.rm=TRUE) # 47.5 (wrong 50.1)
summary(age)

# MFS
MFS_time <- time_length(difftime(mets_date, diag_date), "years")
#MFS <- time_length(difftime(mets_date, brain_mets_sample_date), "years")
MFS_time[!is.na(clin_matched[,"BMFS"])]
MFS_time[!is.na(clin_matched[,"BMFS"])] <- clin_matched[!is.na(
  clin_matched[,"BMFS"]),"BMFS"]/12

mean(MFS_time, na.rm=TRUE) # 4.66 (wrong 3.9)
summary(MFS_time*12)

# SPM
SPM_time <- time_length(difftime(last_date, mets_date), "years")
SPM_time[!is.na(clin_matched[,"SPBM"])]
SPM_time[!is.na(clin_matched[,"SPBM"])] <- clin_matched[!is.na(clin_matched[,"SPBM"]),"SPBM"]/12
summary(SPM_time*12)

SPM_OS_time <- time_length(difftime(last_date, mets_date), "years")
SPM_OS_time[!is.na(clin_matched[,"OS"])]
SPM_OS_time[!is.na(clin_matched[,"OS"])] <- clin_matched[!is.na(clin_matched[,"OS"]),"OS"]/12

# save results
BRM_TIL_clin <- data.frame(ID=common_id, age=age, MFS=MFS_time, race=race,
 insitution=insitution, histology=histology , HR=HR, HER2=HER2, BCsubtype=BCsubtype)
save(BRM_TIL_clin, file="RawData_v2/TILpaper/clinical.RData")

MFS <- Surv(time=MFS_time*12, event=rep(1, length(MFS_time)))
SPM <- Surv(time=SPM_time*12, event=death_ind)
SPM_OS <- Surv(time=SPM_OS_time*12, event=death_ind)
rownames(SPM) <- rownames(MFS) <- rownames(SPM_OS) <- common_id
save(MFS, file="RawData_v2/TILpaper/MFS.RData")
save(SPM, file="RawData_v2/TILpaper/SPM.RData")
save(SPM_OS, file="RawData_v2/TILpaper/SPM_OS.RData")





