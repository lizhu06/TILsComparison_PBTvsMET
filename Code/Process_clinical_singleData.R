rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")

## load data
load("Data_v2/sample_pair_info.RData")

########################
## brain mets survival
########################
data <- read.csv("RawData/primary_and_mets_data/20170223_FINAL_Pitt_BrainMet-clinicalData.csv")
data <- data[c(1:21), ]
dim(data) #21 55

## change name
id <- data[,"Case"]

replace_name <- function(x){
	ins <- strsplit(x, split="_")[[1]][2]
	if(ins == "RCS"){
		seg <- strsplit(x, split="_")[[1]]
		name <- paste(seg[1], "P_", seg[2], sep="")
	}else if(ins == "Pitt"){
		seg <- strsplit(x, split="_")[[1]]
		name <- paste("BP", seg[1], sep="")
	}
	return(name)
}

new_id <- sapply(1:length(id), function(x) replace_name(id[x]))
all(new_id %in% sample_pair_info[,"primary_id"])
setdiff(sample_pair_info[which(sample_pair_info[,"mets_site"]=="Brain"),"primary_id"], new_id) #"7P_RCS"

data[,"Case"] <- new_id
save(data, file="Data_v2/brain_mets_clinical.RData")

########################
## Bone mets survival
########################
data <- read.csv("RawData/primary_and_mets_data/Bone_mets_survival.csv")
dim(data) #11 31
id <- data[,"Case"]
new_id <- paste(id,"P", sep="")
all(new_id %in% sample_pair_info[,"primary_id"])
data[,"Case"] <- new_id
save(data, file="Data_v2/bone_mets_clinical.RData")


########################
## OV and GI survival
########################
data <- read.csv("RawData/primary_and_mets_data/Ova.GI.PairedMet.ClinData.Li.csv")
pri_id <- data[,"ID"]
data <- data[which(pri_id!=""), ]
colnames(data)[2] <- "primary_id"
save(data, file="Data_v2/OV_GI_clinical.RData")





