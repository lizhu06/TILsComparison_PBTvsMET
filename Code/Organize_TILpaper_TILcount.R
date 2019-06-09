rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")

tilcount <- read.csv("RawData_v2/TILpaper/3.23.17-TILbiopsydata-preliminaryanalysis-2.2.17.csv")
dim(tilcount) # 201 12

table(tilcount[,"Sample.Location"]) # brain mets 62, breast primary 114
length(unique(tilcount[, "Sample_."])) # 200
tilcount[which(tilcount[,"Sample_."] %in% 
  tilcount[duplicated(tilcount[, "Sample_."]), "Sample_."]), ]

length(unique(tilcount[,"ID"])) # 100

# match breast primary and brain mets

# 49 unique pairs of brain mets and breast primary
brain_mets_til <- tilcount[which(tilcount[,"Sample.Location"] == "brain met"), ] 
dim(brain_mets_til) # 62
length(unique(brain_mets_til[,"ID"])) # 54
primary_til <- tilcount[which(tilcount[,"Sample.Location"] == "breast primary"), ] 
dim(primary_til) # 114
length(unique(primary_til[,"ID"])) # 95
common_id <- intersect(unique(brain_mets_til[,"ID"]), unique(primary_til[,"ID"]))
length(common_id) # 49

sum(brain_mets_til[,"ID"] %in% common_id) #57
sum(primary_til[,"ID"] %in% common_id) # 51

#save(common_id, file="RawData_v2/TILpaper/common_id.RData")

# BRM: find out the samples with more than 1 slide
brain_mets_matched <- brain_mets_til[which(brain_mets_til[,"ID"] %in% common_id), ]
dim(brain_mets_matched) # 57
dup_brain_id <- brain_mets_matched[duplicated(brain_mets_matched[,"ID"]), "ID"]
length(dup_brain_id) # 8
brain_mets_matched[which(brain_mets_matched[,"ID"] %in% dup_brain_id), ]

# BRM: take the average for multiple to one match
for(i in 1:length(dup_brain_id)){
  brain_mets_matched[which(brain_mets_matched[,"ID"] == dup_brain_id[i]), "Final.TIL.Analysis"] <- 
   rep(mean(brain_mets_matched[
    which(brain_mets_matched[,"ID"] == dup_brain_id[i]), "Final.TIL.Analysis"]), 
    sum(brain_mets_matched[,"ID"] == dup_brain_id[i]))
}
brain_mets_matched[which(brain_mets_matched[,"ID"] %in% dup_brain_id), ]
brain_mets_matched_tilcount <- brain_mets_matched[match(common_id, brain_mets_matched[, "ID"]), ]
save(brain_mets_matched_tilcount, file="RawData_v2/TILpaper/brain_mets_matched_tilcount.RData")

# PBT: find out the samples with more than 1 slide
primary_matched <- primary_til[which(primary_til[,"ID"] %in% common_id), ]
dim(primary_matched) # 51
dup_primary_id <- primary_matched[duplicated(primary_matched[,"ID"]), "ID"]
length(dup_primary_id) # 2
primary_matched[which(primary_matched[,"ID"] %in% dup_primary_id), ]

# PBT take the average for multiple to one match
for(i in 1:length(dup_primary_id)){
  primary_matched[which(primary_matched[,"ID"] == dup_primary_id[i]), "Final.TIL.Analysis"] <- 
   rep(mean(primary_matched[
    which(primary_matched[,"ID"] == dup_primary_id[i]), "Final.TIL.Analysis"]), 
    sum(primary_matched[,"ID"] == dup_primary_id[i]))
}
primary_matched[which(primary_matched[,"ID"] %in% dup_primary_id), ]
primary_matched_tilcount <- primary_matched[match(common_id, primary_matched[, "ID"]), ]
save(primary_matched_tilcount, file="RawData_v2/TILpaper/primary_matched_tilcount.RData")





