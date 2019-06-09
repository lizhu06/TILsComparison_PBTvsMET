rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")

######################
# primary and mets 
######################
## count (primary, mets)
load("RawData_v2/PanMets.Paired.counts.salmon.0.8.2.Rda")
dim(new) #60619   105
count <- new
save(count, file="Data_v2/count.RData")

### log2 cpm file
load("RawData_v2/PanMets.Paired.log2.tmmNorm.cpm.salmon.0.8.2.Rda")
log2cpm <- tmp.2[,1:(ncol(tmp.2)-4)]
dim(log2cpm)  # 60619   105
log2cpm_mean <- apply(log2cpm,1,mean)
summary(log2cpm_mean)
save(log2cpm, file="Data_v2/log2cpm.RData")

### tpm file
load("RawData_v2/txi.salmon082.PanMets.GeneTPM.Rda")

tpm <- tpm[, colnames(count)]
dim(tpm) #60619   105
all(colnames(tpm) == colnames(count)) #60619   105
save(tpm, file="Data_v2/tpm.RData")

### sample annotation file
annot <- read.csv("RawData_v2/SampleAnnotation.Li.csv")
dim(annot) #108 9
load("Data_v2/all_clinical_info.RData")
all_clin_pair <- annot[match(all_clin[,"primary_id"], annot[,"ID"]), "pair"]
all_clin_mets_id <- sapply(1:length(all_clin_pair), function(x) annot[which(
  annot[,"pair"]==all_clin_pair[x] & annot[,"group"]=="MET"), "ID"])
all_clin$mets_id <- all_clin_mets_id

# check if clinical info consistent
annot_prim_ER <- annot[match(all_clin$primary_id, annot[,"ID"]), c("ID", "ER_status")]
which(annot_prim_ER[,2]!=all_clin[,"ER.Prim"])

annot_mets_ER <- annot[match(all_clin$mets_id, annot[,"ID"]), c("ID", "ER_status")]
which(annot_mets_ER[,2]!=all_clin[,"ER.MET"])
all_clin[34, ] # OM2 inconsistent

# add ER, PR, HR, HER2 to annot
id <- annot[,"ID"]
all(id %in% colnames(log2cpm))
dup_id <- id[duplicated(id)]  # "OP5" "GP2" "GP4" matched to two pairs
annot[which(annot[,"ID"] %in% dup_id),]
sample_annot <- annot

sample_annot$ER <- rep(NA, nrow(sample_annot))
sample_annot[which(sample_annot[,"group"]=="Primary"), "ER"] <- all_clin[match(
  sample_annot[which(sample_annot[,"group"]=="Primary"), "ID"], all_clin[,"primary_id"]), "ER.Prim"]
sample_annot[which(sample_annot[,"group"]=="MET"), "ER"] <- all_clin[match(
  sample_annot[which(sample_annot[,"group"]=="MET"), "ID"], all_clin[,"mets_id"]), "ER.MET"]

sample_annot$PR <- rep(NA, nrow(sample_annot))
sample_annot[which(sample_annot[,"group"]=="Primary"), "PR"] <- all_clin[match(
  sample_annot[which(sample_annot[,"group"]=="Primary"), "ID"], all_clin[,"primary_id"]), "PR.Prim"]
sample_annot[which(sample_annot[,"group"]=="MET"), "PR"] <- all_clin[match(
  sample_annot[which(sample_annot[,"group"]=="MET"), "ID"], all_clin[,"mets_id"]), "PR.MET"]

sample_annot$HER2 <- rep(NA, nrow(sample_annot))
sample_annot[which(sample_annot[,"group"]=="Primary"), "HER2"] <- all_clin[match(
  sample_annot[which(sample_annot[,"group"]=="Primary"), "ID"], all_clin[,"primary_id"]), "HER2.Prim"]
sample_annot[which(sample_annot[,"group"]=="MET"), "HER2"] <- all_clin[match(
  sample_annot[which(sample_annot[,"group"]=="MET"), "ID"], all_clin[,"mets_id"]), "HER2.MET"]

sample_annot$HR <- rep("Pos", nrow(sample_annot))
sample_annot$HR[sample_annot$ER=="Neg" & sample_annot$PR=="Neg"] <- "Neg"

save(sample_annot, file="Data_v2/sample_annot.RData")

### sample_pair
pairs <- unique(annot[,"pair"])  # determines the order
mets_site <- sapply(1:length(pairs), function(x) annot[which(
	annot[,"pair"]==pairs[x] & annot[,"group"] == "MET"), "site"])		
primary_id <- sapply(1:length(pairs), function(x) annot[which(
	annot[,"pair"]==pairs[x] & annot[,"group"] == "Primary"), "ID"])
mets_id <- sapply(1:length(pairs), function(x) annot[which(
	annot[,"pair"]==pairs[x] & annot[,"group"] == "MET"), "ID"])
ER_prim <- sapply(1:length(pairs), function(x) sample_annot[which(
  sample_annot[,"pair"]==pairs[x] & 
  sample_annot[,"group"] == "Primary"), "ER"])  	
ER_MET <- sapply(1:length(pairs), function(x) sample_annot[which(
  sample_annot[,"pair"]==pairs[x] & 
  sample_annot[,"group"] == "MET"), "ER"]) 

PR_prim <- sapply(1:length(pairs), function(x) sample_annot[which(
  sample_annot[,"pair"]==pairs[x] & 
  sample_annot[,"group"] == "Primary"), "PR"])    
PR_MET <- sapply(1:length(pairs), function(x) sample_annot[which(
  sample_annot[,"pair"]==pairs[x] & 
  sample_annot[,"group"] == "MET"), "PR"])  

HER2_prim <- sapply(1:length(pairs), function(x) sample_annot[which(
  sample_annot[,"pair"]==pairs[x] & 
  sample_annot[,"group"] == "Primary"), "HER2"])    
HER2_MET <- sapply(1:length(pairs), function(x) sample_annot[which(
  sample_annot[,"pair"]==pairs[x] & 
  sample_annot[,"group"] == "MET"), "HER2"])

HR_prim <- sapply(1:length(pairs), function(x) sample_annot[which(
  sample_annot[,"pair"]==pairs[x] & 
  sample_annot[,"group"] == "Primary"), "HR"])    
HR_MET <- sapply(1:length(pairs), function(x) sample_annot[which(
  sample_annot[,"pair"]==pairs[x] & 
  sample_annot[,"group"] == "MET"), "HR"])      

sample_pair_info <- data.frame(pairs=pairs, mets_site=mets_site, 
	primary_id=primary_id, mets_id=mets_id, 
  ER_prim=ER_prim, ER_MET=ER_MET, PR_prim=PR_prim, PR_MET=PR_MET,
  HER2_prim=HER2_prim, HER2_MET=HER2_MET, HR_prim=HR_prim, HR_MET=HR_MET)
save(sample_pair_info, file="Data_v2/sample_pair_info.RData")

### gene annotation using biomart
load("Data_v2/gene_annot_biomart.RData")

## Normal ids
normal_ov <- c("0089-56L","0090-56L","0091-56L","0092-56L")
normal_breast <- c("0036-50L","0037-50L","0038-50L","0039-50L", "0040-50L","ON5")
save(normal_ov, file="Data_v2/normal_ov_id.RData")
save(normal_breast, file="Data_v2/normal_breast_id.RData")

load("RawData_v2/txi.salmon082.PanMets.GeneCounts.with_normal.rda")
count_normal <- cts[, c(normal_ov, normal_breast)]
dim(count_normal) #60619    10
save(count_normal, file="Data_v2/count_normal.RData")

load("RawData_v2/txi.salmon082.PanMets.GeneTPM.Rda")
tpm_normal <- tpm[, c(normal_ov, normal_breast)]
dim(tpm_normal)
save(tpm_normal, file="Data_v2/tpm_normal.RData")


######################
# normal tissue data
######################
load("Data/normal_median_log2tpm.RData")

######################
# local recurrence
######################
## count
load("RawData/localRecurrences.salmon.cts.Rda")
count_local <- localRecurrences.salmon.cts
save(count_local, file="Data/count_local.RData")

## log2cpm
load("RawData/log2cpm_localRecurrences_TMMnormalized_by_Li.RData") # local recur dara
log2cpm_local <- log2.tmmNorm.cpm
save(log2cpm_local, file="Data/log2cpm_local.RData")

load("Data/local_gene_annot_biomart.RData")

