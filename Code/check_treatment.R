rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")
load("Data_v2/all_clinical_info.RData")

tab2 <- function(x){
	x[is.na(x)] <- "missing"
	print(table(x))
}

tab2(all_clin$PreEndocrine)
tab2(all_clin$PreHER2)
tab2(all_clin$PreChemo)

tab2(all_clin$PostEndocrine)
tab2(all_clin$PostHER2)
tab2(all_clin$PostChemo)

colname <- c("Met.Location","PreEndocrine", "PreChemo", "PreHER2", 
	"PostEndocrine", "PostChemo", "PostHER2")
all_clin[,colname]

tab2_two <- function(x,y){
	x[is.na(x)] <- "missing"
	y[is.na(y)] <- "missing"
	print(table(x,y))
}

tab2_two(all_clin$PreEndocrine, all_clin$ER.Prim)
tab2_two(all_clin$PreHER2, all_clin$HER2.Prim)

all_clin[which(all_clin$PreEndocrine == "Yes" & all_clin$ER.Prim=="Neg"), ]
all_clin[which(all_clin$PreHER2 == "Yes" & all_clin$HER2.Prim=="Neg"), ]