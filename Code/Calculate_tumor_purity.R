rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")

#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)

library(estimate)  #estimate_1.0.13
#help(package="estimate")

### load log2 cpm file ###
load("Data_v2/log2tpm_unique.RData")
load("Data_v2/gene_annot_biomart_unique.RData")
log2tpm <- log2tpm_unique 
rownames(log2tpm) <- gene_annot_biomart_unique[,"external_gene_name_v2"]

# check overlap genes with ESTIMATE database
data(common_genes)
dim(common_genes) #10412     5
symbol_in <- rownames(log2tpm) %in% common_genes[,"GeneSymbol"]
sum(symbol_in)  #10033

##### Estimate purity
#write.table(log2tpm, file="Data_v2/log2tpm_for_ESTIMATE.txt", quote=FALSE,
#	sep="\t")
#filterCommonGenes(input.f="Data_v2/log2tpm_for_ESTIMATE.txt", 
#	output.f="Data_v2/log2tpm_from_ESTIMATE.gct", id="GeneSymbol")

# no purity estimate for illumina paltform
estimateScore("Data_v2/log2tpm_from_ESTIMATE.gct", 
	"Results_v2/log2tpm_from_ESTIMATE_score.gct", platform="illumina") 

estimateScore("Data_v2/log2tpm_from_ESTIMATE.gct", 
	"Results_v2/log2tpm_from_ESTIMATE_score_affy.gct", platform="affymetrix")

#### read in tumor purity
purity_illu <- read.delim("Results_v2/log2tpm_from_ESTIMATE_score.gct")
purity <- read.delim("Results_v2/log2tpm_from_ESTIMATE_score_affy.gct")
purity2 <- purity[3:nrow(purity), -c(1,2)]
rownames(purity2) <- purity[3:nrow(purity),1]

id <- c(unlist(purity[2,-c(1,2)]))
names(id) <- NULL
changename <- function(str){
	str <- gsub("\\.", "-", str)
	fc <- substr(str, 1, 1)
	if(fc=="X"){
		return(substr(str, 2, nchar(str)))
	}else{
		return(str)
	}
}
id2 <- unlist(sapply(1:length(id), function(x) changename(id[x])))

all(id2==colnames(log2tpm))

colnames(purity2) <- id2
estimate_purity_score <- data.matrix(purity2)
save(estimate_purity_score, file="Results_v2/log2tpm_estimate_purity_score.RData")

## compare with old version
load("Data/log2tpm_estimate_purity_score.RData")
old <- estimate_purity_score
load("Results_v2/log2tpm_estimate_purity_score.RData")
new <- estimate_purity_score
dim(old)
dim(new)
cor(old[4,], new[4,])  # 0.998

get_purity <- function(est_score){
	purity <- cos(0.6049872018 + 0.0001467884 * est_score)
	return(purity)
}
get_purity(new[3,1:5])



