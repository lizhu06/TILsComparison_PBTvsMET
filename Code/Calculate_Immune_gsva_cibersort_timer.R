rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")
# R 3.4.0
library(edgeR) # 3.20.4
library(GSVA)  # 1.26.0

#### load data
load("RawData_v2/txi.salmon082.PanMets.GeneTPM.Rda")
load("Data_v2/log2tpm_unique.RData")
log2tpm <- log2tpm_unique
tpm <- tpm[rownames(log2tpm_unique), ] # has more samples!
load("Data_v2/sample_annot.RData")
load("Data_v2/sample_pair_info.RData")
#load("Data/estimate_purity_score.RData")
load("Data_v2/gene_annot_biomart_unique.RData")
rownames(log2tpm) <- rownames(tpm) <- gene_annot_biomart_unique[,"external_gene_name_v2"]
dim(log2tpm) #55557   105

################################
### compute davoli signatures
################################
load("/net/wong05/home/liz86/Steffi/Kevin_IDC_ILC_DE/RawData/Immune_cell_signatures_Davoli2016.RData")

## check if genes duplicated
davoli_gens <- unlist(immune_cell_signatures_davoli)
sapply(1:length(davoli_gens), function(x) sum(davoli_gens[x] %in% rownames(log2tpm)))

## check if any genes are not in the dataset
immune_cell_signatures_davoli[[4]] <- gsub("FAIM3", "FCMR", 
	immune_cell_signatures_davoli[[4]])
sapply(1:length(immune_cell_signatures_davoli), function(x) 
	immune_cell_signatures_davoli[[x]][!(immune_cell_signatures_davoli[[x]] 
		%in% rownames(log2tpm))])
davoli_genes_common <- immune_cell_signatures_davoli
save(davoli_genes_common, file="Data_v2/davoli_genes_common.RData")

## GSVA (use all genes, log2tpm)
fit <- gsva(log2tpm, davoli_genes_common, rnaseq=FALSE)
gsva_score_davoli_log2tpm <- fit
save(gsva_score_davoli_log2tpm, 
	file="Results_v2/gsva_score_davoli_log2tpm.RData")

################################
# Tamborero signature
###############################
load("/net/wong05/home/liz86/Steffi/Kevin_IDC_ILC_DE/RawData/Immune_cell_signatures_Tamborero2017.RData")
sapply(1:length(immune_cell_signatures_tamborero2017), function(x) 
	length(immune_cell_signatures_tamborero2017[[x]]))

## check if genes duplicated
tam_gens <- unlist(immune_cell_signatures_tamborero2017)
sapply(1:length(tam_gens), function(x) sum(tam_gens[x] %in% rownames(log2tpm)))

# To alias
for(i in 1:length(immune_cell_signatures_tamborero2017)){
	out_gene <- immune_cell_signatures_tamborero2017[[i]][!
	(immune_cell_signatures_tamborero2017[[i]] %in% 
		rownames(log2tpm))]
	if(length(out_gene)>0){
		symbol2 <- alias2SymbolTable(out_gene, species="Hs")
		symbol2[is.na(symbol2)] <- "NA2"
		symbol2_in <- symbol2 %in% rownames(log2tpm)
		out_gene_hit <- out_gene[symbol2_in]
		out_gene_hit_symbol <- symbol2[symbol2_in]
		immune_cell_signatures_tamborero2017[[i]][match(out_gene_hit,
			immune_cell_signatures_tamborero2017[[i]])] <- out_gene_hit_symbol
	}
}

sapply(1:length(immune_cell_signatures_tamborero2017), function(x) 
	immune_cell_signatures_tamborero2017[[x]][!(immune_cell_signatures_tamborero2017[[x]] %in% 
		rownames(log2tpm))])

## manually correct gene symbols
immune_cell_signatures_tamborero2017[[5]] <- gsub("6-Mar", "MARCH6",
	immune_cell_signatures_tamborero2017[[5]])

sapply(1:length(immune_cell_signatures_tamborero2017), function(x) 
	sum(!(immune_cell_signatures_tamborero2017[[x]] %in% 
		rownames(log2tpm))))
tamborero_genes_common <- immune_cell_signatures_tamborero2017
save(tamborero_genes_common, file="Data_v2/tamborero_genes_common.RData")

## GSVA (log2tpm)
fit <- gsva(log2tpm, tamborero_genes_common, rnaseq=FALSE)
gsva_score_tamborero_log2tpm <- fit
save(gsva_score_tamborero_log2tpm, file="Results_v2/gsva_score_tamborero_log2tpm.RData")

#################################
# Prepare for TIMER
#################################
tpm_timer <- cbind(rownames(tpm), tpm)
write.csv(tpm_timer, file="Data_v2/tpm_timer.csv", 
	quote=FALSE, row.names=FALSE) # file format doesn't work for cibersort

# load results
timer_raw <- read.csv("Results_v2/TIMER_res_pri_mets.csv")
timer <- timer_raw[, -1]
rownames(timer) <- timer_raw[,1]
timer <- t(timer)
save(timer,file="Results_v2/TIMER_res_pri_mets.RData")


#################################
# Prepare for CIBERSORT
#################################
tpm_for_cibersort <- cbind(rownames(tpm), tpm)
write.table(tpm_for_cibersort, file="Data_v2/tpm_for_cibersort.txt",
	quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

# load results
ciber_raw <- read.csv("Results_v2/CIBERSORT_relative.csv")
ciber_relative <- ciber_raw[,seq(2, 23)]
rownames(ciber_relative) <- ciber_raw[,1]
ciber_relative <- t(ciber_relative)
save(ciber_relative, file="Results_v2/cibersort_relative.RData")

# load results
ciber_raw <- read.csv("Results_v2/CIBERSORT_absolute.csv")
ciber_absolute <- ciber_raw[,seq(2, 23)]
rownames(ciber_absolute) <- ciber_raw[,1]
ciber_absolute <- t(ciber_absolute)
save(ciber_absolute, file="Results_v2/cibersort_absolute.RData")











