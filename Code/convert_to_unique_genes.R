rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")

## load data
load("Data_v2/tpm.RData")
log2tpm <- log2(tpm + 1)
load("Data_v2/gene_annot_biomart.RData")
load("Data_v2/sample_annot.RData")

length(unique(gene_annot_biomart[,"external_gene_name_v2"])) #55644
dup_genes <- unique(gene_annot_biomart[duplicated(gene_annot_biomart[,"external_gene_name_v2"]), "external_gene_name_v2"])
dup_genes <- dup_genes[!is.na(dup_genes)]
length(dup_genes)  # 142
non_dup_gene_index <- which((!(gene_annot_biomart[,"external_gene_name_v2"] %in% dup_genes)) 
	& (!is.na(gene_annot_biomart[,"external_gene_name_v2"])) )
length(non_dup_gene_index) # 55415

select_dup_genes <- function(gene_name, counter, data, method){
	candi_index <- which(gene_annot_biomart[,"external_gene_name_v2"]==gene_name)
	sub_data <- data[candi_index, ]
	if(method=="IQR"){
		mea_vec <- apply(sub_data, 1, IQR)
		}else{
			mea_vec <- apply(sub_data, 1, sd)
		}
	print(counter)
	return(candi_index[which.max(mea_vec)])
}

dup_gene_keep_index <- sapply(1:length(dup_genes), function(x) 
	select_dup_genes(dup_genes[x], x, log2tpm, method="SD"))

gene_keep_index <- c(non_dup_gene_index, dup_gene_keep_index)
gene_annot_biomart_unique <- gene_annot_biomart[gene_keep_index, ]
log2tpm_unique <- log2tpm[gene_keep_index, ]
dim(log2tpm_unique) #55557   105
save(log2tpm_unique, file="Data_v2/log2tpm_unique.RData")
save(gene_annot_biomart_unique, file="Data_v2/gene_annot_biomart_unique.RData")


