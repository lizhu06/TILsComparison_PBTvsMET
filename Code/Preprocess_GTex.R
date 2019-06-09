rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")

## Gtex
## tpm
gtex <- read.table("RawData/GTEx_20171117/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct", 
	skip=2, stringsAsFactors=F, header=T)
gtex_tpm <- gtex[,-c(1,2)]
dim(gtex_tpm) #56202 11688
gtex_annot <- gtex[,c(1,2)]

## count
gtex <- read.table("RawData/GTEx_20171117/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct", 
	skip=2, stringsAsFactors=F, header=T)
gtex_count <- gtex[,-c(1,2)]
dim(gtex_count) #56202 11688
all(gtex_annot == gtex[,c(1,2)]) # TRUE
all(colnames(gtex_count) == colnames(gtex_tpm)) # TRUE

## samples
samples <- read.delim('RawData/GTEx_20171117/GTEx_v7_Annotations_SampleAttributesDS.txt', sep='\t')
dim(samples) #15598    63
samples[,"SAMPID"] <- gsub('-', '\\.', samples[,"SAMPID"])
all(colnames(gtex_tpm) %in% samples[,"SAMPID"]) # TRUE
length(unique(samples[,"SAMPID"])) # 15598
samples_annot <- samples[match(colnames(gtex_tpm), samples[,"SAMPID"]), ]
dim(samples_annot) #11688    63
tissues <- samples_annot[,"SMTSD"]
table(tissues)
names(tissues) <- samples_annot[,"SAMPID"]

unique_genes <- unique(gtex_annot[,2])
length(unique_genes) #54271
dup_genes <- unique(gtex_annot[,2][duplicated(gtex_annot[,2])])
length(dup_genes) # 197

no_dup_genes_index <- which(!(gtex_annot[,2] %in% dup_genes))
length(no_dup_genes_index) # 54074
54074 + 197  #54271

return_largest_iqr_index <- function(gene_id, counter, data, method){
	candi_index <- which(gtex_annot[,2]==gene_id)
	sub_data <- data[candi_index, ]
	if(method=="IQR"){
		mea_vec <- apply(sub_data, 1, IQR)
		}else{
			mea_vec <- apply(sub_data, 1, sd)
		}
	print(counter)
	return(candi_index[which.max(mea_vec)])
}

return_sum <- function(gene_id, data){
	candi_index <- which(gtex_annot[,2]==gene_id)
	sub_data <- data[candi_index, ]
	return(apply(sub_data,2, sum))
}

if(FALSE){
	t0 <- proc.time()
	dup_gene_keep_index <- sapply(1:10, function(x) 
		return_largest_iqr_index(dup_genes[x], x, gtex_tpm, method="SD"))
	proc.time()-t0
	t0 <- proc.time()
	dup_gene_keep_index <- sapply(1:10, function(x) 
		return_largest_iqr_index(dup_genes[x], x, gtex_tpm, method="IQR"))
	proc.time()-t0
}


dup_gene_keep_index <- sapply(1:length(dup_genes), function(x) 
	return_largest_iqr_index(dup_genes[x], x, gtex_tpm, method="SD"))

dup_gene_sum <- sapply(1:length(dup_genes), function(x) 
	return_sum(dup_genes[x], gtex_tpm))

gene_keep_index <- c(no_dup_genes_index, dup_gene_keep_index)

gtex_gene_annot_uni_genes <- gtex_annot[gene_keep_index, ]
gtex_tpm_uni_genes <- gtex_tpm[gene_keep_index, ]
gtex_count_uni_genes <- gtex_count[gene_keep_index, ]
rownames(gtex_tpm_uni_genes) <- rownames(gtex_count_uni_genes) <- gtex_gene_annot_uni_genes[,2]

gtex_tpm_uni_genes_sum <- rbind(gtex_tpm[no_dup_genes_index, ], t(dup_gene_sum))
rownames(gtex_tpm_uni_genes_sum) <- c(gtex_annot[no_dup_genes_index,"Description"], dup_genes)

save(gtex_tpm_uni_genes_sum, file="Data/gtex_tpm_uni_genes_sum.RData")
save(gtex_gene_annot_uni_genes, file="Data/gtex_gene_annot_uni_genes.RData")
save(gtex_tpm_uni_genes, file="Data/gtex_tpm_uni_genes.RData")
save(gtex_count_uni_genes, file="Data/gtex_count_uni_genes.RData")
save(samples_annot, file="Data/gtex_sample_annot.RData")
save(tissues, file="Data/gtex_tissues.RData")

## not run below

if(FALSE){
	load("Data/gtex_tpm_uni_genes_sum.RData")
	load("Data/gtex_tpm_uni_genes.RData")
	all(rownames(gtex_tpm_uni_genes_sum)==rownames(gtex_tpm_uni_genes)
	gtex_tpm_uni_genes_sum[1:5,1:5]
	gtex_tpm_uni_genes[1:5,1:5]

}

gtex_tpm_uni_genes["TRAT1",100:120]
gtex_tpm[which(gtex_annot[,2]=="TRAT1"), "GTEX.13SLW.0226.SM.5S2NA"]

## HPA (bone marrow, only one sample)
hpa <- read.delim("RawData/hpa_tissue_20171117/rna_tissue.tsv")
bone <- hpa[which(hpa[,"Sample"]=="bone marrow"), ]
save(bone, file="Data/hpa_bone_tpm.RData")






