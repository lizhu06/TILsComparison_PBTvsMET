rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")

## Gtex tpm
load("Data/Gtex_n100_combined_tissues_tissues_keep.RData")
load("Data/Gtex_n100_combined_tissues_log2tpm_keep.RData")
all(colnames(log2tpm_keep) == names(tissue_keep))

unique(tissue_keep)
print(paste(unique(tissue_keep), collapse=", "))
tissue_to_keep <- c("Breast", "Brain", "Ovary", "Small Intestine")
gtex <- log2tpm_keep[, tissue_keep %in% tissue_to_keep]
gtex_tissue <- tissue_keep[tissue_keep %in% tissue_to_keep]
dim(gtex) #54271   400
save(gtex, file="Data_v2/gtex_breast_brain_ov_si.RData")
save(gtex_tissue, file="Data_v2/gtex_breast_brain_ov_si_tissue_label.RData")

## load our pairs
load("Data_v2/log2tpm_unique.RData")
exc_id <- c("7M_RCS","7P_RCS")
log2tpm_unique <- log2tpm_unique[, !(colnames(log2tpm_unique) %in% exc_id)]
dim(log2tpm_unique) #55557   103

## normal
load("Data/normal_tpm.RData")
normal_log2tpm <- log2(normal_tpm+1)
normal_log2tpm_unique <- normal_log2tpm[rownames(log2tpm_unique), ]
dim(normal_log2tpm_unique) #55643    10

load("Data_v2/gene_annot_biomart_unique.RData")
rownames(log2tpm_unique) <- rownames(normal_log2tpm_unique) <- gene_annot_biomart_unique[,"external_gene_name_v2"]

load("Data_v2/sample_annot.RData")
pairs_tissue <- sample_annot[match(colnames(log2tpm_unique), sample_annot[,"ID"]), "site"]

## combined all
common_genes <- intersect(rownames(gtex), rownames(log2tpm_unique))
length(common_genes) #33275
setdiff(rownames(gtex), rownames(pairs))[1:10]
setdiff(rownames(pairs), rownames(gtex))[1:10]
"ENSG00000227234" %in% rownames(gene_annot_biomart_unique)

gtex_common <- gtex[common_genes, ]
pri_mets_common <- log2tpm_unique[common_genes, ]
normal_common <- normal_log2tpm_unique[common_genes, ]
normal_names <- gsub("NM", "Normal", colnames(normal_common))
normal_names <- gsub("OV", "ovary", normal_names)
normal_names <- gsub("BT", "breast", normal_names)
normal_names <- sapply(1:length(normal_names), function(x) 
	paste(strsplit(normal_names[x], split="_")[[1]][1:2], collapse="_"))

all_log2tpm <- cbind(gtex_common, pri_mets_common, normal_common)
all_tissues <- c(paste0("GTEx_", gtex_tissue), paste0("Tumor_", pairs_tissue), 
	normal_names)
all_tissues <- gsub(" ", "_", all_tissues)

## store data
save(all_log2tpm, file="Data_v2/gtex_our_pairs_log2tpm_combined.RData")
save(all_tissues, file="Data_v2/gtex_our_pairs_combined_labels.RData")
#save(all_log2tpm_rmBatch, file="Data_v2/gtex_our_pairs_log2tpm_combined_combat.RData")










