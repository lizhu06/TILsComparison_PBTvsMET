rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")
library(DESeq2) 
packageVersion("DESeq2") #1.18.1

#### load data
load("Data_v2/count.RData")
dim(count) #60619   105

load("Data_v2/sample_annot.RData")
load("Data_v2/sample_pair_info.RData")

sample_pair_info <- sample_pair_info[
  which(sample_pair_info[,"mets_id"] != "7M_RCS"), ]
sample_annot <- sample_annot[which(
  !(sample_annot[,"ID"] %in% c("7M_RCS","7P_RCS"))), ]
count <- count[,which(!(colnames(count) %in% c("7M_RCS","7P_RCS")))]

#### gene annotation
load("Data_v2/gene_annot_biomart_unique.RData")
length(unique(rownames(count))) # 60619
length(unique(rownames(gene_annot_biomart_unique))) # 55557
dim(gene_annot_biomart_unique) #55557     8
all(rownames(gene_annot_biomart_unique) %in% rownames(count)) # TRUE

count_unique <- count[rownames(gene_annot_biomart_unique), ]
rownames(count_unique) <- gene_annot_biomart_unique[,"external_gene_name_v2"]
dim(count_unique) #55557   105

### label sample, same order as columns of count_unique
pri_met_label <- sample_annot[match(colnames(count_unique), 
	sample_annot[,"ID"]), "group"]
site_label <- sample_annot[match(colnames(count_unique),
	sample_annot[,"ID"]), "site"]
site_label[which(site_label=="Breast")] <- sample_pair_info[match(
	colnames(count_unique)[which(site_label=="Breast")], 
	sample_pair_info[,"primary_id"]), "mets_site"]
table(pri_met_label) # 53 MET, 50 PRI
table(site_label) # Bone 22, Brain 42, gi 12, ovary 27
names(pri_met_label) <- names(site_label) <- colnames(count_unique)

# test PRI vs MET (site specific)
library(BiocParallel)
register(MulticoreParam(5))
count_unique_round <- round(count_unique)
label <- relevel(as.factor(pri_met_label),ref="Primary")
colData <- data.frame(label=label)
uni_site <- c("Brain", "Bone", "GI", "Ovary")
for(i in 1:length(uni_site)){
	keep_index <- which(site_label==uni_site[i])
	count_site <- count_unique_round[, keep_index]
	colData_site <- colData[keep_index, , drop=FALSE]
	dds <- DESeqDataSetFromMatrix(countData=count_site,
	                                 colData=colData_site,
	                                 design= ~label)
	#t0<-proc.time()
	dds <- DESeq(dds, parallel = T)
	#t1<-proc.time()-t0

	res <- results(dds, contrast=c("label","MET","Primary"))
	save(res,file=paste0("Results_v2/DE_MET_vs_PRI_", uni_site[i], ".RData"))
}



