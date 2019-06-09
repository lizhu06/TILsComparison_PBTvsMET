## Re-annotate the genes using Biomart
## Old annotation is outdated
## Li Zhu, 4/18/2018

rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")

library(biomaRt) # 2.32.1 (wong08 2.34.2)

load("RawData_v2/PanMets.Paired.log2.tmmNorm.cpm.salmon.0.8.2.Rda")
cpm <- tmp.2[,1:(ncol(tmp.2)-4)]
annot <- tmp.2[,(ncol(tmp.2)-3):ncol(tmp.2)]
dim(annot) #60619     4
length(unique(annot[,"external_gene_name"])) #55644

tb <- table(annot[,"gene_biotype"])
tb2 <- tb[order(tb, decreasing=TRUE)]
# protein_coding 19,719
# processed_pseudogene 10,196
# lincRNA 7,403
# antisense 5,448
# unprocessed_pseudogene 2,577

tb <- table(annot[,"external_gene_source"])
tb2 <- tb[order(tb, decreasing=TRUE)]
# HGNC symbol 36103
# clone based (ensembl) 19195

### annot using biomart
#listMarts()
#ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") # sometimes not work
#db <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol',
#  'external_gene_name', "external_gene_source"), mart= ensembl)
#dim(db) #64616     4
#save(db, file="RawData/ensembl_db.RData") # saved on 4/18/2018

load("RawData_v2/ensembl_db.RData")

annot_biomart <- db[match(rownames(annot),db[,"ensembl_gene_id"]), ]
rownames(annot_biomart) <- NULL

### 201 gene symbols are different
all(annot[,"external_gene_name"] == annot_biomart[,"hgnc_symbol"]
	, na.rm=TRUE)
symbol_diff_index <- which(annot[,"external_gene_name"] != 
	annot_biomart[,"external_gene_name"] & (!is.na(annot[,"external_gene_name"]))
	& (annot_biomart[,"external_gene_name"] != ""))
length(symbol_diff_index)  #3759
annot[symbol_diff_index[1:5],]
annot_biomart[symbol_diff_index[1:5], ]

### combine two gene annotation
colnames(annot_biomart) <- paste0(colnames(annot_biomart), "_v2")
gene_annot_biomart <- cbind(annot, annot_biomart)

save(gene_annot_biomart, file="Data_v2/gene_annot_biomart.RData")






