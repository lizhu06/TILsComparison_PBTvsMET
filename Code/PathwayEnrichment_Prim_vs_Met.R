rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")

# load Pathways
source("/home05/liz86/Steffi/Functions/Pathway_Fisher_Function.R")
# load database
library('GSA')
database = GSA.read.gmt("RawData_v2/genesets2531.gmt")
database2 = database$genesets
names(database2) = database$geneset.names

## primary vs mets
load("Results_v2/DE_MET_vs_PRI.RData")
sum(res[,"padj"] < 0.05 & res[,"log2FoldChange"] > 0, na.rm=TRUE)
sum(res[,"padj"] < 0.05 & res[,"log2FoldChange"] < 0, na.rm=TRUE)
sig <- rownames(res)[which(res[,"padj"] < 0.05)]
whole<-rownames(res)
length(sig)   #2695
length(whole)  #55557
pathwayRes <- pathway_fisher_detail_gene(significant=sig,
  whole=whole,fdr=0.05,database=database2)
print(sum(unlist(pathwayRes[,3])<0.05))  #695
pathwayRes_sort <- pathwayRes[order(as.numeric(pathwayRes[,3])),]
write.csv(pathwayRes_sort,
  file="Results_v2/PathwayRes_MET_vs_PRI.csv")

## site specific
uni_site <- c("Brain", "Bone", "GI", "Ovary")
for(i in 1:length(uni_site)){
	load(paste0("Results_v2/DE_MET_vs_PRI_", uni_site[i], ".RData"))
	sig <- rownames(res)[which(res[,"padj"] < 0.05)]
	whole<-rownames(res)
	length(sig)   
	length(whole)  
	pathwayRes <- pathway_fisher_detail_gene(significant=sig,
	  whole=whole,fdr=0.05,database=database2)
	print(sum(unlist(pathwayRes[,3])<0.05))  #695
	pathwayRes_sort <- pathwayRes[order(as.numeric(pathwayRes[,3])),]
	write.csv(pathwayRes_sort,
	  file=paste0("Results_v2/PathwayRes_MET_vs_PRI_", uni_site[i], ".csv"))
}

## brain ER+
load("Results_v2/DE_MET_vs_PRI_Brain_ERpos.RData")
sig <- rownames(res)[which(res[,"padj"] < 0.05)]
whole<-rownames(res)
length(sig)   #2695
length(whole)  #55557
pathwayRes <- pathway_fisher_detail_gene(significant=sig,
  whole=whole,fdr=0.05,database=database2)
print(sum(unlist(pathwayRes[,3])<0.05))  #695
pathwayRes_sort <- pathwayRes[order(as.numeric(pathwayRes[,3])),]
write.csv(pathwayRes_sort,
  file="Results_v2/PathwayRes_MET_vs_PRI_Brain_ERpos.csv")
