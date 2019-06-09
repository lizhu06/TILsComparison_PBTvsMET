rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")
library(ggplot2)
# set color
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n <- 6
ggcols = gg_color_hue(n)

load("Data_v2/sample_annot.RData")
load("Results_v2/log2tpm_estimate_purity_score.RData")
til <- read.csv("RawData_v2/PittTILData.csv")
til_sample <- til[,"Sample_."]
til_sample <- sub("_Pitt", "", til_sample)
til_id <- sapply(1:length(til_sample), function(x) 
	substr(til_sample[x], 1, nchar(til_sample[x])-1))
til_type <- sapply(1:length(til_sample), function(x) 
	substr(til_sample[x], nchar(til_sample[x]), nchar(til_sample[x])))
head_str <- rep("BP", length(til_id))
head_str[til_type=="M"] <- "BM"
til_comp_id <- paste(head_str, til_id, sep="")
all(til_comp_id %in% sample_annot[,"ID"]) 
til_comp_id[which(!(til_comp_id %in% sample_annot[,"ID"]))] # BP/BM 28

rownames(til) <- til_comp_id
common_id <- intersect(til_comp_id, sample_annot[,"ID"])
common_site <- sample_annot[match(common_id, 
  sample_annot[,"ID"]), "site"]

overlap_samples <- data.frame(ID=common_id, Site=common_site)
write.csv(overlap_samples, file="Results_v2/overlap_samples.csv")

data_plot <- data.frame(
  Estimate_immune=estimate_purity_score["ImmuneScore", common_id],
	TIL_count=til[common_id,"Final.TIL.Analysis...."], 
  sample_type=til[common_id, "Sample.Location"],
  norm_immune=estimate_purity_score["ImmuneScore", common_id]/
  (1-estimate_purity_score["TumorPurity", common_id]))
corr <- cor(data_plot[,-3], method="spearman")
cor.test(data_plot[,1], data_plot[,2], method="spearman")
data_plot$sample_type <- gsub("breast primary", "BPT", data_plot$sample_type)
data_plot$sample_type <- gsub("brain met", "BRM", data_plot$sample_type)

data_plot$sample_type <- factor(data_plot$sample_type, 
  levels=c("BPT", "BRM"))

p <- ggplot(data_plot, aes(y=Estimate_immune, x=TIL_count, group=sample_type)) + 
  geom_point(aes(color=sample_type))+
  ggtitle(paste0("corr=", 
  	round(corr[1,2], digits=2))) + 
  scale_colour_manual(values = ggcols[c(4,2)])+
  #annotate("text", x=20, y=0, label =paste0("spearman corr=", 
  #	round(corr[1,2], digits=2)),
  # size=4, color="red") +
  labs(x="TIL count", y="ESTIMATE immune score")

pdf("Results_v2/TIL_vs_ESTMATE_immune_score.pdf", width=3.2, height=2.5)
p
dev.off()


p <- ggplot(data_plot, aes(y=norm_immune, x=TIL_count)) + 
  geom_point(aes(color=sample_type))+
  ggtitle(paste0("corr=", 
    round(corr[2,3], digits=2))) + 
  scale_colour_manual(values = ggcols[c(4,2)])+
  #annotate("text", x=20, y=0, label =paste0("spearman corr=", 
  # round(corr[1,2], digits=2)),
  # size=4, color="red") +
  labs(x="TIL count", y="ESTIMATE immune score/(1-purity)")

pdf("Results_v2/TIL_vs_normalized_immune.pdf", width=3.2, height=2.5)
p
dev.off()
