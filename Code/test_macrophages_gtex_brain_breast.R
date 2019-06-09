rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")

library(tidyr)
library(reshape2)
library(ggplot2)
library(multcomp)
library(ggpubr)

#load("Data_v2/gtex_our_pairs_log2tpm_combined.RData")
load("Data_v2/gtex_breast_brain_ov_si_tissue_label.RData")

load("Results_v3/gtex_raw_immune/gsva_davoli_gtex.RData")
load("Results_v3/gtex_raw_immune/gsva_tamborero_gtex.RData")
load("Results_v3/gtex_raw_immune/cibersort_absolute.RData")
load("Results_v3/gtex_raw_immune/cibersort_relative.RData")

all_macrophages <- data.frame(tissue=gtex_tissue, 
	Cibersort_M2=ciber_relative["Macrophages.M2",],
	GSVA_Davoli_M2=gsva_davoli_gtex["Macrophages.M2", ])
all(colnames(ciber_relative) == colnames(gsva_davoli_gtex))

data_long <- melt(all_macrophages, id.vars="tissue")

# Box plot of all four tissues
p <- ggboxplot(data_long, x = "tissue", y = "value",
          color = "tissue", palette = "jco", scales = "free",
          add = "jitter",
          facet.by = "variable", short.panel.labs = TRUE)
# Use only p.format as label. Remove method name.
#p <- p + stat_compare_means(label = "p.format", 
#  method = "wilcox.test")

#p[[s]] <- ggplot(data_long, aes(x=variable, y=value, fill=group)) + 
#  geom_boxplot() + 
#  facet_wrap(~variable, scale="free")
pdf(file=paste0("Results_v3/gtex_raw_immune/Boxplot_M2_gtex_brain_bone_gi_ov_breast", 
  ".pdf"), width=10, height=4)
print(p)
dev.off()  

# only cibersort 
all_macrophages <- data.frame(tissue=gtex_tissue, 
  Cibersort_M2=ciber_relative["Macrophages.M2",]*100,
  GSVA_Davoli_M2=gsva_davoli_gtex["Macrophages.M2", ])
all(colnames(ciber_relative) == colnames(gsva_davoli_gtex))

data_long <- melt(all_macrophages, id.vars="tissue")

# Box plot of all four tissues
p <- ggboxplot(data_long, x = "tissue", y = "value",
          color = "tissue", palette = "jco", scales = "free",
          add = "jitter",
          facet.by = "variable", short.panel.labs = TRUE)
# Use only p.format as label. Remove method name.
#p <- p + stat_compare_means(label = "p.format", 
#  method = "wilcox.test")

#p[[s]] <- ggplot(data_long, aes(x=variable, y=value, fill=group)) + 
#  geom_boxplot() + 
#  facet_wrap(~variable, scale="free")
pdf(file=paste0("Results_v3/gtex_raw_immune/Boxplot_M2_gtex_brain_bone_gi_ov_breast_percentage*100", 
  ".pdf"), width=10, height=4)
print(p)
dev.off()  


# Box plot of breast & brain
data_long_sub <- data_long[which(data_long$tissue %in% c("Breast", "Brain")), ]
p <- ggboxplot(data_long_sub, x = "tissue", y = "value",
          color = "tissue", palette = "jco", scales = "free",
          add = "jitter",
          facet.by = "variable", short.panel.labs = TRUE)
# Use only p.format as label. Remove method name.
p <- p + stat_compare_means( 
  method = "wilcox.test")

#p[[s]] <- ggplot(data_long, aes(x=variable, y=value, fill=group)) + 
#  geom_boxplot() + 
#  facet_wrap(~variable, scale="free")
pdf(file=paste0("Results_v3/gtex_raw_immune/Boxplot_M2_gtex_brain_breast", 
  ".pdf"), width=8, height=4)
print(p)
dev.off()  