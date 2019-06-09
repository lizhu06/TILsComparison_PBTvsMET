rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")
library(tidyr)
library(reshape2)
library(ggplot2)
library(multcomp)

## load  scores and clinical
load("RawData_v2/TILpaper/brain_mets_matched_tilcount.RData")
load("RawData_v2/TILpaper/primary_matched_tilcount.RData")
load("RawData_v2/TILpaper/clinical.RData")
all(primary_matched_tilcount[,"ID"]==brain_mets_matched_tilcount[,"ID"])
all(primary_matched_tilcount[,"ID"]==BRM_TIL_clin[,"ID"])
BRM_til <- brain_mets_matched_tilcount[,"Final.TIL.Analysis"]
PBT_til <- primary_matched_tilcount[,"Final.TIL.Analysis"]
delta_til <- BRM_til - PBT_til

# set color
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n <- 6
ggcols = gg_color_hue(n)

############################
## plot primary and mets
############################
data_plot <- data.frame(pair_id=primary_matched_tilcount[,"ID"], 
  PBT=primary_matched_tilcount[,"Final.TIL.Analysis"],
  BRM=brain_mets_matched_tilcount[,"Final.TIL.Analysis"])
data_long <- melt(data_plot, id.vars=c("pair_id")) 
colnames(data_long)[2] <- "Group"
p <- ggplot(data_long, aes(x=Group, y=value, fill=Group)) +
  geom_point(size=2, alpha=0.6, aes(color=Group))
p <- p + geom_line(aes(group=pair_id), alpha = 0.6, colour = "grey", data = data_long) 
p <- p + labs(x="", y="sTILs %") + 
  scale_colour_manual(values = ggcols[c(4,2)])

pdf(file="Results_v2/boxplot_sTILs_primary_mets.pdf", width=3, height=2.5)
plot(p)
dev.off()

# plot delta
df_immune_change <- data.frame(immune_change=brain_mets_matched_tilcount[,"Final.TIL.Analysis"] - 
  primary_matched_tilcount[,"Final.TIL.Analysis"])

p_delta <- ggplot(df_immune_change, aes(x="", y=immune_change)) +
  geom_boxplot() + geom_jitter(position=position_jitter(0.2))
#p <- p + geom_line(aes(group = pairs), alpha = 0.6, colour = "grey", data = data_long) 
p_delta <- p_delta + labs(x="", y="sTILs % change") 

pdf(file="Results_v2/boxplot_sTILs_change.pdf", width=2, height=2.5)
plot(p_delta)
dev.off()

wilcox.test(brain_mets_matched_tilcount[,"Final.TIL.Analysis"], 
  primary_matched_tilcount[,"Final.TIL.Analysis"], paired=TRUE)

# P=0.0002

######################
# HR and HER2 subtypes
#######################
table(BRM_TIL_clin[,"BCsubtype"])
HR_HER2 <- BRM_TIL_clin[,"BCsubtype"]
HR_HER2 <- gsub("HR-/HER2-", "TNBC", HR_HER2)
HR_HER2 <- factor(HR_HER2, levels=c("HR+/HER2+", "HR+/HER2-", "HR-/HER2+", 
  "TNBC"))
table(HR_HER2)
#HRpos_HER2pos  HRpos_HER2neg  HRneg_HER2pos      tripleNeg HRpos_HER2equi
#            8             25              3             10              4
uni_HR_HER2 <- levels(HR_HER2)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n <- 11
ggcols2 = gg_color_hue(n)

color_index <- c(c(7,2,9,3))
HR_HER2_vec <-c("HR+/HER2+", "HR+/HER2-", "HR-/HER2+", 
  "TNBC")
HR_HER2_name_vec <- HR_HER2_vec
HR_HER2_name_vec_fileName <- c("HR+_HER2+", "HR+_HER2-", "HR-_HER2+", 
  "TNBC")
#i <- 1
pval <- rep(NA, 4)
range <- range(PBT_til, BRM_til, delta_til)
for(i in 1:4){
  keep_index <- which(HR_HER2 == HR_HER2_vec[i])
  data_wide <- data.frame(pair_id=primary_matched_tilcount[keep_index,"ID"], 
    PBT=PBT_til[keep_index],
    BRM=BRM_til[keep_index], Change=delta_til[keep_index])
  
  data_long <- melt(data_wide, id.vars=c("pair_id")) 
  colnames(data_long)[2] <- "Group"
  data_long$Group <- factor(data_long$Group, levels=paste0(c("PBT", "BRM", "Change")))

  p <- ggplot(data_long, aes(x=Group, y=value, fill=Group)) +
    geom_point(size=2, alpha=0.6, aes(color=Group))
    
  #p <- p + geom_line(aes(group = pairs), alpha = 0.6, colour = "grey", data = data_long)

  p <- p + labs(x="", y="sTILs %") + 
    scale_colour_manual(values = c(ggcols[c(4)],ggcols2[i], "black"))
  for(t in 1:nrow(data_wide)){
    p <- p + geom_segment(x=1, y=data_wide[t, 2], 
        xend=2, yend=data_wide[t, 3], colour="grey") 
  }
  p <- p + ylim(range[1]-5, range[2]+5)
  p <- p+theme(axis.text.x = element_text(angle = 45, hjust = 1))
  pdf(file=paste0("Results_v2/boxplot_sTILs_BRM_TIL_", HR_HER2_name_vec_fileName[i], ".pdf"), 
    width=3, height=2.5)
  plot(p)
  dev.off()
  fit <- wilcox.test(data_wide[,2], 
    data_wide[,3], paired=TRUE)
  pval[i] <- fit$p.value
}

