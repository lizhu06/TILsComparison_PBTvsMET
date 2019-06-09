rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")
library(tidyr)
library(reshape2)
library(ggplot2)
library(multcomp)

## load  scores
load("Results_v2/immune_res.RData")

#all_res_norm <- lapply(1:length(all_res), function(x) sweep(all_res[[x]],
#  1, 1-estimate_purity_score["TumorPurity", id], "/"))
immune <- all_res$estimate[,"ImmuneScore", drop=FALSE]

## load sample data
load("Data_v2/sample_pair_info.RData")
load("Data_v2/sample_annot.RData")
sample_pair_info <- sample_pair_info[
  which(sample_pair_info[,"mets_id"] != "7M_RCS"), ]
sample_annot <- sample_annot[which(
  !(sample_annot[,"ID"] %in% c("7M_RCS","7P_RCS"))), ]

# create data for plot
sample_pair_info$primary_immune <- immune[sample_pair_info$primary_id, 1]
sample_pair_info$mets_immune <- immune[sample_pair_info$mets_id, 1]
sample_pair_info$immune_change <- sample_pair_info$mets_immune - 
  sample_pair_info$primary_immune


# check immune scores of primaries that matched to 2 mets
dup_pri_id <- sample_pair_info[duplicated(sample_pair_info$primary_id), "primary_id"]
sample_pair_info[which(sample_pair_info[,"primary_id"] %in% dup_pri_id), ]
# remove one of the mets
#sample_pair_info <- sample_pair_info[!duplicated(sample_pair_info$primary_id), ]

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
data_plot <- sample_pair_info[, c("pairs", "mets_site", "primary_immune",
 "mets_immune")]
colnames(data_plot)[c(3,4)] <- c("PBT", "MET")

data_long <- melt(data_plot, id.vars=c("pairs", "mets_site")) 
colnames(data_long)[3] <- "Group"
p <- ggplot(data_long, aes(x=Group, y=value, fill=Group)) +
  geom_point(size=2, alpha=0.6, aes(color=Group))
p <- p + geom_line(aes(group = pairs), alpha = 0.6, colour = "grey", data = data_long) 
p <- p + labs(x="", y="Immune ssGSEA score") + 
  scale_colour_manual(values = ggcols[c(4,1)])

pdf(file="Results_v3/test_plot_raw_total_immune_treatAsDiffPair/boxplot_raw_totalImmune_primary_mets.pdf", width=3, height=2.5)
plot(p)
dev.off()

# plot delta
df_immune_change <- sample_pair_info[, c("pairs", "mets_site", "immune_change")]

p_delta <- ggplot(df_immune_change, aes(x="", y=immune_change)) +
  geom_boxplot() + geom_jitter(position=position_jitter(0.2))
#p <- p + geom_line(aes(group = pairs), alpha = 0.6, colour = "grey", data = data_long) 
p_delta <- p_delta + labs(x="", y="Immune ssGSEA score change") 

pdf(file="Results_v3/test_plot_raw_total_immune_treatAsDiffPair/boxplot_raw_totalImmune__change.pdf", width=2, height=2.5)
plot(p_delta)
dev.off()

wilcox.test(sample_pair_info$primary_immune, sample_pair_info$mets_immune, paired=TRUE)

# P=0.0001

################
# Site specific
#################
## plot primary and mets
color_index <- c(2,3,5,6)
site_vec <- c("Brain", "Ovary", "Bone", "GI")
site_name_vec <- c("BRM", "OVM", "BOM", "GIM")
#i <- 1
pval <- rep(NA, 4)
rangee <- range(sample_pair_info[, c("primary_immune", "mets_immune", "immune_change")])
for(i in 1:4){
  data_plot <- sample_pair_info[which(sample_pair_info$mets_site==site_vec[i]), 
    c("pairs", "mets_site", "primary_immune", "mets_immune", "immune_change")]

  colnames(data_plot)[c(3,4,5)] <- c("PBT", site_name_vec[i], "Change")

  data_long <- melt(data_plot, id.vars=c("pairs", "mets_site")) 
  colnames(data_long)[3] <- "Group"
  data_long$Group <- factor(data_long$Group, levels=c("PBT", site_name_vec[i], "Change"))

  p <- ggplot(data_long, aes(x=Group, y=value, fill=Group)) +
    geom_point(size=2, alpha=0.6, aes(color=Group))
    
  #p <- p + geom_line(aes(group = pairs), alpha = 0.6, colour = "grey", data = data_long)

  p <- p + labs(x="", y="Immune ssGSEA score") + 
    scale_colour_manual(values = c(ggcols[c(4,color_index[i])], "black"))
  for(t in 1:nrow(data_plot)){
    p <- p + geom_segment(x=1, y=data_plot[t, 3], 
        xend=2, yend=data_plot[t, 4], colour="grey") 
  }
  p <- p + ylim(rangee[1]-100, rangee[2]+100)
  pdf(file=paste0("Results_v3/test_plot_raw_total_immune_treatAsDiffPair/boxplot_raw_totalImmune_primary_mets_", site_vec[i], ".pdf"), 
    width=3, height=2.5)
  plot(p)
  dev.off()
  fit <- wilcox.test(data_plot[,3], 
    data_plot[,4], paired=TRUE)
  pval[i] <- fit$p.value
}

######################
# HR and HER2 subtypes
#######################
sample_pair_info$HR_HER2 <- rep(NA, nrow(sample_pair_info))

sample_pair_info$HR_HER2[which(sample_pair_info$HR_prim=="Pos" & 
  sample_pair_info$HER2_prim=="Pos")] <- "HR+/HER2+"

sample_pair_info$HR_HER2[which(sample_pair_info$HR_prim=="Pos" & 
  sample_pair_info$HER2_prim=="Neg")] <- "HR+/HER2-"

sample_pair_info$HR_HER2[which(sample_pair_info$HR_prim=="Neg" & 
  sample_pair_info$HER2_prim=="Pos")] <- "HR-/HER2+"

sample_pair_info$HR_HER2[which(sample_pair_info$HR_prim=="Neg" & 
  sample_pair_info$HER2_prim=="Neg")] <- "TNBC"

sample_pair_info$HR_HER2[which(sample_pair_info$HR_prim=="Pos" & 
  sample_pair_info$HER2_prim=="Equi")] <- "HR+/HER2="

HR_HER2 <- sample_pair_info$HR_HER2
HR_HER2 <- factor(HR_HER2, levels=c("HR+/HER2+", "HR+/HER2-", "HR-/HER2+", 
  "TNBC", "HR+/HER2="))
table(HR_HER2)
#HRpos_HER2pos  HRpos_HER2neg  HRneg_HER2pos      tripleNeg HRpos_HER2equi
#            8             25              3             10              4
uni_HR_HER2 <- levels(HR_HER2)

table(sample_pair_info[which(sample_pair_info$mets_site=="Brain"), "HR_HER2"])
#HRneg_HER2pos HRpos_HER2neg HRpos_HER2pos     tripleNeg
#            3             5             5             8
uni_HR_HER2_brain <- c("HR+/HER2+", "HR+/HER2-", "HR-/HER2+", 
  "TNBC")


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
rangee <- range(sample_pair_info[, c("primary_immune", "mets_immune", "immune_change")])
for(i in 1:4){
  data_plot <- sample_pair_info[which(sample_pair_info$HR_HER2==HR_HER2_vec[i]), 
    c("pairs", "mets_site", "primary_immune", "mets_immune", "immune_change")]

  colnames(data_plot)[c(3,4,5)] <- paste0(c(" PBT", " MET", " Change"))

  data_long <- melt(data_plot, id.vars=c("pairs", "mets_site")) 
  colnames(data_long)[3] <- "Group"
  data_long$Group <- factor(data_long$Group, levels=paste0(c(" PBT", " MET", " Change")))

  p <- ggplot(data_long, aes(x=Group, y=value, fill=Group)) +
    geom_point(size=2, alpha=0.6, aes(color=Group))
    
  #p <- p + geom_line(aes(group = pairs), alpha = 0.6, colour = "grey", data = data_long)

  p <- p + labs(x="", y="Immune ssGSEA score") + 
    scale_colour_manual(values = c(ggcols[c(4)],ggcols2[i], "black"))
  for(t in 1:nrow(data_plot)){
    p <- p + geom_segment(x=1, y=data_plot[t, 3], 
        xend=2, yend=data_plot[t, 4], colour="grey") 
  }
  p <- p + ylim(rangee[1]-100, rangee[2]+100)
  p <- p+theme(axis.text.x = element_text(angle = 45, hjust = 1))
  pdf(file=paste0("Results_v3/test_plot_raw_total_immune_treatAsDiffPair/boxplot_raw_totalImmune_primary_mets_", HR_HER2_name_vec_fileName[i], ".pdf"), 
    width=3, height=2.5)
  plot(p)
  dev.off()
  fit <- wilcox.test(data_plot[,3], 
    data_plot[,4], paired=TRUE)
  pval[i] <- fit$p.value
}

