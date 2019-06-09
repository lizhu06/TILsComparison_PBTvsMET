rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")
library(ggplot2)
library(reshape2)
library(gridExtra)
load("Results_v2/OV_mets_staining_score.RData")
ovstain <- score
load("Results_v2/brain_mets_staining_score.RData")
brainstain <- score
load("Results_v2/OV_mets_stainedSample_bioinfo_score.RData")
ovciber <- ciber_rel_sub
load("Results_v2/Brain_mets_stainedSample_bioinfo_score.RData")
brainciber <- ciber_rel_sub

# calculate delta for staining results
ov_pair_id <- unique(ovstain[,1])
ov_stain_delta <- matrix(NA, 3, 4)
for(i in 1:length(ov_pair_id)){
  for(j in 1:4){
    mets_score <- as.numeric(ovstain[which(
      ovstain[,"pair_id"]==ov_pair_id[i] & 
      ovstain[,"group_id"]=="MET"), j+3])
    pri_score <- as.numeric(ovstain[which(
      ovstain[,"pair_id"]==ov_pair_id[i] & 
      ovstain[,"group_id"]=="PRI"), j+3])
    ov_stain_delta[i,j] <- mean(mets_score) - mean(pri_score)
  }
}
ov_pair_id
rownames(ov_stain_delta) <- c("656M-3150T", "OM1-OP1", "OM3-OP3")
colnames(ovstain)
colnames(ov_stain_delta) <- c("Tcell", "Treg", 
  "Macrophages", "Bcell")

brain_pair_id <- unique(brainstain[,1])
brain_stain_delta <- matrix(NA, 3, 4)
for(i in 1:length(brain_pair_id)){
  for(j in 1:4){
    mets_score <- as.numeric(brainstain[which(
      brainstain[,"pair_id"]==brain_pair_id[i] & 
      brainstain[,"group_id"]=="MET"), j+3])
    pri_score <- as.numeric(brainstain[which(
      brainstain[,"pair_id"]==brain_pair_id[i] & 
      brainstain[,"group_id"]=="PRI"), j+3])
    brain_stain_delta[i,j] <- mean(mets_score) - mean(pri_score)
  }
}
brain_pair_id
rownames(brain_stain_delta) <- c("BM72-BP72", "BM71-BP71", "BM52-BP52")
colnames(brainstain)
colnames(brain_stain_delta) <- c("Bcell", "Tcell", 
  "Macrophages", "Treg")
brain_stain_delta <- brain_stain_delta[, match(colnames(ov_stain_delta), 
  colnames(brain_stain_delta))]

# calculate delta in bioinfo results
ov_ciber_delta <- ovciber[, seq(4,6)] - ovciber[, seq(1,3)]
colnames(ov_ciber_delta) <- c("656M-3150T", "OM1-OP1", "OM3-OP3")
ov_ciber_delta <- t(ov_ciber_delta)
ov_ciber_delta <- ov_ciber_delta[match(rownames(ov_stain_delta), 
  rownames(ov_ciber_delta)), match(colnames(ov_stain_delta), 
  colnames(ov_ciber_delta))]

brain_ciber_delta <- brainciber[, seq(4,6)] - brainciber[, seq(1,3)]
colnames(brain_ciber_delta) <- c("BM52-BP52", "BM71-BP71", "BM72-BP72")
brain_ciber_delta <- t(brain_ciber_delta)
brain_ciber_delta <- brain_ciber_delta[match(rownames(brain_stain_delta), 
  rownames(brain_ciber_delta)), match(colnames(brain_stain_delta),
  colnames(brain_ciber_delta))]

## calculate correlation
ov_cor_vec <- rep(NA, 4)
ov_cor_all <- cor(c(ov_ciber_delta), c(ov_stain_delta), method="spearman")
for(j in 1:ncol(ov_ciber_delta)){
  ov_cor_vec[j] <- cor(ov_ciber_delta[,j], ov_stain_delta[,j], 
    method="spearman")
}

brain_cor_vec <- rep(NA, 4)
brain_cor_all <- cor(c(brain_ciber_delta), c(brain_stain_delta), method="spearman")
for(j in 1:ncol(brain_ciber_delta)){
  brain_cor_vec[j] <- cor(brain_ciber_delta[,j], brain_stain_delta[,j], 
    method="spearman")
}
names(ov_cor_vec) <- names(brain_cor_vec) <- colnames(ov_ciber_delta)
ov_cor_vec 
#      Tcell        Treg Macrophages       Bcell
#        1.0          NA         0.5         1.0
brain_cor_vec
#      Tcell        Treg Macrophages       Bcell
#       -1.0         0.0         1.0         0.5

## plot scatter plot
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 8
cols = gg_color_hue(n)
cols_sel <- cols[c(1, 3, 6, 7)]

groups <- factor(rep(paste0(names(ov_cor_vec), " cor=", 
  round(ov_cor_vec, digits=2)), each=3),
  level=c("Bcell cor=1", "Macrophages cor=0.5", "Tcell cor=1", "Treg cor=NA"))
ov_plot <- data.frame(ciber=c(ov_ciber_delta), stain=c(ov_stain_delta),
  groups=groups)

p1 <- ggplot(ov_plot, aes(x=ciber, y=stain, group=groups, color=groups)) +
  geom_point()+
  scale_color_manual(values=cols_sel)+
  labs(x="Percentage delta by CIBERSORT", y="Percentage delta by staining")+
  theme(legend.title=element_blank())

groups <- factor(rep(paste0(names(brain_cor_vec), " cor=", 
  round(brain_cor_vec, digits=2)), each=3),
  level=c("Bcell cor=0.5", "Macrophages cor=1", "Tcell cor=-1", "Treg cor=0"))
brain_plot <- data.frame(ciber=c(brain_ciber_delta), stain=c(brain_stain_delta),
  groups=groups)

p2 <- ggplot(brain_plot, aes(x=ciber, y=stain, group=groups, color=groups)) +
  geom_point()+
  scale_color_manual(values=cols_sel)+
  labs(x="Percentage delta by CIBERSORT", y="Percentage delta by staining")+
  theme(legend.title=element_blank())

pdf("Results_v2/OV_stain_ciber_corr.pdf", width=8.5, height=2.5)
grid.arrange(p1, p2, nrow=1)
dev.off()

pdf(file="Results_v2/OV_stain_ciber_corr.pdf", width=8, height=4)
par(mfrow=c(1,2), mar=c(5.1, 4.1, 0.5, 2.1), mgp=c(2.5, 1, 0), las=0)
plot(c(ov_ciber_delta), c(ov_stain_delta), 
  col=as.factor(rep(colnames(ov_ciber_delta), each=3)),
  pch=as.numeric(as.factor(rep(colnames(ov_ciber_delta), each=3))),
  ylab="Percentage delta by staining",
  xlab="Percentage delta by CIBERSORT",
  xlim=c(-0.25, 0.4), ylim=c(-0.5, 1))
legend_vec <- levels(as.factor(rep(colnames(ov_ciber_delta), each=3)))
legend_vec
ov_cor_vec[legend_vec]
legend_text <- c("Bcell cor=1", "Macrophages cor=0.5", "Tcell cor=1", "Treg cor=NA")
legend("topleft", legend=legend_text,
  col=1:4, pch=1:4, ncol=1)

plot(c(brain_ciber_delta), c(brain_stain_delta), 
  col=as.factor(rep(colnames(brain_ciber_delta), each=3)),
  pch=as.numeric(as.factor(rep(colnames(brain_ciber_delta), each=3))),
  ylab="Percentage delta by staining",
  xlab="Percentage delta by CIBERSORT",
   xlim=c(-0.25, 0.4), ylim=c(-0.3, 0.8))
legend_vec <- levels(as.factor(rep(colnames(brain_ciber_delta), each=3)))
brain_cor_vec[legend_vec]
legend_text <- c("Bcell cor=0.5", "Macrophages cor=1", "Tcell cor=-1", "Treg cor=0")
legend("topleft", legend=legend_text,
  col=1:4, pch=1:4, ncol=1)
dev.off()

















