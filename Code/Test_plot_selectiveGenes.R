rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")
library(tidyr)
library(reshape2)
library(ggplot2)
library(multcomp)

## load  scores
load("Data_v2/log2tpm_unique.RData")
log2tpm <- log2tpm_unique
rm(log2tpm_unique)
load("Data_v2/gene_annot_biomart_unique.RData")
rownames(log2tpm) <- gene_annot_biomart_unique[,"external_gene_name_v2"]

## load sample data
load("Data_v2/sample_pair_info.RData")
load("Data_v2/sample_annot.RData")
load("Results_v2/immune_res_aveDup.RData")

sample_pair_info <- sample_pair_info[
  which(sample_pair_info[,"mets_id"] != "7M_RCS"), ]
sample_annot <- sample_annot[which(
  !(sample_annot[,"ID"] %in% c("7M_RCS","7P_RCS"))), ]

selective_genes <- c("CD274", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "VSIR")
all(selective_genes %in% rownames(log2tpm))
#rownames(log2tpm)[grepl("HAVCR2", rownames(log2tpm))]

sample_pair_info2 <- sample_pair_info
for(i in 1:length(selective_genes)){

  sample_pair_info2$primary_gene <- log2tpm[selective_genes[i],
    sample_pair_info2$primary_id]
  sample_pair_info2$mets_gene <- log2tpm[selective_genes[i],
    sample_pair_info2$mets_id]
  sample_pair_info2$gene_change <- sample_pair_info2$mets_gene - 
    sample_pair_info2$primary_gene

  colnames(sample_pair_info2)[(ncol(sample_pair_info2)-2):
    ncol(sample_pair_info2)] <- paste0(c("primary_", "mets_", "change_"), 
      selective_genes[i])
}

sample_pair_info2$primary_immune <- immune_res_aveDup[[1]][sample_pair_info2$primary_id,"ImmuneScore"]
sample_pair_info2$mets_immune <- immune_res_aveDup[[1]][sample_pair_info2$mets_id,"ImmuneScore"]
sample_pair_info2$change_immune <- sample_pair_info2$mets_immune - sample_pair_info2$primary_immune

sample_pair_info2$primary_purity <- immune_res_aveDup[[1]][sample_pair_info2$primary_id,"TumorPurity"]
sample_pair_info2$mets_purity <- immune_res_aveDup[[1]][sample_pair_info2$mets_id,"TumorPurity"]
sample_pair_info2$change_purity <- sample_pair_info2$mets_purity - sample_pair_info2$primary_purity


## check outliers
if(FALSE){
  change_col <- c("pairs","change_CD274","change_PDCD1", "change_CTLA4","change_immune")
  summary(sample_pair_info2[,c("change_CD274","change_PDCD1", "change_CTLA4","change_immune")])
  sample_pair_info2[which(sample_pair_info2[,"change_CD274"]>0.5), change_col]
  sample_pair_info2[which(sample_pair_info2[,"change_PDCD1"]>0.5), change_col]
  sample_pair_info2[which(sample_pair_info2[,"change_CTLA4"]>0.5), change_col]
  sample_pair_info2[which(sample_pair_info2[,"change_immune"]>500), change_col]
  cor(sample_pair_info2[,c("change_CD274","change_PDCD1", "change_CTLA4","change_immune")], method="spearman")
  cor.test(sample_pair_info2[,"change_CD274"],
    sample_pair_info2[,"change_immune"], method="spearman")
  cor.test(sample_pair_info2[,"change_PDCD1"],
    sample_pair_info2[,"change_immune"], method="spearman")
  cor.test(sample_pair_info2[,"change_CTLA4"],
    sample_pair_info2[,"change_immune"], method="spearman")

  # select a few interesting pairs
  interesting_pairs <- c("Bo43", "GiGM1", "Bo34")
  interesting_pairs <- sample_pair_info2[match(interesting_pairs, 
    sample_pair_info2[,"pairs"]),]
}

# load survival 
load("Data_v2/MFS.RData")
load("Data_v2/SPM.RData")


# for multiple mets matched to same primary, take average 
#     and remove the duplicated mets 
dup_pri_id <- sample_pair_info2[duplicated(
  sample_pair_info2$primary_id), "primary_id"]
for(du in 1:length(dup_pri_id)){
  row_index <- which(sample_pair_info2[,"primary_id"] == dup_pri_id[du])
  col_index <- (ncol(sample_pair_info2)-length(selective_genes)*3+1):
    ncol(sample_pair_info2)
  mean_exp <- apply(sample_pair_info2[row_index, col_index], 2, mean)
  sample_pair_info2[row_index[1], col_index] <- mean_exp
  sample_pair_info2 <- sample_pair_info2[-row_index[-1],]
}
dim(sample_pair_info2) #50 21

# set color
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n <- 6
ggcols = gg_color_hue(n)

########################################################
## plot primary and mets without normalization
########################################################
pval <- matrix(NA, length(selective_genes), 2)
rownames(pval) <- selective_genes
colnames(pval) <- c("median", "pvalue")
for(i in 1:length(selective_genes)){
  # plot primary and mets separately
  data_plot <- sample_pair_info2[, 
    c("pairs", 
      paste0(c("primary_", "mets_"), selective_genes[i]))]
  colnames(data_plot)[c(2,3)] <- c("PBT", "MET")
  data_long <- melt(data_plot, id.vars=c("pairs")) 
  colnames(data_long)[2] <- "Group"

  # normalization


  #data_long$Group <- factor(data_long$Group, 
  #  levels=c("PBT", site_name_vec[i], "Change"))

  p <- ggplot(data_long, aes(x=Group, y=value, fill=Group)) +
    geom_point(size=2, alpha=0.6, aes(color=Group))
    
  #p <- p + geom_line(aes(group = pairs), alpha = 0.6, colour = "grey", data = data_long)

  p <- p + labs(x="", y=selective_genes[i]) + 
    scale_colour_manual(values = c(ggcols[c(4,1)]))

  # add segmentation to connect pairs  
  for(t in 1:nrow(data_plot)){
    p <- p + geom_segment(x=1, y=data_plot[t, 2], 
        xend=2, yend=data_plot[t, 3], colour="grey") 
  }
  p <- p + labs(y=paste0(selective_genes[i])) 
  #p <- p + ylim(rangee[1]-100, rangee[2]+100)
  pdf(file=paste0("Results_v2/boxplot_selective_genes/boxplot_genes_", selective_genes[i], 
    "_primary_mets.pdf"), 
    width=3, height=2.5)
  plot(p)
  dev.off()

  # test p-value 
  fit <- wilcox.test(data_plot[,2], 
    data_plot[,3], paired=TRUE)
  pval[i,2] <- fit$p.value

  # plot delta
  data_plot <- sample_pair_info2[,  c("pairs", 
      paste0(c("change_"), selective_genes[i]))]
  colnames(data_plot)[2] <- c("Change")
  pval[i,1] <- median(data_plot[,"Change"])
  p_delta <- ggplot(data_plot, aes(x="", y=Change)) +
    geom_boxplot() + geom_jitter(position=position_jitter(0.2))
  #p <- p + geom_line(aes(group = pairs), alpha = 0.6, colour = "grey", data = data_long) 
  p_delta <- p_delta + labs(x="", y=paste0(selective_genes[i], " Change")) 

  pdf(file=paste0("Results_v2/boxplot_selective_genes/boxplot_genes_", 
    selective_genes[i], 
    "_change.pdf"), width=2, height=2.5)
  plot(p_delta)
  dev.off()
}

#           median       pvalue
#CD274  -0.5057083 1.649725e-03
#PDCD1  -0.2696317 4.018867e-03
#CTLA4  -0.9305997 4.229706e-07
#LAG3   -0.1779306 1.296294e-01
#HAVCR2 -0.5348288 5.765332e-03
#VSIR   -0.1380703 5.495042e-01


log2tpm["LAG3", c("19P", "19M")]



########################################################
## plot primary and mets with normalization
########################################################
pval_norm <- matrix(NA, length(selective_genes), 2)
rownames(pval_norm) <- selective_genes
colnames(pval_norm) <- c("median", "pvalue")
for(i in 1:length(selective_genes)){
  # plot primary and mets separately
  data_plot <- sample_pair_info2[, 
    c("pairs", 
      paste0(c("primary_", "mets_"), selective_genes[i]), 
      "primary_purity", "mets_purity")]
  colnames(data_plot)[c(2,3)] <- c("PBT", "MET")

  # normalization
  data_plot[,"PBT"] <- data_plot[,"PBT"]/(1-data_plot[,"primary_purity"])
  data_plot[,"MET"] <- data_plot[,"MET"]/(1-data_plot[,"mets_purity"])
  data_plot <- data_plot[, -c(4,5)]

  data_long <- melt(data_plot, id.vars=c("pairs")) 
  colnames(data_long)[2] <- "Group"
  #data_long$Group <- factor(data_long$Group, 
  #  levels=c("PBT", site_name_vec[i], "Change"))

  p <- ggplot(data_long, aes(x=Group, y=value, fill=Group)) +
    geom_point(size=2, alpha=0.6, aes(color=Group))
    
  #p <- p + geom_line(aes(group = pairs), alpha = 0.6, colour = "grey", data = data_long)

  p <- p + labs(x="", y=selective_genes[i]) + 
    scale_colour_manual(values = c(ggcols[c(4,1)]))

  # add segmentation to connect pairs  
  for(t in 1:nrow(data_plot)){
    p <- p + geom_segment(x=1, y=data_plot[t, 2], 
        xend=2, yend=data_plot[t, 3], colour="grey") 
  }
  p <- p + labs(y=paste0(selective_genes[i])) 
  #p <- p + ylim(rangee[1]-100, rangee[2]+100)
  pdf(file=paste0("Results_v2/boxplot_selective_genes_normalized/boxplot_genes_", selective_genes[i], 
    "_primary_mets.pdf"), 
    width=3, height=2.5)
  plot(p)
  dev.off()

  # test p-value 
  fit <- wilcox.test(data_plot[,2], 
    data_plot[,3], paired=TRUE)
  pval_norm[i, 2] <- fit$p.value

  # plot delta
  data_plot <- data.frame(pairs=data_plot[,"pairs"], Change=data_plot[,"MET"] - data_plot[,"PBT"])
  pval_norm[i, 1] <- median(data_plot[,"Change"])
  p_delta <- ggplot(data_plot, aes(x="", y=Change)) +
    geom_boxplot() + geom_jitter(position=position_jitter(0.2))
  #p <- p + geom_line(aes(group = pairs), alpha = 0.6, colour = "grey", data = data_long) 
  p_delta <- p_delta + labs(x="", y=paste0(selective_genes[i], " Change")) 

  pdf(file=paste0("Results_v2/boxplot_selective_genes_normalized/boxplot_genes_", 
    selective_genes[i], 
    "_change.pdf"), width=2, height=2.5)
  plot(p_delta)
  dev.off()
}
pval_norm
#            median       pvalue
#CD274   1.16359042 5.765332e-03
#PDCD1  -0.03260769 8.848694e-01
#CTLA4  -1.61888789 7.496180e-03
#LAG3    0.83288347 1.850278e-02
#HAVCR2  2.70669955 2.351629e-04
#VSIR    3.25593385 1.121878e-05


################
# Site specific
#################
## plot primary and mets
color_index <- c(2,3,5,6)
site_vec <- c("Brain", "Ovary", "Bone", "GI")
site_name_vec <- c("BRM", "OVM", "BOM", "GIM")
#i <- 1
pval <- matrix(NA, 4, length(selective_genes))

for(j in 1:length(selective_genes)){
  for(i in 1:4){
    data_plot <- sample_pair_info2[which(sample_pair_info2$mets_site==site_vec[i]), 
      c("pairs", "mets_site", 
        paste0(c("primary_", "mets_","change_"), selective_genes[j]))]

    colnames(data_plot)[c(3,4,5)] <- c("PBT", site_name_vec[i], "Change")

    rangee <- range(data_plot[, c("PBT", site_name_vec[i], "Change")])

    data_long <- melt(data_plot, id.vars=c("pairs", "mets_site")) 
    colnames(data_long)[3] <- "Group"
    data_long$Group <- factor(data_long$Group, levels=c("PBT", 
      site_name_vec[i], "Change"))

    p <- ggplot(data_long, aes(x=Group, y=value, fill=Group)) +
      geom_point(size=2, alpha=0.6, aes(color=Group))
      
    #p <- p + geom_line(aes(group = pairs), alpha = 0.6, colour = "grey", data = data_long)

    p <- p + labs(x="", y=paste0(selective_genes[j])) + 
      scale_colour_manual(values = c(ggcols[c(4,color_index[i])], "black"))
    for(t in 1:nrow(data_plot)){
      p <- p + geom_segment(x=1, y=data_plot[t, 3], 
          xend=2, yend=data_plot[t, 4], colour="grey") 
    }
    p <- p + ylim(rangee[1]-1, rangee[2]+1)

    pdf(file=paste0("Results_v2/boxplot_selective_genes/boxplot_genes_", 
      selective_genes[j], "_",
      site_vec[i], ".pdf"), 
      width=3, height=2.5)
    plot(p)
    dev.off()

    fit <- wilcox.test(data_plot[,3], 
      data_plot[,4], paired=TRUE)
    pval[i,j] <- fit$p.value
  }
}

