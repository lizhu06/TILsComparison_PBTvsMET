rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")
library(ggplot2)
library(reshape2)

## load  scores
load("Results_v2/gsva_score_davoli_log2tpm.RData")
load("Results_v2/gsva_score_tamborero_log2tpm.RData")
load("Results_v2/cibersort_relative.RData")
load("Results_v2/cibersort_absolute.RData")
load("Results_v2/TIMER_res_pri_mets.RData")
load("Results_v2/log2tpm_estimate_purity_score.RData")
rownames(ciber_absolute) <- gsub("T.cells.regulatory..Tregs.", 
	"T.cells.regulatory", rownames(ciber_absolute))
rownames(ciber_relative) <- gsub("T.cells.regulatory..Tregs.", 
	"T.cells.regulatory", rownames(ciber_relative))

pri_id <- c("0031-50T", "OP1", "OP3")
mets_id <- c("0006-56M", "OM1", "OM3")
id <- c(pri_id, mets_id)
pair_id <- rep(c("31-50T/6-56M", "OP1/OM1", "OP3/OM3"), 2)
group_id <- c(rep("PRI", 3), rep("MET", 3))

### set colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 8
cols = gg_color_hue(n)
#names(cols) <- c("Bcell1", "Bcell2", "macro0", "macro1", "macro2",
#	"cd8", "treg", "tumor")

ciber_cell_names <- c("B.cells.naive", "B.cells.memory", 
	 "Macrophages.M0", "Macrophages.M1", 
	"Macrophages.M2", "T.cells.CD8", "T.cells.regulatory")

ciber_rel_sub <- rbind(
	apply(ciber_relative[
		c("B.cells.naive", "B.cells.memory"), id],2, sum),
	apply(ciber_relative[
		c("Macrophages.M0", "Macrophages.M1", 
			"Macrophages.M2"), id],2, sum),
	ciber_relative[c("T.cells.CD8", "T.cells.regulatory"), id])
rownames(ciber_rel_sub) <- c("Bcell", "Macrophages", "Tcell", "Treg")	
save(ciber_rel_sub, file="Results_v2/OV_mets_stainedSample_bioinfo_score.RData")
## transform measurements
#purity <- estimate_purity_score["TumorPurity",id]
#stromal_score <- estimate_purity_score["StromalScore",id]

## plot
all_res <- list(
	ciber_rel_sub
)
names(all_res) <- c("ciber_rel")
ylab_tex <- c("percent")
cols_sel <- cols[c(1, 3, 6, 7)]

plot_data <- function(data_index){
	if(FALSE){
		data_index <- 1
	}
	data <- cbind(pair_id, group_id, t(all_res[[data_index]]))
	data <- cbind(rownames(data), data)
	colnames(data)[1] <- "ID"
	rownames(data) <- NULL
	data2 <- melt(data.frame(data), id.vars=c("ID", "pair_id", "group_id"))
	colnames(data2)[c(4,5)] <- c("Cell", "score")
	data3 <- data.frame(ID=factor(data2$ID), Pair=factor(data2$pair_id),
		Group=factor(data2$group_id, levels=c("PRI","MET")), Cell=data2$Cell, 
		score=as.numeric(data2$score)*100)
	data3$Cell <- factor(data3$Cell, 
		levels=c("Bcell", "Macrophages", "Tcell", "Treg"))

	p <- ggplot(data3, aes(x=Group, y=score, group=Cell, color=Cell)) + 
		geom_line() +
	  geom_point()+
	  scale_color_manual(values=cols_sel)+
	  theme(text = element_text(size=10))
	p <- p + labs(y=ylab_tex[data_index], x="")
	p <- p + facet_grid(. ~ Pair)

	pdf(file=paste0("Results_v2/OV_staining_samples_", 
		names(all_res)[data_index], "_v2.pdf"), 
		width=4.5, height=2)
	print(p)
	dev.off()
	return(p)
}
rs <- lapply(1:length(all_res), function(x) plot_data(x))


