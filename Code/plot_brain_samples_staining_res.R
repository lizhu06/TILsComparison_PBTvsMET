rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")
library(ggplot2)
library(reshape2)

res <- read.csv("Results_v2/RawData_BrainMets_8.13.18.csv")
res <- res[1:51,1:11]

pri_id <- c("BP52", "BP72", "BP71")
mets_id <- c("BM52", "BM72", "BM71")
pair_id <- c("BP/BM52", "BP/BM72", "BP/BM71")
id <- c(pri_id, mets_id)
uni_pair_id <- rep(pair_id, 2)
uni_group_id <- c(rep("PRI", 3), rep("MET", 3))
pair_id <- uni_pair_id[match(res[,"Slide.ID"], id)] # expressed as mets name
group_id <- uni_group_id[match(res[,"Slide.ID"], id)] # primary or mets

immune_perc <- res[, c("CD20", "CD8", "CD68", "FOXP3")]
immune_perc <- do.call(cbind, lapply(1:ncol(immune_perc), function(x) 
	as.numeric(gsub("%", "", immune_perc[,x]))))/100
immune_num <- sweep(immune_perc, 1, res[,"Number.of.Cells"], "*")
total_immune <- apply(immune_num, 1, sum)
immune_prop_rel <- sweep(immune_num, 1, total_immune, "/")

score <- as.data.frame(cbind(pair_id, group_id, res[,3], immune_prop_rel))
colnames(score) <- c("pair_id", "group_id", "Slide.ID",
 "CD20", "CD8", "CD68", "FOXP3")
save(score, file="Results_v2/brain_mets_staining_score.RData")
score2 <- melt(score, id.vars=c("pair_id", "group_id","Slide.ID"))
colnames(score2)[c(4,5)] <- c("cell", "percent")


# calculate average and SD (and convert to long data)
uni_cell_type <- as.character(unique(score2[,"cell"]))
sum_score2 <- matrix(NA, length(id)*length(uni_cell_type), 6)
temp <- 0
for(i in 1:length(id)){
	for(j in 1:length(uni_cell_type)){
		temp <- temp + 1
		vec <- as.numeric(score2[which(score2[,"Slide.ID"]==id[i] & 
			score2[,"cell"] == uni_cell_type[j]),"percent"])*100 # convert to percent
		sum_score2[temp, 1] <- id[i]
		sum_score2[temp, 2] <- uni_pair_id[i]
		sum_score2[temp, 3] <- uni_group_id[i]
		sum_score2[temp, 4] <- uni_cell_type[j]
		sum_score2[temp, 5] <- mean(vec)
		sum_score2[temp, 6] <- sd(vec)/sqrt(length(vec))
		print(j)
	}
}
df <- data.frame(ID=as.factor(sum_score2[,1]), 
	Pair=as.factor(sum_score2[,2]), 
	Group=factor(sum_score2[,3], levels=c("PRI","MET")),
	Marker=as.factor(sum_score2[,4]), mean=as.numeric(sum_score2[,5]),
		sd=as.numeric(sum_score2[,6]))

#df$Region <- factor(rep("tumor", nrow(df)), levels=c("tumor","stroma"))
#df$Region[sapply(1:length(df$MarkerLong), function(x) 
#	grepl("stroma", df$MarkerLong[x]))] <- "stroma"

#df$Marker <- factor(sapply(1:length(df$MarkerLong), function(x) 
#	strsplit(as.character(df$MarkerLong[x]), split="\\.")[[1]][1]))

levels(df$Marker)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 8
cols = gg_color_hue(n)
#names(cols) <- c("Bcell1", "Bcell2", "macro0", "macro2", "macro3",
#	"cd8", "treg", "tumor")

cols_sel <- cols[c(1, 3, 6, 7)]

## plot
p <- ggplot(df, aes(x=Group, y=mean, color=Marker, group=Marker)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                 position=position_dodge(0.05))+
  scale_color_manual(values=cols_sel) +
  theme(text = element_text(size=10))
p <- p + labs(y="percent", x="")
p <- p + facet_grid(. ~ Pair)

pdf(file="Results_v2/Brain_staining_res_rel_totalImmune.pdf", 
	width=4, height=2)
print(p)
dev.off()





