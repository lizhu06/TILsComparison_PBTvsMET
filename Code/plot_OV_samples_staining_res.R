rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")
library(ggplot2)
library(reshape2)

res <- read.csv("Results_v2/RawData_Ovarian_staining_5_29_18.csv")
res[which(res[,"Slide.ID"] == "31-50 (primary)"), "Slide.ID"] <- "0031-50T"
res[which(res[,"Slide.ID"] == "6-56 (met)"), "Slide.ID"] <- "0006-56M"

pri_id <- c("0031-50T", "OP1", "OP3")
mets_id <- c("0006-56M", "OM1", "OM3")
pair_id <- c("31-50T/6-56M", "OP1/OM1", "OP3/OM3")
id <- c(pri_id, mets_id)
uni_pair_id <- rep(pair_id, 2)
uni_group_id <- c(rep("PRI", 3), rep("MET", 3))
pair_id <- uni_pair_id[match(res[,"Slide.ID"], id)]
group_id <- uni_group_id[match(res[,"Slide.ID"], id)]

# convert stroma and tumor percentage
stroma_perc <- as.numeric(sapply(1:nrow(res), function(x) 
	gsub("%", "", res[x, "Stroma.Area"])))/100
tumor_perc <- as.numeric(sapply(1:nrow(res), function(x) 
	gsub("%", "", res[x, "Tumor.Area"])))/100
# normalize to 1
total_perc <- stroma_perc + tumor_perc
stroma_perc <- stroma_perc/total_perc
tumor_perc <- tumor_perc/total_perc
#tumor_perc[tumor_perc==0] <- 0.0001

immune_perc <- res[, 6:15]		
for(i in 1:nrow(immune_perc)){
	for(j in 1:ncol(immune_perc)){
		immune_perc[i,j] <- as.numeric(gsub("%", "", immune_perc[i,j]))/100
	}
}
immune_perc <- data.matrix(immune_perc)

# get cell number
stroma_tumor_perc <- cbind(stroma_perc, tumor_perc)[, rep(c(1,2), times=5)]
immune_num <- sweep(immune_perc * stroma_tumor_perc, 
	1, res[,"X..of.Cells"], "*")
immune_num <- immune_num[, -c(7, 8)]

# calculate combined cell number
comb_num <- matrix(NA, nrow(res), 4)
for(i in 1:4){
	comb_num[, i] <- immune_num[, (i-1)*2+1] + 
		immune_num[, (i-1)*2+2]
}
colnames(comb_num) <- c("CD8", "FOXP3", "CD68", "CD20")

total_immune <- apply(comb_num, 1, sum)
comb_prop_rel <- sweep(comb_num, 1, total_immune, "/")
# remove one slide doesn't have immune 
rm_index <- which(apply(is.na(comb_prop_rel),1, max)==1)
comb_prop_rel <- comb_prop_rel[-rm_index, ] 
pair_id <- pair_id[-rm_index]
group_id <- group_id[-rm_index]
dim(comb_prop_rel) #50  4

# merge to long data
score <- cbind(pair_id, group_id, res[-rm_index,5], comb_prop_rel)
colnames(score)[3] <- "Slide.ID"
score <- as.data.frame(score)
save(score, file="Results_v2/OV_mets_staining_score.RData")
score2 <- melt(score, id.vars=c("pair_id", "group_id","Slide.ID"))
colnames(score2)[c(4,5)] <- c("cell", "percentage")
score2[,"percentage"] <-  as.numeric(score2[,"percentage"])*100
uni_cell_type <- colnames(comb_num)

# calculate average and SD
uni_cell_type <- as.character(unique(score2[,"cell"]))
sum_score2 <- matrix(NA, length(id)*length(uni_cell_type), 6)
temp <- 0
for(i in 1:length(id)){
	for(j in 1:length(uni_cell_type)){
		temp <- temp + 1
		vec <- score2[which(score2[,"Slide.ID"]==id[i] & 
			score2[,"cell"] == uni_cell_type[j]),"percentage"]
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
p <- ggplot(df, aes(x=Group, y=mean, group=Marker, color=Marker)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                 position=position_dodge(0.05))+
  scale_color_manual(values=cols_sel) +
  theme(text = element_text(size=10))
p <- p + labs(y="percent", x="")
p <- p + facet_grid(. ~ Pair)

pdf(file="Results_v2/OV_staining_res_regions_combined_rel_totalImmune_trash.pdf", 
	width=4, height=2)
print(p)
dev.off()





