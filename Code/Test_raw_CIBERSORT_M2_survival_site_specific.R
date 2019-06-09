rm(list=ls())
options(stringsAsFactors = FALSE)
library(survival)
#library(rms)
library(gridExtra)
setwd("/net/wong05/home/liz86/Steffi/primary_vs_mets/")
library("survminer")

#### load data
load("Data_v2/sample_annot.RData")
load("Data_v2/sample_pair_info.RData")
load("Results_v2/immune_res_aveDup.RData")

load("Data_v2/all_clinical_info.RData")
load("Data_v2/MFS.RData")
load("Data_v2/SPM.RData")

surv_list <- list(MFS, SPM)
names(surv_list) <- c("MFS", "SPM")

# get immune
immune <- immune_res_aveDup
check_subj_id <- sapply(1:length(immune), function(x) 
	all(rownames(immune[[1]])==rownames(immune[[x]])))
immune[[4]] <- immune[[4]][,"Macrophages.M2", drop=FALSE]
immune <- immune[-5]

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

ggcols = gg_color_hue(2)

### plot KM for all patients
width <- 10
height_vec <- c(4, 50, 80, 4, 30)
plot_km <- 1
site_vec <- c("Bone", "Brain", "Ovary", "GI")

for(si in 1:length(site_vec)){
	for(su in 1:length(surv_list)){
		surv_obj <- surv_list[[su]]
		surv_name <- names(surv_list)[su]

		#for(t in 1:length(immune)){
			t <- 4
			data <- immune[[t]]
			var_names <- colnames(data)
			fileName <- paste0("raw_immune_survival/Raw_", 
				names(immune)[t], "_macrophage_m2_", surv_name, "_", site_vec[si])		

			plot_index <- 0
			choice_vec <- c("primary", "mets", "delta")
			pval_mat <- matrix(NA, ncol(data), 3)
			rownames(pval_mat) <- colnames(data)
			colnames(pval_mat) <- choice_vec

			pl <- list()
			for(c in 1:length(choice_vec)){
				if(choice_vec[c]=="primary"){
					id <- rownames(surv_obj)
					score <- data[id, , drop=FALSE]
					mets_site <- sample_pair_info[match(id, 
						sample_pair_info$primary_id), "mets_site"]  
				}else if(choice_vec[c]=="mets"){
					id <- sample_pair_info[match(rownames(surv_obj),
						sample_pair_info[,"primary_id"]), "mets_id"]
					score <- data[id, , drop=FALSE]
					mets_site <- sample_pair_info[match(id, 
						sample_pair_info$mets_id), "mets_site"] 
				}else{
					primary_id <- rownames(surv_obj)
					mets_id <- sample_pair_info[match(primary_id,
							sample_pair_info[,"primary_id"]), "mets_id"]
					score <- data[mets_id, , drop=FALSE] - 
						data[primary_id, , drop=FALSE]
					mets_site <- sample_pair_info[match(mets_id, 
						sample_pair_info$mets_id), "mets_site"] 
				}
				keep_index <- which(mets_site==site_vec[si])
				score <- score[keep_index, , drop=FALSE]
				surv_obj <- surv_obj[keep_index, ]

				for(i in 1:ncol(score)){
					    
			    value <- score[,i, drop=FALSE]
			    
			    # remove missing
			    miss_index <- is.na(value)
			    value <- value[!miss_index]
			    value <- as.factor((value > median(value))*1)
			    value_factor <- rep(NA, length(value))
			    value_factor[value==0] <- "<=median"
			    value_factor[value==1] <- ">median"
			    value <- factor(value_factor, levels=c("<=median", 
			    	">median"))
			    surv_obj <- surv_obj[!miss_index, ]
			    
			    if(length(unique(value))>1){
			    	# log-rank test
			    	fit <- survdiff(surv_obj ~ value)
			    	p1 <- round(1 - pchisq(fit$chisq, 
			    		length(levels(value))-1), digits=4)
			    	if(p1==0){
			    	  p1 <- "<0.0001"
			    	}
			    	pval_mat[i, c] <- p1
			    	fit <- survfit(surv_obj~ value)
			    	
			    	survdata <- data.frame(time=surv_obj[,1], 
			    		status=surv_obj[, 2], 
			    	                       group=value)
			    	if(plot_km){
			    		plot_index <- plot_index + 1
			    	  pl[[plot_index]] <- ggsurvplot(
			    	    fit,                     # survfit object with calculated statistics.
			    	    data = survdata,             # data used to fit survival curves.
			    	    risk.table = TRUE,       # show risk table.
			    	    fontsize=4,
			    	    pval = TRUE,             # show p-value of log-rank test.
			    	    conf.int = FALSE,         # show confidence intervals for 
			    	    # point estimates of survival curves.
			    	    xlim = c(0,max(surv_obj[,1], na.rm=TRUE)+5),         # present narrower X axis, but not affect
			    	    # survival estimates.
			    	    xlab = "Month",   # customize X axis label.
			    	    ylab = surv_name,
			    	    break.time.by = 30,     # break X axis in time intervals by 500.
			    	    ggtheme = theme_light(), # customize plot and risk table with a theme.
			    	    risk.table.y.text.col = TRUE, # colour risk table text annotations.
			    	    risk.table.y.text = FALSE, # show bars instead of names in text annotations
			    	    # in legend of risk table                   
			    	    legend.labs = paste0(var_names[i], " ",levels(value)),
			    	    font.legend=8,
			    	    #pval.method = TRUE,
			    	    title=choice_vec[c], 
			    	    palette=ggcols[c(2,1)]                  
			    	  )
			    	}
			    }
				}
			}
		  if(plot_km){
		    pdf(file=paste0("Results_v2/", fileName, ".pdf"), width, height=height_vec[t])
		    arrange_ggsurvplots(pl, print = TRUE, 
		                        title = NA, ncol = 3, nrow = ncol(data), 
		                        surv.plot.height=0.7, risk.table.height =0.25)
		    dev.off()
		  }
			write.csv(pval_mat, file=paste0("Results_v2/", fileName, ".csv"))
		#}
	}
}



