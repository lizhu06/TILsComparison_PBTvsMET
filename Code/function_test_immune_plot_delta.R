gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n <- 6
ggcols = gg_color_hue(n)

paired_test <- function(score, fileName, ylab, num_row=3, 
	width=9, height=5, titleFont=4){

	## test for p-value and qvalue
	test_immune_site_combined <- matrix(NA, ncol(score), 3)
	for(i in 1:ncol(score)){
		mets_score <- as.numeric(score[mets_id, i])
		primary_score <- as.numeric(score[pri_id, i])
		test_immune_site_combined[i, 1] <- mean(mets_score - primary_score)
		test_immune_site_combined[i, 2] <- wilcox.test(primary_score, mets_score, 
			paired=T, conf.int=T)$p.value
	}
	test_immune_site_combined[,3] <- p.adjust(
		test_immune_site_combined[,2], method="fdr")
	rownames(test_immune_site_combined) <- colnames(score)
	colnames(test_immune_site_combined) <- c("ave delta", "pval", "qval")
	write.csv(test_immune_site_combined, 
		file=paste0(fileName, ".csv"))

	# plot boxplot of delta
	delta_score <- score[mets_id, ] - score[pri_id, ]
	delta_score <- as.data.frame(delta_score)
	delta_score$mets_id <- mets_id

	data_plot <- melt(delta_score, 
		variable.name="cellType", value.name="score", 
		id.vars=c("mets_id"))
	data_plot$cellType <- factor(as.character(data_plot$cellType))

	p <- ggplot(data_plot, aes(x=factor(0), y=score)) + 
  		geom_boxplot() + 
  		theme(axis.title.x=element_blank(),
  		   axis.text.x=element_blank(),
  		   axis.ticks.x=element_blank(),
         strip.text = element_text(size=titleFont)) + 
  		#geom_point(aes(fill=cellType), size=5, shape=21, 
  		#	position=position_jitterdodge()) + 
  		facet_wrap(~cellType, scale="free_y", nrow=num_row) +
      labs(x="", y=ylab)

  pdf(file=paste0(fileName, ".pdf"), width=width, height=height)
  plot(p)
  dev.off()		
}

paired_test_one_box <- function(score, fileName, 
  width=9, height=5){

  ## test for p-value and qvalue
  test_immune_site_combined <- matrix(NA, ncol(score), 3)
  for(i in 1:ncol(score)){
    mets_score <- as.numeric(score[mets_id, i])
    primary_score <- as.numeric(score[pri_id, i])
    test_immune_site_combined[i, 1] <- mean(mets_score - primary_score)
    test_immune_site_combined[i, 2] <- wilcox.test(primary_score, mets_score, 
      paired=T, conf.int=T)$p.value
  }
  test_immune_site_combined[,3] <- p.adjust(
    test_immune_site_combined[,2], method="fdr")
  rownames(test_immune_site_combined) <- colnames(score)
  colnames(test_immune_site_combined) <- c("ave delta", "pval", "qval")
  write.csv(test_immune_site_combined, 
    file=paste0(fileName, ".csv"))

  # plot boxplot of delta
  delta_score <- score[mets_id, ] - score[pri_id, ]
  delta_score <- as.data.frame(delta_score)
  delta_score$mets_id <- mets_id

  data_plot <- melt(delta_score, 
    variable.name="cellType", value.name="score", 
    id.vars=c("mets_id"))
  #data_plot$cellType <- factor(as.character(data_plot$cellType))

  p <- ggplot(data_plot, aes(x=cellType, y=score)) + 
      geom_boxplot() + 
      theme(axis.text.x=element_text(angle=45,hjust=1))+
#      theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank()) + 
      #geom_point(aes(fill=cellType), size=5, shape=21, 
      # position=position_jitterdodge()) + 
      labs(x="", y="")

  pdf(file=paste0(fileName, ".pdf"), width=width, height=height)
  plot(p)
  dev.off()   
}

### tissue specific 
paired_test_tissue_sep_one_box <- function(score, fileName, width, height){
	uni_site_color <- ggcols[c(2,3,5,6)]
  res <- matrix(NA, ncol(score), 21)
	# get pvalue table
  for(i in 1:ncol(score)){
    delta <- (score[mets_id, i] - score[pri_id, i])
    fit <- lm(delta ~ mets_site) 
    sum_fit <- summary(fit)$coefficients
    # slope for tumor purity
    #res[i, 1] <- sum_fit[2,1]
    #res[i, 2] <- sum_fit[2,4]

    #brain
    res[i, 1] <- sum_fit[1,1]
    res[i, 2] <- sum_fit[1,2]
    res[i, 3] <- sum_fit[1,4]
    
    # other site
    for(j in 1:3){
      K <- matrix(rep(0, nrow(sum_fit)), nrow=1)
      K[1, 1] <- K[1, j+1] <- 1
      lc_sum <- summary(glht(fit, linfct = K))
      res[i, 3+(j-1)*3+1] <- lc_sum$test$coefficients
      #res[i, 4+(j-3)*2+2] <- round(lincom[[6]][2], digits=2)
      res[i, 3+(j-1)*3+2] <- lc_sum$test$sigma
      res[i, 3+(j-1)*3+3] <- lc_sum$test$pvalue[1]
    }

    # delta
    for(j in 1:3){
      res[i, 12+(j-1)*3+1] <- sum_fit[j+1, 1]
      res[i, 12+(j-1)*3+2] <- sum_fit[j+1, 2]
      res[i, 12+(j-1)*3+3] <- sum_fit[j+1, 4]
    }
  }

  colnames(res) <- c("brain_est","brain_se", "brain_qval",
    "OV_est","OV_se", "OV_qval",
    "bone_est","bone_se", "bone_qval", 
    "GI_est", "GI_se", "GI_qval",  
    "OV_delta_est","OV_delta_se", "OV_delta_qval",
    "bone_delta_est","bone_delta_se", "bone_delta_qval", 
    "GI_delta_est", "GI_delta_se", "GI_delta_qval")
  rownames(res) <- colnames(score)

  # multiple test correction
  p_val <- res[, seq(3, 21, by=3)]
  p_adj <- matrix(p.adjust(p_val, method="fdr"), ncol=7, byrow=FALSE)
  res[, seq(3, 21, by=3)] <- p_adj
	write.csv(res, file=paste0(fileName, ".csv"))

	# boxplot
	delta <- score[mets_id, ] - score[pri_id, ]
	delta <- as.data.frame(delta)
	delta$mets_id <- mets_id
	delta$site <- mets_site
	data_plot <- melt(delta, 
		variable.name="cellType", value.name="delta_score", 
		id.vars=c("mets_id", "site"))
	data_plot$site <- factor(data_plot$site, levels=uni_mets_site)

	p <- ggplot(data_plot, aes(x=cellType, y=delta_score)) + 
  		geom_boxplot(aes(fill=site)) + 
      theme(axis.text.x=element_text(angle=45,hjust=1))+
  		#geom_point(aes(fill=cellType), size=5, shape=21, 
  		#	position=position_jitterdodge()) + 
  		#facet_wrap(~cellType, scale="free", nrow=num_row)+
  		scale_fill_manual(values=ggcols[c(2,3,5,6)]) +
      labs(x="", y="")
  pdf(file=paste0(fileName, ".pdf"), width=width, height=height)
  plot(p)
  dev.off()	
}

### HR and HER2
paired_test_HR_HER2 <- function(score, fileName,
  width, height, onlyBrain=FALSE){
  
  ggcols = gg_color_hue(11)
  uni_HR_HER2_color <- ggcols[c(4,7,2,9,3)]

  if(!onlyBrain){
    uni_HR_HER2_p <- uni_HR_HER2
    }else{
      uni_HR_HER2_p <- uni_HR_HER2_brain
    }
  
  if(!onlyBrain){
    keep_index <- seq(1, length(pri_id))
  }else{
    keep_index <- which(
      sample_pair_info$mets_site=="Brain")
  }

  # get pvalue table
  num_other_levels <- length(uni_HR_HER2_p)-1
  res <- matrix(NA, ncol(score), (num_other_levels+1)*4+
    num_other_levels*4)
  # get pvalue table
  for(i in 1:ncol(score)){
    delta <- score[mets_id[keep_index], i] - 
      score[pri_id[keep_index], i]
    group <- HR_HER2[keep_index]
    fit <- lm(delta ~group) 
    sum_fit <- summary(fit)$coefficients

    #brain
    res[i, 1] <- sum_fit[1,1]
    res[i, 2] <- sum_fit[1,2]
    res[i, 3] <- sum_fit[1,4]
    
    # other site
    for(j in 1:num_other_levels){
      K <- matrix(rep(0, num_other_levels+1), nrow=1)
      K[1, 1] <- K[1, j+1] <- 1
      lc_sum <- summary(glht(fit, linfct = K))
      res[i, 4+(j-1)*4+1] <- lc_sum$test$coefficients
      #res[i, 4+(j-3)*2+2] <- round(lincom[[6]][2], digits=2)
      res[i, 4+(j-1)*4+2] <- lc_sum$test$sigma
      res[i, 4+(j-1)*4+3] <- lc_sum$test$pvalue[1]
    }

    # delta
    for(j in 1:num_other_levels){
      res[i, (num_other_levels+1)*4+(j-1)*4+1] <- sum_fit[j+1, 1]
      res[i, (num_other_levels+1)*4+(j-1)*4+2] <- sum_fit[j+1, 2]
      res[i, (num_other_levels+1)*4+(j-1)*4+3] <- sum_fit[j+1, 4]
    }
  }

  colnames(res) <- c(sapply(1:length(uni_HR_HER2_p), function(x)
    paste0(uni_HR_HER2_p[x], c("_est", "_se", "_pval", "_qval"))),
    sapply(1:(length(uni_HR_HER2_p)-1), function(x)
      paste0(uni_HR_HER2_p[x+1], 
      c("_delta_est", "_delta_se", "_delta_pval", "_delta_qval"))))
  rownames(res) <- colnames(score)

  # multiple test correction
  p_val <- res[, seq(3, ncol(res), by=4)]
  p_adj <- matrix(p.adjust(p_val, method="fdr"), 
    ncol=2*(num_other_levels+1)-1, byrow=FALSE)
  res[, seq(4, ncol(res), by=4)] <- p_adj
  if(onlyBrain){
    write.csv(res, file=paste0(fileName, "_onlyBrain.csv"))
  }else{
    write.csv(res, file=paste0(fileName, ".csv"))
  }
  
  # boxplot
  delta <- score[mets_id[keep_index], ] - score[pri_id[keep_index], ]
  delta <- as.data.frame(delta)
  delta$mets_id <- mets_id[keep_index]
  delta$HR_HER2 <- HR_HER2[keep_index]
  #mets_site[keep_index]
  delta <- delta[which(!is.na(delta$HR_HER2)), ]
  data_plot <- melt(delta, 
    variable.name="cellType", value.name="delta_score", 
    id.vars=c("mets_id", "HR_HER2"))
  #data_plot$HR_HER2 <- factor(data_plot$HR_HER2, levels=uni_HR_HER2)

  p <- ggplot(data_plot, aes(x=cellType, y=delta_score)) + 
      geom_boxplot(aes(fill=HR_HER2)) + 
      theme(axis.text.x=element_text(angle=45,hjust=1))+
      #geom_point(aes(fill=cellType), size=5, shape=21, 
      # position=position_jitterdodge()) + 
      #facet_wrap(~cellType, scale="free", nrow=num_row)+
      scale_fill_manual(values=uni_HR_HER2_color)+
      labs(x="", y="")

  if(onlyBrain){
    pdf(file=paste0(fileName, "_onlyBrain.pdf"), width=width, height=height)
    plot(p)
    dev.off() 
  }else{
    pdf(file=paste0(fileName, ".pdf"), width=width, height=height)
    plot(p)
    dev.off() 
  }

}


### HR and HER2
paired_test_HR_HER2_ESTIMATE <- function(score, fileName,
  width, height, onlyBrain=FALSE){
  
  ggcols = gg_color_hue(11)
  uni_HR_HER2_color <- ggcols[c(4,7,2,9,3)]

  if(!onlyBrain){
    uni_HR_HER2_p <- uni_HR_HER2
    }else{
      uni_HR_HER2_p <- uni_HR_HER2_brain
    }
  
  if(!onlyBrain){
    keep_index <- seq(1, length(pri_id))
  }else{
    keep_index <- which(
      sample_pair_info$mets_site=="Brain")
  }

  # get pvalue table
  num_other_levels <- length(uni_HR_HER2_p)-1
  res <- matrix(NA, ncol(score), (num_other_levels+1)*4+
    num_other_levels*4)
  # get pvalue table
  for(i in 1:ncol(score)){
    delta <- score[mets_id[keep_index], i] - 
      score[pri_id[keep_index], i]
    group <- HR_HER2[keep_index]
    fit <- lm(delta ~group) 
    sum_fit <- summary(fit)$coefficients

    #brain
    res[i, 1] <- sum_fit[1,1]
    res[i, 2] <- sum_fit[1,2]
    res[i, 3] <- sum_fit[1,4]
    
    # other site
    for(j in 1:num_other_levels){
      K <- matrix(rep(0, num_other_levels+1), nrow=1)
      K[1, 1] <- K[1, j+1] <- 1
      lc_sum <- summary(glht(fit, linfct = K))
      res[i, 4+(j-1)*4+1] <- lc_sum$test$coefficients
      #res[i, 4+(j-3)*2+2] <- round(lincom[[6]][2], digits=2)
      res[i, 4+(j-1)*4+2] <- lc_sum$test$sigma
      res[i, 4+(j-1)*4+3] <- lc_sum$test$pvalue[1]
    }

    # delta
    for(j in 1:num_other_levels){
      res[i, (num_other_levels+1)*4+(j-1)*4+1] <- sum_fit[j+1, 1]
      res[i, (num_other_levels+1)*4+(j-1)*4+2] <- sum_fit[j+1, 2]
      res[i, (num_other_levels+1)*4+(j-1)*4+3] <- sum_fit[j+1, 4]
    }
  }

  colnames(res) <- c(sapply(1:length(uni_HR_HER2_p), function(x)
    paste0(uni_HR_HER2_p[x], c("_est", "_se", "_pval", "_qval"))),
    sapply(1:(length(uni_HR_HER2_p)-1), function(x)
      paste0(uni_HR_HER2_p[x+1], 
      c("_delta_est", "_delta_se", "_delta_pval", "_delta_qval"))))
  rownames(res) <- colnames(score)

  # multiple test correction
  p_val <- res[, seq(3, ncol(res), by=4)]
  p_adj <- matrix(p.adjust(p_val, method="fdr"), 
    ncol=2*(num_other_levels+1)-1, byrow=FALSE)
  res[, seq(4, ncol(res), by=4)] <- p_adj
  if(onlyBrain){
    write.csv(res, file=paste0(fileName, "_onlyBrain.csv"))
  }else{
    write.csv(res, file=paste0(fileName, ".csv"))
  }
  
  # boxplot
  delta <- score[mets_id[keep_index], ] - score[pri_id[keep_index], ]
  delta <- as.data.frame(delta)
  delta$mets_id <- mets_id[keep_index]
  delta$HR_HER2 <- HR_HER2[keep_index]
  #mets_site[keep_index]
  delta <- delta[which(!is.na(delta$HR_HER2)), ]
  data_plot <- melt(delta, 
    variable.name="cellType", value.name="delta_score", 
    id.vars=c("mets_id", "HR_HER2"))
  #data_plot$HR_HER2 <- factor(data_plot$HR_HER2, levels=uni_HR_HER2)

  p <- ggplot(data_plot, aes(x=HR_HER2, y=delta_score)) + 
      geom_boxplot(aes(fill=HR_HER2)) + 
      theme(axis.text.x=element_text(angle=45,hjust=1))+
      #geom_point(aes(fill=cellType), size=5, shape=21, 
      # position=position_jitterdodge()) + 
      #facet_wrap(~cellType, scale="free", nrow=num_row)+
      scale_fill_manual(values=uni_HR_HER2_color)+
      labs(x="", y="")

  if(onlyBrain){
    pdf(file=paste0(fileName, "_onlyBrain.pdf"), width=width, height=height)
    plot(p)
    dev.off() 
  }else{
    pdf(file=paste0(fileName, ".pdf"), width=width, height=height)
    plot(p)
    dev.off() 
  }

}

