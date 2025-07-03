testCelltype_paired <- function(Mat, celltype, pairsOfTimes=list(c('Y1 SA','Y1 Day 0'),
                                                                 c('Y1 Day 0','Y1 Day 7'),
                                                                 c('Y1 Day 7', 'Y1 Day 90')),
                                size=5){
  cMat = Mat[Cell==celltype & visit %in% unique(unlist(pairsOfTimes))]
  
  cMat$Value = cMat$cell_type_frac_total
  
  # Get all consecutive time point pairs
  time_points <- levels(cMat$visit)
  
  pairsOfTimes = pairsOfTimes
  
  # Loop through consecutive time points
  pairwise_results <- data.frame(
    Time1 = character(),
    Time2 = character(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  for (i in 1:length(pairsOfTimes)) {
    # Subset data for the two consecutive time points
    pairwise_data <- subset(cMat, visit %in% pairsOfTimes[[i]])
    pairwise_data = dcast(pairwise_data, subject ~ visit, value.var = 'Value')
    pairwise_data = pairwise_data[ !is.na(pairwise_data[,2][[1]]) & 
                                     !is.na(pairwise_data[,3][[1]])]   
    
    # Perform a paired t-test
    test_result <- wilcox.test(
      pairwise_data[,2][[1]], 
      pairwise_data[,3][[1]],
      paired = TRUE
    )
    
    medDiff = median( pairwise_data[,2][[1]]-pairwise_data[,3][[1]])
    # Store the results
    pval =test_result$p.value
    pairwise_results <- rbind(pairwise_results, data.frame(
      Time1 = pairsOfTimes[[i]][1],
      Time2 = pairsOfTimes[[i]][2],
      pval = pval,
      p_value = ifelse(pval < 0.05,'p<0.05', ''),
      p_value2 = ifelse(pval < 0.05,'*',''),
      Time1Val = median(pairwise_data[,2][[1]]),
      Time2Val = median(pairwise_data[,3][[1]]),
      Log2FC = log2(median(pairwise_data[,3][[1]]) /median(pairwise_data[,2][[1]]))
    ))
    #pairwise_results =data.table(pairwise_results)
  }
  p <- ggplot(cMat, aes(x = visit, y = Value)) +
    geom_boxplot(fill = "skyblue", color = "black") +geom_point(alpha=0.3)+
    geom_line(data=cMat,
              aes(x=visit,
                  y=Value,
                  group=subject),
              alpha=0.3)+
    labs(title = celltype,
         x = "",
         y = "") +
    theme_minimal() +ylim(NA, max(cMat$Value) + 0.7)
    theme(text=element_text(size=14),
          axis.text.x = element_text(angle=90))
    #scale_y_log10()
  
  for (i in 1:nrow(pairwise_results)) {
    p <- p + annotate(
      "text",
      x = i + 0.5, # Position between the two time points
      y = max(cMat$Value) + 0.5, # Position above the boxplots
      label = pairwise_results$p_value[i],
      size =size ,
      color = "red"
    )
  }
  pairwise_results$Celltype=celltype
  return(list(Res=pairwise_results,
              Plot =p))
  
}