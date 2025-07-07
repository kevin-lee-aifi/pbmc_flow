theme_by_timepoint <- function(base_size = 14,
                                   tp_color_map = c(
                                     "Y1 SA" = "#babaa6",
                                     "Y1 Day 0" = "#bfa600", 
                                     "Y1 Day 7" = "#a94b23",
                                     "Y1 Day 90" = "#a0204f",
                                     "Y2 SA" = "#f4f4c3",
                                     "Y2 Day 0" = "#ffec6e",
                                     "Y2 Day 7" = "#f6a27e",
                                     "Y2 Day 90" = "#f679a7"
                                   ),
                                   color_map = c(
                                     "FH1002" = "#ff7f0e", "FH1003" = "#279e68",
                                     "FH1004" = "#d62728", "FH1005" = "#aa40fc", "FH1006" = "#8c564b",
                                     "FH1007" = "#e377c2", "FH1008" = "#b5bd61", "FH1009" = "#17becf",
                                     "FH1011" = "#ffbb78", "FH1012" = "#98df8a",
                                     "FH1014" = "#ff9896", "FH1016" = "#c5b0d5", "FH1017" = "#c49c94"
                                   )) {
  list(
    theme_bw(base_size = base_size) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 4),
        axis.title = element_text(size = base_size, face = "bold"),
        axis.text = element_text(size = base_size - 2),
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none"
      ),
    scale_fill_manual(values = tp_color_map),
    scale_color_manual(values = color_map)
  )
}

testCelltype_paired <- function(Mat, celltype, pairsOfTimes, size=5) {
  cMat = Mat[cell==celltype & visit %in% unique(unlist(pairsOfTimes))]
  
  cMat$Value = cMat$cell_type_frac_total
  
  # Get all consecutive time point pairs
  time_points <- levels(cMat$visit)
  
  pairsOfTimes = pairsOfTimes
  
  # Identify subjects with paired data for plotting
  subjects_with_pairs <- c()
  
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
    
    subjects_with_pairs <- c(subjects_with_pairs, pairwise_data$subject)
    
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
  }
  
  # Filter cMat to only subjects with paired data
  cMat_plot <- cMat[subject %in% unique(subjects_with_pairs)]
  
    p <- ggplot(cMat_plot, aes(x = visit, y = Value, fill = visit)) +
      geom_boxplot(color = "black") +
      geom_point(aes(color = subject), alpha = 0.3) +
      geom_line(aes(group = subject, color = subject), alpha = 0.3) +
      theme_by_timepoint() +
      labs(title = celltype, x = "", y = "") +
      ylim(NA, max(cMat_plot$Value) + 0.7)
  
  for (i in 1:nrow(pairwise_results)) {
    p <- p + annotate(
      "text",
      x = i + 0.5,
      y = max(cMat_plot$Value) + 0.5,
      label = pairwise_results$p_value[i],
      size =size ,
      color = "red"
    )
  }
  pairwise_results$Celltype=celltype
  return(list(Res=pairwise_results,
              Plot =p))
}