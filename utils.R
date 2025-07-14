suppressPackageStartupMessages({
    library(ggplot2)
    library(ggalluvial)
    library(dplyr)
    library(RColorBrewer)
    library(gridExtra)
    library(rlang)
    library(tidyr)
    library(patchwork)
    library(data.table)
    library(ggpubr)
    library(tidyverse)
    library(compositions)
})



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
    scale_fill_manual(values = tp_color_map, name = "visit"),
    scale_color_manual(values = color_map, name = "subject"),
    guides(
      fill = guide_legend(order = 1),
      color = guide_legend(order = 2)
    )
  )
}



testCelltype_paired <- function(Mat, celltype, pairsOfTimes, size=5, use_duplicates = FALSE) {
  if (!is.data.table(Mat)) {
    Mat <- as.data.table(Mat)
  }
  
  cMat = Mat[cell==celltype & visit %in% unique(unlist(pairsOfTimes))]
  
  if (use_duplicates) {
    cMat <- cMat[, .(cell_type_frac_total = mean(cell_type_frac_total)), 
                 by = .(subject, visit, cell)]
  } else {
    duplicate_counts <- cMat[, .N, by = .(subject, visit, cell)]
    duplicate_combos <- unique(duplicate_counts[N > 1, .(subject, visit)])
    if (nrow(duplicate_combos) > 0) {
      cMat <- cMat[!paste(subject, visit) %in% paste(duplicate_combos$subject, duplicate_combos$visit)]
    }
  }
  
  cMat$Value = cMat$cell_type_frac_total
  
  time_points <- levels(cMat$visit)
  
  pairsOfTimes = pairsOfTimes
  
  subjects_with_pairs <- c()
  
  pairwise_results <- data.frame(
    Time1 = character(),
    Time2 = character(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  for (i in 1:length(pairsOfTimes)) {
    pairwise_data <- cMat[visit %in% pairsOfTimes[[i]]]
    pairwise_data = dcast(pairwise_data, subject ~ visit, value.var = 'Value')
    pairwise_data = pairwise_data[!is.na(pairwise_data[[2]]) & !is.na(pairwise_data[[3]])]   
    
    subjects_with_pairs <- c(subjects_with_pairs, pairwise_data$subject)
    
    test_result <- wilcox.test(
      pairwise_data[[2]], 
      pairwise_data[[3]],
      paired = TRUE,
      exact = FALSE
    )
    
    medDiff = median(pairwise_data[[2]] - pairwise_data[[3]])
    pval = test_result$p.value
    pairwise_results <- rbind(pairwise_results, data.frame(
      Time1 = pairsOfTimes[[i]][1],
      Time2 = pairsOfTimes[[i]][2],
      pval = pval,
      p_value = ifelse(pval < 0.05,'p<0.05', ''),
      p_value2 = ifelse(pval < 0.05,'*',''),
      Time1Val = median(pairwise_data[[2]]),
      Time2Val = median(pairwise_data[[3]]),
      Log2FC = log2(median(pairwise_data[[3]]) / median(pairwise_data[[2]]))
    ))
  }
  
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



freq_clr <- function(freq_table, sample_col, freq_col, celltype_col){
    freq_selected <- select(freq_table, all_of(c(freq_col, celltype_col, sample_col)))
    freq_grouped <- group_by(freq_selected, across(all_of(c(sample_col, celltype_col))))
    freq_clean <- summarise(freq_grouped, across(all_of(freq_col), function(x) mean(x, na.rm = TRUE)), .groups = 'drop')
    
    freq <- pivot_wider(freq_clean, 
                       id_cols = all_of(sample_col), 
                       names_from = all_of(celltype_col),
                       values_from = all_of(freq_col))
    
    freq_mx_cols <- select(freq, -all_of(sample_col))
    freq_mx <- as.matrix(freq_mx_cols)
    rownames(freq_mx) <- freq[[sample_col]]
    
    freq_clr_matrix <- compositions::clr(freq_mx)
    freq_clr_tibble <- as_tibble(freq_clr_matrix, rownames = sample_col)
    freq_clr <- pivot_longer(freq_clr_tibble,
                           cols = -all_of(sample_col), 
                           names_to = celltype_col, 
                           values_to = paste0(freq_col, '_clr'))
    
    freq_table_grouped <- group_by(freq_table, across(all_of(c(sample_col, celltype_col))))
    freq_table_clean <- summarise(freq_table_grouped,
                                 across(where(is.numeric), function(x) mean(x, na.rm = TRUE)),
                                 across(where(is.character), first),
                                 .groups = 'drop')
    
    freq_meta_clr <- full_join(freq_table_clean, freq_clr, by = c(sample_col, celltype_col))
    
    return(freq_meta_clr)
}