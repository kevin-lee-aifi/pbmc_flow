plot_sample_availability <- function(df, 
                                   timepoint_labels,
                                   subject_col = "subject",
                                   visit_col = "visit",
                                   plot_width = 10,
                                   plot_height = 10,
                                   title = "Sample Availability by Subject and Timepoint",
                                   x_label = "Visit/Timepoint",
                                   y_label = "Subject",
                                   present_color = "steelblue",
                                   absent_color = "lightgray") {
  
  # Rename columns to standard names
  names(df)[names(df) == subject_col] <- "subject"
  names(df)[names(df) == visit_col] <- "visit"
  
  # Set plot dimensions
  options(repr.plot.width = plot_width, repr.plot.height = plot_height)
  
  # Calculate timepoints per subject (for potential future use)
  timepoints_per_subject <- aggregate(visit ~ subject, df, function(x) length(unique(x)))
  names(timepoints_per_subject)[2] <- "num_timepoints"
  
  # Create presence matrix
  timepoint_data <- unique(df[, c("subject", "visit")])
  timepoint_data$present <- 1
  
  all_subjects <- unique(timepoint_data$subject)
  all_visits <- unique(timepoint_data$visit)
  
  complete_grid <- expand.grid(subject = all_subjects, visit = all_visits)
  plot_data <- merge(complete_grid, timepoint_data, all.x = TRUE)
  plot_data$present[is.na(plot_data$present)] <- 0
  
  # Apply factor levels and labels
  plot_data$visit <- factor(plot_data$visit, levels = names(timepoint_labels), labels = timepoint_labels)
  plot_data$subject <- factor(plot_data$subject, levels = sort(all_subjects))
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = visit, y = subject, fill = factor(present))) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c("0" = absent_color, "1" = present_color), 
                      labels = c("0" = "False", "1" = "True"),
                      name = "Data Available") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, x = x_label, y = y_label)
  
  return(p)
}




plot_duplicate_analysis <- function(df,
                                  plot_width = 28,
                                  plot_height = 8,
                                  point_size = 3,
                                  alpha = 0.7,
                                  jitter_width = 0.1,
                                  ncol = 3,
                                  nrow = 1,
                                  percent_col = "percent",
                                  subject_col = "subject",
                                  visit_col = "visit", 
                                  cell_col = "cell") {
  
  # Set plot dimensions
  options(repr.plot.width = plot_width, repr.plot.height = plot_height)
  
  # Find duplicates
  duplicate_counts <- aggregate(df[[percent_col]], 
                              by = list(subject = df[[subject_col]], 
                                       visit = df[[visit_col]], 
                                       cell = df[[cell_col]]), 
                              FUN = length)
  names(duplicate_counts)[4] <- "N"
  
  duplicates_only <- duplicate_counts[duplicate_counts$N > 1, ]
  duplicate_combos <- unique(duplicates_only[, c("subject", "visit")])
  
  duplicate_data <- merge(df, duplicate_combos, by = c(subject_col, visit_col))
  duplicate_data$combo_id <- paste(duplicate_data[[subject_col]], duplicate_data[[visit_col]], sep = " - ")
  
  # Calculate duplicate statistics
  duplicate_stats <- aggregate(duplicate_data[[percent_col]],
                             by = list(cell = duplicate_data[[cell_col]], 
                                     combo = duplicate_data$combo_id),
                             FUN = function(x) c(mean = mean(x), 
                                               sd = sd(x),
                                               min = min(x),
                                               max = max(x),
                                               diff = max(x) - min(x)))
  
  duplicate_stats <- data.frame(
    cell = duplicate_stats$cell,
    combo = duplicate_stats$combo,
    mean_percent = duplicate_stats$x[,"mean"],
    sd_percent = duplicate_stats$x[,"sd"],
    min_percent = duplicate_stats$x[,"min"],
    max_percent = duplicate_stats$x[,"max"],
    diff_percent = duplicate_stats$x[,"diff"]
  )
  
  # Create plots
  p1 <- ggplot(duplicate_stats, aes(x = reorder(cell, sd_percent), y = sd_percent, fill = combo)) +
    geom_col(position = "dodge", alpha = alpha) +
    coord_flip() +
    labs(title = "Standard Deviation of Duplicate Measurements",
         x = "",
         y = "Standard Deviation (%)") +
    theme_minimal()
  
  p2 <- ggplot(duplicate_data, aes(x = !!sym(cell_col), y = !!sym(percent_col), color = combo_id)) +
    geom_point(size = point_size, alpha = alpha, position = position_jitter(width = jitter_width)) +
    coord_flip() +
    labs(title = "Individual Duplicate Measurements",
         x = "",
         y = "Percent of Live Cells") +
    theme_minimal()
  
  p3 <- ggplot(duplicate_stats, aes(x = reorder(cell, diff_percent), y = diff_percent, fill = combo)) +
    geom_col(position = "dodge", alpha = alpha) +
    coord_flip() +
    labs(title = "Difference Between Duplicate Measurements",
         x = "",
         y = "Difference (Max - Min) %") +
    theme_minimal()
  
  plots <- list(p1, p2, p3)
  
  combined_plot <- ggpubr::ggarrange(plotlist = plots, common.legend = TRUE, ncol = ncol, nrow = nrow)
  
  return(list(
    combined_plot = combined_plot,
    individual_plots = plots,
    duplicate_stats = duplicate_stats,
    duplicate_data = duplicate_data
  ))
}



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
                               


process_data <- function(df,
                        subject_col, 
                        visit_col, 
                        population_col, 
                        percent_of_col, 
                        percent_col, 
                        plate_col,
                        celltypes, 
                        viable_filter) {
    timepoint_labels <- c(
        "Flu Year 1 Stand-Alone" = "Y1 SA",
        "Flu Year 1 Day 0" = "Y1 Day 0", 
        "Flu Year 1 Day 7" = "Y1 Day 7",
        "Flu Year 1 Day 90" = "Y1 Day 90",
        "Flu Year 2 Stand-Alone" = "Y2 SA",
        "Flu Year 2 Day 0" = "Y2 Day 0",
        "Flu Year 2 Day 7" = "Y2 Day 7", 
        "Flu Year 2 Day 90" = "Y2 Day 90"
    )
    df <- df[df[[visit_col]] %in% names(timepoint_labels), ]
    df$Sample_ID <- paste(df[[subject_col]], df[[visit_col]], sep = "_")
    df <- df[df[[population_col]] %in% celltypes & df[[percent_of_col]] == viable_filter, ]
    df$total_percents <- ave(df[[percent_col]], df[[subject_col]], df[[visit_col]], 
                            df[[percent_of_col]], FUN = function(x) sum(x, na.rm = TRUE))
    df$frequency <- df[[percent_col]] / df$total_percents
    df <- df[df[[plate_col]] != "B196", ]
    return(df)
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


                                  
testCelltype_paired_clr <- function(clr_data, 
                                      celltype, 
                                      pairsOfTimes,
                                      size = 5,
                                      population_col = "population",
                                      show_pval = "sig") {
  
  # Timepoint mapping
  timepoint_labels <- c(
    "Flu Year 1 Stand-Alone" = "Y1 SA",
    "Flu Year 1 Day 0" = "Y1 Day 0", 
    "Flu Year 1 Day 7" = "Y1 Day 7",
    "Flu Year 1 Day 90" = "Y1 Day 90",
    "Flu Year 2 Stand-Alone" = "Y2 SA",
    "Flu Year 2 Day 0" = "Y2 Day 0",
    "Flu Year 2 Day 7" = "Y2 Day 7", 
    "Flu Year 2 Day 90" = "Y2 Day 90"
  )
  
  # Filter for specific celltype
  cMat <- clr_data[clr_data[[population_col]] == celltype, ]
    
  # Extract visit from Sample_ID (assuming format: Subject_Visit)
  cMat$visit_full <- sub(".*_", "", cMat$Sample_ID)
  cMat$subject_extracted <- sub("_.*", "", cMat$Sample_ID)
  
  # Map full timepoint names to shortened versions
  cMat$visit_extracted <- timepoint_labels[cMat$visit_full]
  
  # Filter for visits in pairsOfTimes
  all_visits <- unique(unlist(pairsOfTimes))
  cMat <- cMat[cMat$visit_extracted %in% all_visits, ]
  
  # Results storage
  pairwise_results <- data.frame(
    Time1 = character(),
    Time2 = character(),
    pval = numeric(),
    p_value = character(),
    p_value_all = character(),
    p_value2 = character(),
    Time1Val = numeric(),
    Time2Val = numeric(),
    Log2FC = numeric(),
    stringsAsFactors = FALSE
  )
  
  subjects_with_pairs <- c()
  
  # Loop through each timepoint pair
  for (i in 1:length(pairsOfTimes)) {
    pair_visits <- pairsOfTimes[[i]]
      
    # Filter data for this pair
    pair_data <- cMat[cMat$visit_extracted %in% pair_visits, ]
    
    # Only include subjects with data at both timepoints
    subject_counts <- aggregate(pair_data$frequency_clr, 
                               by = list(subject = pair_data$subject_extracted), 
                               FUN = length)
    complete_subjects <- subject_counts$subject[subject_counts$x == 2]
    pair_data <- pair_data[pair_data$subject_extracted %in% complete_subjects, ]
    
    subjects_with_pairs <- c(subjects_with_pairs, complete_subjects)
    
    if (nrow(pair_data) > 0) {
      # Create timepoint factor for regression
      pair_data$timepoint <- factor(pair_data$visit_extracted, levels = pair_visits)
      
      # Perform linear regression: CLR ~ timepoint
      if (length(unique(pair_data$timepoint)) == 2 && nrow(pair_data) > 3) {
        model <- lm(frequency_clr ~ timepoint, data = pair_data)
        model_summary <- summary(model)
        
        # Extract p-value for timepoint effect
        timepoint_coeff_name <- paste0("timepoint", pair_visits[2])
        pval <- if(timepoint_coeff_name %in% rownames(model_summary$coefficients)) {
          model_summary$coefficients[timepoint_coeff_name, "Pr(>|t|)"]
        } else { 
          NA 
        }
        
        # Calculate median values for each timepoint
        time1_data <- pair_data[pair_data$timepoint == pair_visits[1], ]
        time2_data <- pair_data[pair_data$timepoint == pair_visits[2], ]
        
        time1_val <- median(time1_data$frequency_clr, na.rm = TRUE)
        time2_val <- median(time2_data$frequency_clr, na.rm = TRUE)
        
      } else {
        pval <- NA
        time1_val <- NA
        time2_val <- NA
      }
    } else {
      pval <- NA
      time1_val <- NA
      time2_val <- NA
    }
    
    # Store results
    pairwise_results <- rbind(pairwise_results, data.frame(
      Time1 = pair_visits[1],
      Time2 = pair_visits[2], 
      pval = pval,
      p_value = ifelse(!is.na(pval) && pval < 0.05, 'p<0.05', ''),
      p_value_all = ifelse(!is.na(pval), paste0('p=', round(pval, 3)), ''),
      p_value2 = ifelse(!is.na(pval) && pval < 0.05, '*', ''),
      Time1Val = time1_val,
      Time2Val = time2_val,
      Log2FC = ifelse(!is.na(time1_val) && !is.na(time2_val) && time1_val > 0, 
                     log2(time2_val / time1_val), NA)
    ))
  }
  
  # Prepare plot data using CLR values
  plot_data <- cMat[cMat$subject_extracted %in% unique(subjects_with_pairs), ]
  plot_data$Value <- plot_data$frequency_clr
  
  # Define proper timepoint order - only use timepoints present in the data
  all_timepoint_order <- c("Y1 SA", "Y1 Day 0", "Y1 Day 7", "Y1 Day 90", 
                          "Y2 SA", "Y2 Day 0", "Y2 Day 7", "Y2 Day 90")
  present_timepoints <- unique(plot_data$visit_extracted)
  timepoint_order <- all_timepoint_order[all_timepoint_order %in% present_timepoints]
  plot_data$visit_extracted <- factor(plot_data$visit_extracted, levels = timepoint_order)
  
  # Create plot
  p <- ggplot(plot_data, aes(x = visit_extracted, y = Value, fill = visit_extracted)) +
    geom_boxplot(color = "black") +
    geom_point(aes(color = subject_extracted), alpha = 0.3) +
    geom_line(aes(group = subject_extracted, color = subject_extracted), alpha = 0.3) +
    theme_by_timepoint() +
    labs(title = celltype, x = "", y = "CLR") +
    ylim(NA, max(plot_data$Value, na.rm = TRUE) + 1.0)
  
  # Add p-value annotations between timepoint pairs
  for (i in 1:nrow(pairwise_results)) {
    # Determine which p-value to show based on show_pval parameter
    pval_to_show <- if (show_pval == "all") {
      pairwise_results$p_value_all[i]
    } else {
      pairwise_results$p_value[i]
    }
    
    if (pval_to_show != "") {  # Only show if there's a p-value to display
      time1 <- pairwise_results$Time1[i]
      time2 <- pairwise_results$Time2[i]
      
      # Get positions of the two timepoints in the ACTUAL plot data
      time1_pos <- which(timepoint_order == time1)
      time2_pos <- which(timepoint_order == time2)
      
      # Only add annotation if both timepoints are present in current data
      if (length(time1_pos) > 0 && length(time2_pos) > 0) {
        # Calculate position between the two timepoints
        x_pos <- (time1_pos + time2_pos) / 2
        
        p <- p + annotate(
          "text",
          x = x_pos,
          y = max(plot_data$Value, na.rm = TRUE) + 0.3,
          label = pval_to_show,
          size = size,
          color = "red"
        )
      }
    }
  }
  
  pairwise_results$Celltype <- celltype
  return(list(Res = pairwise_results, Plot = p))
}


                                  
create_longitudinal_cell_plots_clr <- function(clr_data,
                                               celltypes,
                                               timepoint_labels, 
                                               timepoints_pairs,
                                               population_col = "population",
                                               cell_type_name = "cell",
                                               subplot_width = 4,
                                               subplot_height = 4,
                                               ncol = NULL,
                                               stat_test_size = 4,
                                               show_pval = "sig",
                                               title_size = 20,
                                               axis_label_size = 20,
                                               year_titles = c("Year 1", "Year 2"),
                                               x_axis_label = "Clinical Visits",
                                               y_axis_label = "Centered Log Ratio",
                                               plot_path = NULL,
                                               filename_prefix = "longitudinal_clr",
                                               save_width = 20,
                                               save_height = 12,
                                               dpi = 300) {
  
  # Calculate grid dimensions based on number of cell types
  num_cells <- length(celltypes)
  
  # Auto-determine ncol if not specified
  if (is.null(ncol)) {
    ncol <- min(7, max(3, ceiling(sqrt(num_cells))))
  }
  
  # Calculate number of rows needed and total plot dimensions
  nrow <- ceiling(num_cells / ncol)
  plot_width <- ncol * subplot_width
  plot_height <- nrow * subplot_height
  
  # Set plot dimensions
  options(repr.plot.width = plot_width, repr.plot.height = plot_height)
  
  # Function to create multi-panel plot with one subplot per cell type
  create_figure <- function(pairs, title) {
    plotlist <- lapply(celltypes, function(x) {
      result <- testCelltype_paired_clr(clr_data = clr_data, 
                                       celltype = x, 
                                       pairsOfTimes = pairs, 
                                       size = stat_test_size,
                                       population_col = population_col,
                                       show_pval = show_pval)
      return(result$Plot)
    })
    # Remove NULL plots
    plotlist <- plotlist[!sapply(plotlist, is.null)]
    
    if (length(plotlist) == 0) {
      return(NULL)
    }
    
    fig <- ggpubr::ggarrange(plotlist = plotlist, common.legend = TRUE, ncol = ncol, nrow = nrow)
    fig <- annotate_figure(fig, top = text_grob(title, size = title_size))
    annotate_figure(fig, bottom = text_grob(x_axis_label, size = axis_label_size),
                   left = text_grob(y_axis_label, size = axis_label_size, rot = 90))
  }
  
  # Generate figures for each year
  plots <- list()
  year_names <- names(timepoints_pairs)
  
  for (i in seq_along(timepoints_pairs)) {
    year_name <- year_names[i]
    full_title <- paste(year_titles[i], "Relative", cell_type_name, "abundances in PBMCs")
    plots[[year_name]] <- create_figure(timepoints_pairs[[i]], full_title)
    
    # Save plot if plot_path is provided
    if (!is.null(plot_path) && !is.null(plots[[year_name]])) {
      filename <- paste0(gsub(" ", "_", year_name), "_", filename_prefix, ".png")
      ggsave(file.path(plot_path, filename), 
             plot = plots[[year_name]], 
             width = save_width, 
             height = save_height, 
             dpi = dpi,
             bg = "white")  # Set explicit white background
    }
  }
  
  return(plots)
}