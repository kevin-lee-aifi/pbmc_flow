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



create_celltype_plots <- function(input, 
                                 plot_path, 
                                 celltype_level = "population", 
                                 timepoint_labels = NULL, 
                                 response_var = "Response",
                                 visit_var = "Flu.Visit",
                                 plot_title_prefix = "Celltype Analysis", 
                                 filename_prefix = "celltype_composition",
                                 x_label = "Response Groups", 
                                 y_label = "Centered Log Ratio",
                                 width = 20, 
                                 height = 12, 
                                 dpi = 300) {
  
  # Set plot options
  options(repr.plot.width = 26, repr.plot.height = 16)
  
  # Default color map
  color_map <- c(
    "FH1001" = "#1f77b4", "FH1002" = "#ff7f0e", "FH1003" = "#279e68",
    "FH1004" = "#d62728", "FH1005" = "#aa40fc", "FH1006" = "#8c564b",
    "FH1007" = "#e377c2", "FH1008" = "#b5bd61", "FH1009" = "#17becf",
    "FH1010" = "#aec7e8", "FH1011" = "#ffbb78", "FH1012" = "#98df8a",
    "FH1014" = "#ff9896", "FH1016" = "#c5b0d5", "FH1017" = "#c49c94",
    "FH1018" = "#f7b6d2", "FH1021" = "#dbdb8d"
  )
  
  # Default timepoint labels if not provided
  if (is.null(timepoint_labels)) {
    # data_sample <- fread(input_csv, nrows = 100)
    data_sample <- input
    unique_visits <- unique(data_sample[[visit_var]])
    timepoint_labels <- setNames(unique_visits, unique_visits)
  }
  
  all_results <- list()
  
  # Loop through each timepoint
  for (visit in names(timepoint_labels)) {
    
    # Read and filter data
    # frequency_plots <- fread(input_csv)
    frequency_plots <- input
    frequency_plots <- frequency_plots[frequency_plots[[visit_var]] == visit, ]
    frequency_plots <- frequency_plots[frequency_plots[[response_var]] != "", ]
    frequency_plots$facet_label <- factor(gsub("_", " ", frequency_plots[[celltype_level]]), 
                                     levels = gsub("_", " ", levels(input[[celltype_level]])))
    
    # Perform linear regression for each cell type
    cell_types <- unique(frequency_plots[[celltype_level]])
    regression_results <- data.frame()
    
    for (celltype in cell_types) {
      subset_data <- frequency_plots[frequency_plots[[celltype_level]] == celltype, ]
      
      if (nrow(subset_data) > 3) {
        formula_str <- paste("frequency_clr ~ Sex +", response_var)
        model <- lm(as.formula(formula_str), data = subset_data)
        model_summary <- summary(model)
        
        response_coeff_name <- paste0(response_var, "Responder")
        response_pval <- if(response_coeff_name %in% rownames(model_summary$coefficients)) {
          model_summary$coefficients[response_coeff_name, "Pr(>|t|)"]
        } else { NA }
        
        sex_pval <- if("SexM" %in% rownames(model_summary$coefficients)) {
          model_summary$coefficients["SexM", "Pr(>|t|)"]
        } else { NA }
        
        response_coeff <- if(response_coeff_name %in% rownames(model_summary$coefficients)) {
          model_summary$coefficients[response_coeff_name, "Estimate"]
        } else { NA }
        
        regression_results <- rbind(regression_results, data.frame(
          celltype = celltype,
          facet_label = gsub("_", " ", celltype),
          response_pval = response_pval,
          sex_pval = sex_pval,
          response_coeff = response_coeff,
          adj_r_squared = model_summary$adj.r.squared
        ))
      }
    }
    
    # Format p-values
    regression_results$response_p_text <- ifelse(
      is.na(regression_results$response_pval) | regression_results$response_pval > 0.05, "",
      paste("p=", round(regression_results$response_pval, 3))
    )
    regression_results$sex_p_text <- ifelse(
      is.na(regression_results$sex_pval) | regression_results$sex_pval > 0.05, "",
      paste("Sex: p=", round(regression_results$sex_pval, 3))
    )

    max_y <- aggregate(frequency_clr ~ facet_label, data = frequency_plots, max)
    max_y$facet_label <- factor(max_y$facet_label, levels = levels(frequency_plots$facet_label))
    regression_results$facet_label <- factor(regression_results$facet_label, levels = levels(frequency_plots$facet_label))
    regression_results <- merge(regression_results, max_y, by = "facet_label")
    regression_results$facet_label <- factor(regression_results$facet_label, levels = levels(frequency_plots$facet_label))
    regression_results$y_response <- regression_results$frequency_clr + 0.3
    regression_results$y_sex <- regression_results$frequency_clr + 0.6
        
    # Create plot
    p <- ggplot(frequency_plots, aes(x = .data[[response_var]], y = frequency_clr, fill = .data[[response_var]])) +
      geom_boxplot(alpha = 0.8, outlier.shape = NA, color = "black") +
      geom_jitter(aes(color = Subject), width = 0.2, size = 2.2, stroke = 1, alpha = 0.8) +
      geom_text(data = regression_results, 
                aes(x = 2, y = y_response, label = response_p_text),
                inherit.aes = FALSE, size = 3, color = "red", fontface = "bold") +
      geom_text(data = regression_results, 
                aes(x = 2, y = y_sex, label = sex_p_text),
                inherit.aes = FALSE, size = 2.5, color = "blue") +
      facet_wrap(vars(facet_label), scales = "free_y", ncol = 4) +
      labs(
        title = paste0(visit, ": ", plot_title_prefix),
        subtitle = visit,
        x = x_label,
        y = y_label
      ) +
      theme_bw(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 10)),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.position = "right",
        legend.title = element_blank(),
        strip.text = element_text(size = 14, face = "bold")
      ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
      scale_fill_manual(values = c("Responder" = "#00b2a8", "Non responder" = "#d81463")) +
      scale_color_manual(values = color_map)
    
    # Save plot
    filename <- paste0(gsub(" ", "_", visit), "_", filename_prefix, ".png")
    ggsave(file.path(plot_path, filename), plot = p, width = width, height = height, dpi = dpi)
    
    all_results[[visit]] <- p
  }
  
  return(all_results)
}