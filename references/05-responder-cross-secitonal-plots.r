# Suppress package startup messages for cleaner output
suppressPackageStartupMessages({
  library(tidyr)
  library(tibble)
  library(stringr)
  library(parallel)
  library(DESeq2)
  library(ggplot2)
  library(data.table)
  library(qvalue)
  library(ggrepel)
  library(reshape2)
  library(grid)
})

#' Volcano Plot by Cell Type
#'
#' Generates a volcano plot for DESeq2 results by cell type.
#'
#' @param res_file Path to the DESeq2 results CSV file.
#' @param file_save_name Name for the output PDF (without extension).
#' @export
#' @details The input CSV file (`input_deg_res`) must contain the following columns:
#' - `Qvalue`: Numeric values representing the adjusted p-values.
#' - `log2FoldChange`: Numeric values representing the log2 fold change.
#' - `gene`: Character values representing gene names.
#' - `pvalue`: Numeric values representing the raw p-values.
plot_volcano_deg <- function(input_deg_res,
                             title,
                             subtitle,
                             file_save_name = "volcano_plot") {
  if (!file.exists(input_deg_res)) {
    stop("Error: The file does not exist. Please check the file path.")
  }
  res <- tryCatch(
    data.table::fread(input_deg_res),
    error = function(e) {
      stop("Error: Unable to read the file. Please ensure it is a valid CSV file.")
    }
  )
  # Remove rows with NA values in Qvalue or log2FoldChange
  res <- res[!is.na(Qvalue) & !is.na(log2FoldChange), ]
  res$facet_label <- stringr::str_to_title(gsub("_", " ", res$celltype))

  height <- ceiling(length(unique(res$celltype)) / 4) * 7
  png(paste0("files/", file_save_name, "_volcano.png"), width = 20, height = height, units = "in", res = 300)
  print(
    ggplot2::ggplot(
      res,
      ggplot2::aes(
        x = log2FoldChange,
        y = -log10(pvalue),
        color = ifelse(Qvalue < 0.1 & log2FoldChange > 0.5, "Upregulated in Responders",
          ifelse(Qvalue < 0.1 & log2FoldChange < -0.5, "Upregulated in NonResponders",
            "Not Significant"
          )
        ),
        label = gene
      )
    ) +
      ggplot2::geom_point(size = 1.5) +
      ggrepel::geom_text_repel(
        data = subset(res, Qvalue < 0.1),
        force = 1,
        max.overlaps = 20, 
        show.legend = FALSE
      ) +
      ggplot2::geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50", linewidth = 0.5) +
      ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.5) +
      ggplot2::facet_wrap(~facet_label, ncol = 4, scales = "free") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = title,
        subtitle = subtitle,
        x = expression(paste("Log"[2], " Fold Change")),
        y = expression(paste("-Log"[10], " (p-value)"))
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 18, margin = ggplot2::margin(b = 10)),
        axis.title = ggplot2::element_text(size = 16, face = "bold"),
        axis.text = ggplot2::element_text(size = 12),
        legend.position = "bottom",
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 12),
        strip.text = ggplot2::element_text(size = 14, hjust = 0.05),
        panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.2),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = "white", color = NA),
        panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.8),
        panel.spacing = grid::unit(0.5, "cm")
      ) +
      ggplot2::scale_color_manual(
        name = NULL,
        values = c(
          "Upregulated in Responders" = "#00b2a8",
          "Upregulated in NonResponders" = "#d81463",
          "Not Significant" = "gray70"
        )
      ) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5)))
  )
  dev.off()
}

#' Plot Cell Type CLR by Status
#'
#' This function reads a CSV file, filters for pre-treatment and non-dara samples,
#' assigns responder status, and plots cell type CLR as boxplots with jittered points.
#'
#' @param input_csv Path to input CSV file.
#' @param file_save_name Output file name (without extension).
#' @param responder_guids Vector of responder subject GUIDs.
#' @param NonResponder_guids Vector of NonResponder subject GUIDs.
#' @param color_map Named vector mapping subject GUIDs to colors.
plot_celltype_clr <- function(
    input_csv,
    treatment_label,
    celltype_level,
    title,
    subtitle,
    file_save_name = "boxplot",
    color_map = c(
      "FH1001" = "#1f77b4", "FH1002" = "#ff7f0e", "FH1003" = "#279e68",
      "FH1004" = "#d62728", "FH1005" = "#aa40fc", "FH1006" = "#8c564b",
      "FH1007" = "#e377c2", "FH1008" = "#b5bd61", "FH1009" = "#17becf",
      "FH1010" = "#aec7e8", "FH1011" = "#ffbb78", "FH1012" = "#98df8a",
      "FH1014" = "#ff9896", "FH1016" = "#c5b0d5", "FH1017" = "#c49c94",
      "FH1018" = "#f7b6d2", "FH1021" = "#dbdb8d"
    )) {
  # Read input data
  frequency_plots <- fread(input_csv)

  # Filter for pre-treatment and non-dara samples
  frequency_plots <- frequency_plots[sample.visitDetails == treatment_label]
  frequency_plots <- frequency_plots[manual.treatment == "non-dara"]

  # Assign responder status
  frequency_plots$Status <- ""
  responder_guids <- c("FH1002", "FH1005", "FH1006", "FH1008", "FH1012", "FH1014", "FH1017")
  NonResponder_guids <- c("FH1021", "FH1016", "FH1011", "FH1009", "FH1007", "FH1004", "FH1003")
  frequency_plots[subject.subjectGuid %in% NonResponder_guids, ]$Status <- "NonResponder"
  frequency_plots[subject.subjectGuid %in% responder_guids, ]$Status <- "Responder"

  # Keep only rows with assigned status
  frequency_plots <- frequency_plots[Status != ""]
  frequency_plots$facet_label <- stringr::str_to_title(gsub("_", " ", frequency_plots[[celltype_level]]))

  height <- ceiling(length(unique(frequency_plots[[celltype_level]])) / 4) * 7
  # Open PNG device for saving plot
  png(paste0("files/", file_save_name, "_boxfreq.png"), width = 24, height = height, units = "in", res = 300)

  # Create and print ggplot
  print(
    ggplot(frequency_plots, aes(x = Status, y = cell_type_clr, fill = Status)) +
      geom_boxplot(alpha = 0.8, outlier.shape = NA, color = "black") + # Boxplot
      geom_jitter(aes(color = subject.subjectGuid), width = 0.2, size = 2.2, stroke = 1, alpha = 0.8) + # Jittered points
      facet_wrap(vars(facet_label), scales = "free_y", ncol = 4) +
      labs(
        title = title,
        subtitle = subtitle,
        x = "Cell Types per Subject split by Flu Response",
        y = "Centered Log Ratio"
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
      scale_fill_manual(values = c("Responder" = "#00b2a8", "NonResponder" = "#d81463")) +
      scale_color_manual(values = color_map)
  )
  # Close PNG device
  dev.off()
}
