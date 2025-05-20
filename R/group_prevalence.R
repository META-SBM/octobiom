#' Generate Prevalence Comparison Plots for Time Points
#'
#' @param physeq A phyloseq object
#' @param output_dir directory where figures will be
#' @param time_column Name of the column in sample_data containing time points
#' @param top_n Number of top taxa to label in zoom plots (default: 10)
#' @param point_color Color for scatter plot points (default: "#d5b3ff")
#' @param smooth_color Color for LOESS smooth line (default: "#d62828")
#' @param text_size Size for text labels (default: 3)
#' @param highlight_top Create highlight for top 10 taxa
#' @param width figure's width
#' @param height figure's height
#' @param dpi figure's quality
#' @return A list of ggplot objects for each time point comparison
#' @export
compare_prevalence_across_group <- function(physeq, column, top_n = 10,
                                           point_color = "#d5b3ff",
                                           smooth_color = "#d62828",
                                           text_size = 3,
                                           highlight_top = FALSE) {
  # Input validation
  if (!inherits(physeq, "phyloseq")) {
    stop("❌ Input must be a phyloseq object")
  }

  if (!column %in% colnames(sample_data(physeq))) {
    stop("❌ Column not found in sample data")
  }
  # Get unique time points
  time_points <- unique(as.character(sample_data(physeq)[[column]]))

  if (length(time_points) != 2) {
    stop("❌ Need to  2 points for comparison")
  }
  # Calculate prevalence for each time point
  prevalence_list <- list()
  for (tp in time_points) {
    # Subset phyloseq for this time point
    metadata <- as(physeq@sam_data, "data.frame")
    samples_tp <- rownames(metadata[metadata[[Column]] == tp, ])
    # Subsample phyloseq
    ps_sub <- prune_samples(samples_tp, physeq)
    # Calculate prevalence
    otu_mat <- as(otu_table(ps_sub), "matrix")
    if (!taxa_are_rows(ps_sub)) {
      otu_mat <- t(otu_mat)
    }
    prev <- rowSums(otu_mat > 0) / nsamples(ps_sub)
    # Save results
    prevalence_list[[tp]] <- data.frame(
      Taxa = names(prev),
      Prevalence = prev,
      stringsAsFactors = FALSE
    )
    # Check information
    message("Time point: ", tp)
    message("Samples: ", length(samples_tp))
    message("Taxa: ", ntaxa(ps_sub))
}
  # Create all possible pairwise combinations
  time_combinations <- combn(time_points, 2, simplify = FALSE)
  print(time_combinations)
  # Store plots
  plot_list <- list()
  for (pair in time_combinations) {
    time1 <- pair[1]
    time2 <- pair[2]
    # Merge prevalence data
    merged <- merge(prevalence_list[[time1]],
                    prevalence_list[[time2]],
                    by = "Taxa",
                    suffixes = paste0("_", c(time1, time2)),
                    all = TRUE) %>%
      replace(is.na(.), 0)
    # Get taxonomy information if available
    if (!is.null(tax_table(physeq, errorIfNULL = FALSE))) {
      tax_info <- as.data.frame(tax_table(physeq))[merged$Taxa,
                                                   "Species", drop = FALSE]
      merged$Species <- tax_info$Species
    } else {
      merged$Species <- merged$Taxa
    }
    # Identify top taxa by average prevalence
    top_taxa <- merged %>%
      mutate(avg_prevalence = (.[[2]] + .[[3]])/2) %>%
      arrange(desc(avg_prevalence)) %>%
      head(top_n)
    # Main scatter plot
    main_plot <- ggplot(merged, aes_string(x = paste0("Prevalence_", time1),
                                           y = paste0("Prevalence_", time2))) +
      geom_point(alpha = 0.6, color = point_color, size = 3) +
      geom_smooth(method = "loess", color = smooth_color,
                  se = TRUE, fill = "indianred", alpha = 0.2) +
      scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(
        x = paste("Prevalence in", time1),
        y = paste("Prevalence in", time2),
        title = paste("Prevalence Comparison:", time1, "vs", time2)
      ) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.4),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
      )
    # Zoom plot for top taxa
    zoom_plot <- ggplot(top_taxa, aes_string(x = paste0("Prevalence_", time1),
                                             y = paste0("Prevalence_", time2))) + # nolint
      geom_point(color = point_color, size = 4, alpha = 0.8) +
      geom_text_repel(
        aes(label = sprintf("%s\n(%.1f%%, %.1f%%)", Species,
                            get(paste0("Prevalence_", time1)) * 100,
                            get(paste0("Prevalence_", time2)) * 100)),
        size = text_size,
        box.padding = 0.7,
        point.padding = 0.3,
        segment.size = 0.3,
        min.segment.length = 0,
        bg.color = "white",
        bg.r = 0.15,
        max.iter = 10000
      ) +
      labs(
        x = paste("Prevalence in", time1, "(zoom)"),
        y = paste("Prevalence in", time2, "(zoom)")
      ) +
      scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.4),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
      )
    if (highlight_top == TRUE) {
      top_taxa <- merged %>%
        mutate(avg_prevalence = (.[[2]] + .[[3]])/2) %>%
        arrange(desc(avg_prevalence)) %>%
        head(10)

      zoom_plot <- zoom_plot +
        coord_cartesian(clip = "off") + # To prevent label clipping
        ggforce::geom_mark_hull(
          aes(filter = Species %in% top_taxa$Species),
          description = "Top 10 with the most average prevalence",
          con.cap = 3,
          # Line aesthetics
          con.colour = "gray50",  # Gray color
          con.linetype = "dashed",  # Dashed line
          con.size = 0.8,  # Line thickness
          # Description positioning
          description.vjust = -3,  # Move description up slightly
          description.hjust = 6,   # Move description right slightly
          description.colour = "black",  # Description text color
          linetype = 2, color = "gray50"
        )
    }
    # Combine plots
    combined_plot <- (main_plot / zoom_plot) +
      plot_layout(heights  = c(1, 1))
    plot_list[[paste(time1, "vs", time2)]] <- combined_plot
  }
  return(plot_list)
}
