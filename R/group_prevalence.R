#' Generate Prevalence Comparison Plots for All Pairwise Group Combinations
#'
#' @param physeq A phyloseq object
#' @param column Name of the column in sample_data containing groups to compare
#' @param top_n Number of top taxa to label in zoom plots (default: 10)
#' @param point_color Color for scatter plot points (default: "#d5b3ff")
#' @param smooth_color Color for LOESS smooth line (default: "#d62828")
#' @param text_size Size for text labels (default: 3)
#' @param zoom_window list of x and y range (example : list(x = c(0.1, 0.5), y = c(0.1, 0.5)))
#' @param zoom_plot_flag add zoom for top taxa
#' @return A list of ggplot objects for each pairwise group comparison
#' @export
compare_prevalence_across_group <- function(physeq, column,
                                            point_color = "#d5b3ff",
                                            smooth_color = "#d62828",
                                            text_size = 3,
                                            zoom_plot_flag = T,
                                            zoom_window = list(x = c(0.1, 0.5), y = c(0.1, 0.5))) {
  # Input validation
  if (!inherits(physeq, "phyloseq")) {
    stop("❌ Input must be a phyloseq object")
  }

  if (!column %in% colnames(sample_data(physeq))) {
    stop(paste("❌ Column", column, "not found in sample data"))
  }

  # Get unique groups
  groups <- unique(as.character(sample_data(physeq)[[column]]))

  if (length(groups) < 2) {
    stop("❌ Need at least 2 groups for comparison")
  }

  # Calculate prevalence for each group
  prevalence_list <- list()
  for (group in groups) {
    # Subset phyloseq for this group
    metadata <- as(physeq@sam_data, "data.frame")
    samples_group <- rownames(metadata[metadata[[column]] == group, ])
    ps_sub <- prune_samples(samples_group, physeq)

    # Calculate prevalence
    otu_mat <- as(otu_table(ps_sub), "matrix")
    if (!phyloseq::taxa_are_rows(ps_sub)) {
      otu_mat <- t(otu_mat)
    }
    prev <- rowSums(otu_mat > 0) / phyloseq::nsamples(ps_sub)

    # Save results
    prevalence_list[[group]] <- data.frame(
      Taxa = names(prev),
      Prevalence = prev,
      stringsAsFactors = FALSE
    )

    # Print diagnostic information
    message("Group: ", group)
    message("  Samples: ", length(samples_group))
    message("  Taxa: ", phyloseq::ntaxa(ps_sub))
  }

  # Create all possible pairwise combinations
  group_combinations <- combn(groups, 2, simplify = FALSE)

  # Store plots
  plot_list <- list()
  for (pair in group_combinations) {
    group1 <- pair[1]
    group2 <- pair[2]

    # Merge prevalence data
    merged <- merge(prevalence_list[[group1]],
                    prevalence_list[[group2]],
                    by = "Taxa",
                    suffixes = paste0("_", c(group1, group2)),
                    all = TRUE) %>%
      replace(is.na(.), 0)

    # Get taxonomy information if available
    if (!is.null(tax_table(physeq, errorIfNULL = FALSE))) {
      tax_info <- as.data.frame(tax_table(physeq))[merged$Taxa, , drop = FALSE]
      # Use the lowest available taxonomic level for labeling
      merged$Label <- apply(tax_info, 1, function(x) {
        x <- x[!is.na(x) & x != "" & !is.na(x)]
        if (length(x) == 0) return(merged$Taxa[1])
        x[length(x)]
      })
    } else {
      merged$Label <- merged$Taxa
    }
    x_range <- zoom_window$x
    y_range <- zoom_window$y
    # Identify top taxa by average prevalence
    top_taxa <- merged %>%
      mutate(avg_prevalence = (.[[2]] + .[[3]])/2) %>%
      arrange(desc(avg_prevalence)) #%>%

    top_taxa <- top_taxa %>%
      filter(.[[2]]>x_range[1] & .[[2]] < x_range[2]) %>%
      filter(.[[3]]>y_range[1] & .[[3]] < y_range[2])
      #head(top_n)

    # Main scatter plot
    main_plot <- ggplot2::ggplot(merged, ggplot2::aes_string(x = paste0("Prevalence_", group1),
                                           y = paste0("Prevalence_", group2))) +
      ggplot2::geom_point(alpha = 0.6, color = point_color, size = 3) +
      ggplot2::geom_smooth(method = "loess", color = smooth_color,
                  se = TRUE, fill = "indianred", alpha = 0.2) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
      ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      ggplot2::labs(
        x = paste("Prevalence in", group1),
        y = paste("Prevalence in", group2),
        title = paste("Prevalence Comparison:", group1, "vs", group2)
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.4),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
      )
    if(zoom_plot_flag != F) {
    # Zoom plot for top taxa
    zoom_plot <- ggplot2::ggplot(top_taxa, aes_string(x = paste0("Prevalence_", group1),
                                             y = paste0("Prevalence_", group2))) +
      ggplot2::geom_point(color = point_color, size = 4, alpha = 0.8) +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = sprintf("%s\n(%.1f%%, %.1f%%)", Label,
                            .data[[paste0("Prevalence_", group1)]] * 100,
                            .data[[paste0("Prevalence_", group2)]] * 100)),
        size = text_size,
        box.padding = 0.7,
        point.padding = 0.3,
        segment.size = 0.3,
        min.segment.length = 0,
        bg.color = "white",
        bg.r = 0.15,
        max.iter = 10000
      ) +
      ggplot2::labs(
        x = paste("Prevalence in", group1, "(zoom)"),
        y = paste("Prevalence in", group2, "(zoom)")
      ) +
      ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.4),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
      )


    # Combine plots
    combined_plot <- (main_plot / zoom_plot) +
      plot_layout(heights = c(1, 1))

    plot_list[[paste(group1, "vs", group2)]] <- combined_plot
  } else {
      plot_list[[paste(group1, "vs", group2)]] <- main_plot
      }
  }
  return(plot_list)
}
