#' Calculate and Visualize Taxa Retention Across Detection/Prevalence Thresholds
#'
#' This function evaluates how many microbial taxa are
#' retained under different combinations
#' of detection and prevalence thresholds, and visualizes
#' the results as a line plot.
#'
#' @param phyloseq A phyloseq object containing microbiome data
#' @param detection_values Numeric vector of detection
#' thresholds to test (0-100 for percentages)
#' @param prevalence_values Numeric vector of prevalence
#' thresholds to test (0-1 for proportions)
#' @param prevalence_colors Character vector of colors for each prevalence level
#'
#' @return A ggplot object showing taxa retention patterns
#' @export
#'
#' @examples
#' # Example with mock data
#' data(GlobalPatterns)
#' det_vals <- seq(0, 0.1, length.out = 5)
#' prev_vals <- seq(0.1, 0.9, by = 0.2)
#' colors <- viridis::viridis(length(prev_vals))
#' p <- retention_plot(
#'   phyloseq = GlobalPatterns,
#'   detection_values = det_vals,
#'   prevalence_values = prev_vals,
#'   taxonomic_level = "Genus",
#'   prevalence_colors = colors
#' )
#' print(p)
retention_plot <- function(phyloseq,
                                  detection_values,
                                  prevalence_values,
                                  prevalence_colors = NULL) {
  # Input validation
  if (!inherits(phyloseq, "phyloseq")) {
    stop("Input must be a phyloseq object.")
  }
  if (any(prevalence_values < 0 | prevalence_values > 1)) {
    stop("Prevalence values must be between 0 and 1.")
  }
  if (length(prevalence_colors) < length(unique(prevalence_values))) {
    stop("Insufficient set of colors.")
  }
  # Initialize results dataframe
  results <- base::expand.grid(
    detection = detection_values,
    prevalence = prevalence_values,
    n_taxa = NA_integer_
  )
  # Calculate retention for each combination
  for (i in seq_len(nrow(results))) {
    n_taxa <- tryCatch({
      core_taxa <- microbiome::core_members(
        phyloseq,
        detection = results$detection[i],
        prevalence = results$prevalence[i]
      )
      length(core_taxa)
    }, error = function(e) NA)
    results$n_taxa[i] <- n_taxa
  }
  if (is.null(prevalence_colors)) {
    n_colors <- length(unique(prevalence_values))
    prevalence_colors <- scales::hue_pal()(n_colors)
    names(prevalence_colors) <- sort(unique(prevalence_values))
  }
  # Prepare plotting data
  results$prevalence_factor <- factor(results$prevalence,
                                      levels = sort(unique(prevalence_values)))
  # Generate plot
  p <- ggplot2::ggplot(results, ggplot2::aes(x = detection * 100, y = n_taxa,
                           color = prevalence_factor)) +
    ggplot2::geom_point(size = 3, na.rm = TRUE) +
    ggplot2::geom_line(aes(group = prevalence_factor), linewidth = 1) +
    ggplot2::geom_label(
      ggplot2::aes(label = ifelse(is.na(n_taxa), "", n_taxa)),
      size = 2.5,
      color = "black",
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(values = prevalence_colors, name = "Prevalence") +
    ggplot2::labs(
      x = "Detection Threshold (%)",
      y = "Number of remaining bacteria"
    ) +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey80", linewidth = 0.2),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.1),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.title = element_text(face = "bold", size = 8),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 6),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  return(p)
}
