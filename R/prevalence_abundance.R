#' Generate Prevalence-Abundance Relationship Plot
#'
#' Creates a scatter plot showing the relationship
#' between taxon prevalence
#' and abundance (mean relative abundance).
#' Includes options for zoomed-in view of top taxa.
#'
#' @param phyloseq A phyloseq object containing OTU/ASV table and taxonomy
#' @param top_n Number of top taxa to highlight (default: 15)
#' @param zoom_plot Logical, whether to include zoomed plot (default: TRUE)
#' @param output Output plot file folder path
#' @param point_color Color for coloring scatter points
#' @param line_color Color for coloring geom_smooth
#' @param method method for geom_smooth
#' @param formula formula for geom_smooth
#' @return A ggplot object or objects if zoom_plot == TRUE
#' @examples
#' # p <- prevalence_abundance_plot(phyloseq, top_n = 10, zoom_plot = TRUE)
#' @export
prevalence_abundance_plot <- function(phyloseq, top_n = 15, zoom_plot = TRUE,
                                      line_color = NULL,
                                      point_color = NULL,
                                      method = "loess",
                                      formula = y ~ x) {
  # Validate input
  if (!inherits(phyloseq, "phyloseq")) {
    stop("Input must be a phyloseq object")
  }

  # Convert OTU table to matrix
  otu_mat <- methods::as(phyloseq::otu_table(phyloseq), "matrix")
  if (!phyloseq::taxa_are_rows(phyloseq)) {
    otu_mat <- base::t(otu_mat)
  }

  # Check if OTU table is empty
  if (base::nrow(otu_mat) == 0 || base::ncol(otu_mat) == 0) {
    stop("OTU table is empty - no taxa or samples found")
  }

  # Calculate metrics
  occupancy <- base::rowSums(otu_mat > 0) / base::ncol(otu_mat)
  mean_relative_abundance <- base::rowMeans(otu_mat)
  log_mean_relative_abundance <- base::log10(mean_relative_abundance + 1e-5)

  # Create data frame
  df <- base::data.frame(
    Occupancy = occupancy,
    Abundance = log_mean_relative_abundance,
    stringsAsFactors = FALSE
  )

  # Remove NA/NaN/Inf values
  df <- df[base::is.finite(df$Occupancy) & base::is.finite(df$Abundance), ]

  # Check if any data remains
  if (base::nrow(df) == 0) {
    stop("No valid taxa remaining after filtering NA/Inf values")
  }

  # Add taxonomy info if available
  if (!base::is.null(phyloseq::tax_table(phyloseq, errorIfNULL = FALSE))) {
    taxa_names <- base::as.data.frame(phyloseq::tax_table(phyloseq))
    df$Taxon <- taxa_names$Species
  } else {
    df$Taxon <- base::rownames(df)
  }

  # Identify top taxa
  top_taxa <- df %>%
    dplyr::mutate(
      occupancy_rank = base::rank(-Occupancy),
      abundance_rank = base::rank(-Abundance),
      combined_rank = occupancy_rank + abundance_rank
    ) %>%
    dplyr::arrange(combined_rank) %>%
    dplyr::slice_head(n = top_n)
  if (is.null(line_color)) {
    line_color <- "#d62828"
  }
  if (is.null(point_color)) {
    point_color <- "#d5b3ff"
  }
  # Create main plot
  main_plot <- ggplot2::ggplot(df, ggplot2::aes(x = Abundance, y = Occupancy)) +
    ggplot2::geom_point(alpha = 0.6, color = point_color, size = 3) +
    ggplot2::geom_smooth(
      method = method,
      formula = formula,
      color = line_color,
      se = TRUE,
      fill = "indianred",
      alpha = 0.2
    ) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      x = expression(paste("Mean relative abundance (log"[10], ")")),
      y = "Prevalence",
      title = "Prevalence-Abundance Relationship"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "gray90", linewidth = 0.2),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.4),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "gray40"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5)
    )

  if (!zoom_plot) {
    # Return just main plot with annotations
    main_plot <- main_plot +
      ggplot2::geom_vline(
        data = top_taxa,
        ggplot2::aes(xintercept = min(Abundance)),
        color = "gray40",
        linetype = "dashed",
        linewidth = 0.4
      ) +
      ggplot2::geom_hline(
        data = top_taxa,
        ggplot2::aes(yintercept = min(Occupancy)),
        color = "gray40",
        linetype = "dashed",
        linewidth = 0.4
      ) +
      ggrepel::geom_text_repel(
        data = top_taxa,
        ggplot2::aes(label = Taxon),
        size = 3,
        box.padding = 0.3,
        point.padding = 0.3,
        segment.size = 0.3,
        min.segment.length = 0.1,
        bg.color = "white",
        bg.r = 0.1,
        max.overlaps = 20
      )
    return(list(plot = main_plot,df_top = top_taxa))
  } else {
    # Create zoomed plot
    zoom_plot <- ggplot2::ggplot(top_taxa, ggplot2::aes(x = Abundance, y = Occupancy)) +
      ggplot2::geom_point(color = point_color, size = 4, alpha = 0.8) +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = Taxon),
        size = 3,
        box.padding = 0.3,
        point.padding = 0.3,
        segment.size = 0.3,
        min.segment.length = 0.1,
        bg.color = "white",
        bg.r = 0.1,
        max.overlaps = 20
      ) +
      ggplot2::labs(
        x = expression(paste("Mean relative abundance (log"[10], ") zoom")),
        y = "Prevalence zoom(%)",
        title = "Top Taxa by Prevalence and Abundance"
      ) +
      ggplot2::scale_y_continuous(
        labels = scales::percent_format(accuracy = 1),
        breaks = scales::pretty_breaks(n = 4)
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(color = "gray90", linewidth = 0.2),
        panel.grid.minor = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(color = "black", linewidth = 0.4),
        axis.title = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5)
      )

    # Combine plots using patchwork
    combined_plot <- patchwork::wrap_plots(main_plot, zoom_plot, ncol = 1)
    return(list(plot = combined_plot,df_top = top_taxa))
  }
}
