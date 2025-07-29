#' Analyze Beta Diversity Across Filtering Parameters
#'
#' This function calculates beta diversity distances
#' across different detection/prevalence thresholds
#' and compares within-subject vs between-subject distances.
#'
#' @param ps_obj A phyloseq object
#' @param detection_levels Numeric vector of detection thresholds to test
#' @param prevalence_levels Numeric vector of
#' prevalence thresholds to test (0-1)
#' @param dist_method Distance method (e.g., "bray", "jaccard")
#' @param sample_id_col Column name containing sample IDs
#' @param subject_id_col Column name containing subject IDs
#' @param time_col Column name containing time points
#' @param within_label Label for within-subject comparisons
#' @param between_label Label for between-subject comparisons
#' @param color_palette Vector of two colors for plot
#'
#' @return A ggplot object of results
#' @export
analyze_beta_diversity_time_point <- function(ps_obj,
                                   detection_levels,
                                   prevalence_levels,
                                   dist_method = "bray",
                                   subject_id_col,
                                   time_col,
                                   within_label = "Within-subject",
                                   between_label = "Between-subjects",
                                   color_palette) {

  # Validate inputs
  if (!inherits(ps_obj, "phyloseq")) stop("Input must be a phyloseq object")
  required_cols <- c(subject_id_col, time_col)
  if (!all(required_cols %in% sample_variables(ps_obj))) {
    stop("Missing required columns in sample data: ",
         paste(setdiff(required_cols, sample_variables(ps_obj)), collapse = ", "))
  }
  core_taxa_counts <- tidyr::expand_grid(
    detection = detection_levels,
    prevalence = prevalence_levels
  ) %>%
    dplyr::mutate(
      n_taxa = purrr::map2_int(detection, prevalence, ~{
        core_taxa <- microbiome::core_members(ps_obj, detection = .x, prevalence = .y)
        length(core_taxa)
      })
    )
  # Main processing
  results <- tidyr::expand_grid(
    detection = detection_levels,
    prevalence = prevalence_levels
  ) %>%
    dplyr::left_join(core_taxa_counts, by = c("detection", "prevalence")) %>%
    dplyr::mutate(
      data = purrr::map2(detection, prevalence, ~{
        # Filter taxa
        core_taxa <- microbiome::core_members(ps_obj, detection = .x, prevalence = .y)
        if (length(core_taxa) == 0) return(NULL)
        ps_filtered <- phyloseq::prune_taxa(core_taxa, ps_obj)
        # Calculate distance matrix
        dist_matrix <- phyloseq::distance(ps_filtered, method = dist_method)
        # Convert to long format
        meta <- phyloseq::sample_data(ps_filtered) %>%
          methods::as("data.frame") %>%
          tibble::rownames_to_column("Sample") %>%
          dplyr::select(Sample, dplyr::all_of(c(subject_id_col, time_col)))

        methods::as(dist_matrix, "matrix") %>%
          base::as.data.frame() %>%
          tibble::rownames_to_column("Sample1") %>%
          tidyr::pivot_longer(-Sample1, names_to = "Sample2",
                              values_to = "distance") %>%
          dplyr::filter(Sample1 < Sample2) %>%
          dplyr::left_join(meta, by = c("Sample1" = "Sample")) %>%
          dplyr::left_join(meta, by = c("Sample2" = "Sample"),
                           suffix = c("_1", "_2")) %>%
          dplyr::mutate(
            comparison = ifelse(
              !!rlang::sym(paste0(subject_id_col, "_1")) ==
                !!rlang::sym(paste0(subject_id_col, "_2")),
              within_label,
              between_label
            )
          )
      })
    ) %>%
    tidyr::unnest(data) %>%
    dplyr::filter(!is.na(distance))

  # Calculate mean values for text output
  mean_values <-
    results %>%
    dplyr::group_by( detection, prevalence, comparison) %>%
    dplyr::summarise(
      mean_distance = base::mean(distance),
      sd_distance = stats::sd(distance),
      n_taxa = dplyr::first(n_taxa),
      .groups = "drop"
    )
  facet_labels <- tidyr::expand_grid(
    detection = detection_levels,
    prevalence = prevalence_levels
  ) %>%
    dplyr::left_join(core_taxa_counts, by = c("detection", "prevalence")) %>%
    dplyr::mutate(
      detection_label = paste0("Detection: ", detection),
      prevalence_label = paste0("Prevalence: ", prevalence, "\nTaxa: ", n_taxa)
    )
  # Create named vectors for labeller
  detection_labs <- setNames(facet_labels$detection_label,
                             facet_labels$detection)
  prevalence_labs <- setNames(facet_labels$prevalence_label,
                              facet_labels$prevalence)
  results <- results %>%
    mutate(
      facet_label = paste0(
        "Detection: ", detection,
        "\nPrevalence: ", prevalence,
        "\nTaxa: ", n_taxa
      )
    )

  p <- ggplot2::ggplot(results, ggplot2::aes(y = distance, x = comparison, fill = comparison)) +
    ggplot2::geom_violin(alpha = 0.7, trim = FALSE, scale = "width",
                         width = 0.8, linewidth = 0.8) +
    ggplot2::geom_boxplot(width = 0.15, color = "black", alpha = 0.8,
                          outlier.shape = NA, linewidth = 0.6) +
    ggplot2::facet_wrap(~ facet_label) +
    ggplot2::scale_fill_manual(
      values = color_palette,
      name = "Comparison Type"
    ) +
    ggplot2::labs(
      title = "Beta Diversity Distance Comparisons",
      x = "",
      y = "Bray-Curtis Distance",
      caption = paste("Analysis performed on", length(detection_levels),
                      "detection and", length(prevalence_levels),
                      "prevalence levels")
    ) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 10),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = ggplot2::element_text(size = 8),
      axis.title.y = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 8),
      legend.title = ggplot2::element_text(size = 10),
      plot.caption = ggplot2::element_text(size = 8)
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.2))) +
    ggpubr::stat_compare_means(
      comparisons = list(c("Within-subject", "Between-subjects")),
      method = "wilcox.test",
      size = 4,
      label.y = max(results$distance) * 1.1
    )

  results_table <- dplyr::select(mean_values, detection, prevalence, comparison, mean_distance, sd_distance) %>%
    dplyr::mutate(dplyr::across(c(mean_distance, sd_distance), ~base::round(., 2))) %>%
    tidyr::pivot_wider(
      names_from = comparison,
      values_from = c(mean_distance, sd_distance)
    ) %>%
    dplyr::arrange(detection, prevalence) %>%
    dplyr::rename_with(~stringr::str_remove(., "_distance"), dplyr::starts_with(c("mean_", "sd_")))

  #table_plot <- gridExtra::tableGrob(
  #  results_table,
  #  rows = base::rep("", nrow(results_table)),
  #  theme = gridExtra::ttheme_minimal(
  #    core = list(
  #      bg_params = list(fill = "white", col = "gray"),
  #      fg_params = list(hjust = 0, x = 0.05)
  #    ),
  #    colhead = list(
  #      bg_params = list(fill = "lightgray"),
  #      fg_params = list(fontface = "bold")
  #    )
  #  )
  #)

  return(list(violin_plot = p,
              results = results_table))
}
