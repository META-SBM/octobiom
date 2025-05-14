#' Build a Feature-Taxon Correlation Network
#'
#' This function creates a network visualization showing significant correlations between
#' microbial taxa and clinical/environmental features.
#'
#' @param phyloseq_obj A phyloseq object containing OTU/ASV table and sample metadata
#' @param features Character vector of metadata features to analyze (e.g., c("age", "bmi"))
#' @param tax_rank Taxonomic rank to use for analysis ("Phylum", "Genus", "Species", etc.)
#' @param cor_method Correlation method ("pearson", "spearman", or "kendall")
#' @param p_adjust_method P-value adjustment method ("BH", "bonferroni", etc.)
#' @param p_cutoff Significance threshold for adjusted p-values
#' @param edge_width_range Numeric vector specifying min/max edge widths (e.g., c(0.8, 3))
#' @param node_size_range Numeric vector specifying min/max node sizes (e.g., c(4, 15))
#' @param layout Network layout algorithm ("fr", "kk", "dh", etc.)
#' @param feature_colors Optional named vector of colors for features (e.g., c("age" = "blue", "bmi" = "red"))
#' @param edge_palette Optional vector of colors for edges (positive/negative correlations)
#' @param min_corr_strenght Optional value for filtering out weak correlations (default: 0)
#'
#' @return A ggraph object visualizing the feature-taxon correlation network
#'
#' @examples
#' \dontrun{
#' data(GlobalPatterns)
#' net <- build_multi_correlation_network(
#'   phyloseq_obj = GlobalPatterns,
#'   features = c("SampleType", "pH"),
#'   tax_rank = "Genus"
#' )
#' print(net)
#' }
#'
#' @importFrom phyloseq otu_table sample_data tax_table
#' @importFrom igraph graph_from_data_frame V E
#' @importFrom ggraph ggraph geom_edge_arc geom_node_point geom_node_text geom_node_label
#' @importFrom dplyr select mutate group_by group_modify ungroup filter distinct summarise
#' @importFrom tidyr pivot_longer
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble
#' @importFrom stats cor.test p.adjust
#' @importFrom scales hue_pal
#' @export
build_multi_correlation_network <- function(
    phyloseq_obj,
    features = c("age_at_baseline", "bmi"),
    tax_rank = "Species",
    cor_method = "spearman",
    p_adjust_method = "BH",
    p_cutoff = 0.05,
    edge_width_range = c(0.8, 3),
    node_size_range = c(4, 15),
    layout = "fr",
    feature_colors,
    edge_palette = c(positive = "#ff1b6b", negative = "#45caff"),
    min_corr_strenght = 0
) {
  # 1. Data extraction and preparation
  otu <- as(phyloseq::otu_table(phyloseq_obj), 'matrix') %>%
    as.data.frame()
  meta <- as(phyloseq::sample_data(phyloseq_obj), 'data.frame')[, features, drop = FALSE]
  tax <- as.data.frame(phyloseq::tax_table(phyloseq_obj))[[tax_rank]]

  # 2. Convert to long format
  analysis_data <- cbind(meta, otu) %>%
    pivot_longer(
      cols = all_of(colnames(otu)),
      names_to = "taxon",
      values_to = "abundance"
    ) %>%
    mutate(taxon = tax[match(taxon, rownames(tax_table(phyloseq_obj)))])

  # 3. Calculate correlations
  results <- analysis_data %>%
    group_by(taxon) %>%
    group_modify(~ {
      cors <- map_dfr(features, function(feature) {
        test <- cor.test(.x$abundance, .x[[feature]], method = cor_method)
        tibble(
          feature = feature,
          cor = test$estimate,
          p_value = test$p.value,
          mean_abundance = mean(.x$abundance, na.rm = TRUE)
        )
      })
    }) %>%
    ungroup() %>%
    mutate(p_adj = p.adjust(p_value, method = p_adjust_method)) %>%
    filter(p_adj < p_cutoff, !is.na(cor)) %>%
    filter(abs(cor) > min_corr_strenght, !is.na(cor)) %>%
    mutate(
      direction = ifelse(cor > 0, "positive", "negative"),
      taxon = factor(taxon)
    )

  if (nrow(results) == 0) stop("No significant correlations found")

  # 4. Create graph structure
  results_df <- results %>%
    select(feature, taxon, cor, direction, mean_abundance)

  g <- igraph::graph_from_data_frame(results_df, directed = FALSE)

  # 5. Vertex attributes - modified version
  g <- igraph::set_vertex_attr(g, "type",
                               value = ifelse(igraph::V(g)$name %in% features, "Feature", "Taxa"))
  all_taxa_in_graph <- igraph::vertex_attr(g, "name")[igraph::vertex_attr(g, "type") == "Taxa"]
  missing_taxa <- setdiff(all_taxa_in_graph, results_df$taxon)
  if(length(missing_taxa) > 0) message("Lost taxa: ", paste(missing_taxa, collapse = ", "))

  taxon_sizes <- results_df %>%
    dplyr::group_by(taxon) %>%
    dplyr::summarise(
      size = mean(mean_abundance, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    dplyr::right_join(
      tibble::tibble(taxon = all_taxa_in_graph),
      by = "taxon"
    ) %>%
    dplyr::mutate(
      size = dplyr::case_when(
        is.na(size) ~ min(size, na.rm = TRUE),
        size <= 0 ~ 0.1, # Заменяем отрицательные и нулевые значения
        TRUE ~ size
      )
    )

  # Set feature node size to mean of taxa sizes
  g <- igraph::set_vertex_attr(g, "size",
                               index = igraph::V(g)[type == "Taxa"],
                               value = taxon_sizes$size)
  feature_size <- mean(igraph::vertex_attr(g, "size")[igraph::vertex_attr(g, "type") == "Taxa"], na.rm = TRUE)
  g <- igraph::set_vertex_attr(
    g,
    "size",
    index = which(igraph::vertex_attr(g, "type") == "Feature"),
    value = feature_size
  )
  # Color handling
  if (is.null(feature_colors)) {
    feature_colors <- scales::hue_pal()(length(features))
    names(feature_colors) <- features
  }

  # Create taxon-to-feature color mapping
  taxon_to_feature <- results_df %>%
    dplyr::distinct(taxon, feature) %>%
    dplyr::mutate(color = feature_colors[feature])

  vertex_colors <- ifelse(
    igraph::vertex_attr(g, "type") == "Feature",
    feature_colors[igraph::vertex_attr(g, "name")],
    taxon_to_feature$color[match(igraph::vertex_attr(g, "name"), taxon_to_feature$taxon)]
  )

  g <- igraph::set_vertex_attr(g, "color", value = vertex_colors)

  # 6. Edge attributes - переписанный блок
  g <- igraph::set_edge_attr(g, "edge_direction", value = results_df$direction)
  g <- igraph::set_edge_attr(g, "edge_cor_abs", value = abs(results_df$cor))
  g <- igraph::set_edge_attr(g, "edge_cor_value", value = results_df$cor)

  # 7. Visualization
  set.seed(123)
  p <- ggraph::ggraph(g, layout = layout) +
    ggraph::geom_edge_arc(
      ggplot2::aes(
        width = edge_cor_abs,
        color = edge_direction
      ),
      strength = 0.1,
      end_cap = ggraph::circle(5, "mm")) +
    ggraph::geom_node_point(
      ggplot2::aes(size = size, fill = color),
      shape = 21,
      color = "black") +
    ggraph::geom_node_text(
      ggplot2::aes(label = name, filter = type == "Taxa"),
      repel = TRUE,
      size = 3.5) +
    ggraph::geom_node_label(
      ggplot2::aes(label = name, filter = type == "Feature"),
      nudge_y = 0.1,
      size = 4.5) +
    ggraph::scale_edge_width_continuous(
      range = edge_width_range,
      name = "Correlation strength") +
    ggraph::scale_edge_color_manual(
      values = edge_palette,
      name = "Direction") +
    ggplot2::scale_fill_identity(
      guide = "legend",
      labels = c(names(feature_colors), "Taxa"),
      breaks = c(feature_colors, "grey70"),
      name = "Feature/Taxa"
    ) +
    ggplot2::scale_size_continuous(
      range = node_size_range,
      limits = c(0, max(taxon_sizes$size, na.rm = TRUE)),
      name = "Mean abundance"
    )+
    ggplot2::theme_void() +
    ggplot2::labs(
      title = "Feature-Taxon Correlation Network",
      subtitle = paste("Taxonomic level:", tax_rank, "| Correlation method:", cor_method)
    )

  return(list(results_df,p))
}
