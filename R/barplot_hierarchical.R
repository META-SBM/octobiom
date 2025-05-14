#' Create Taxonomy Table for "Other" Category
#'
#' Generates a taxonomy table row for aggregated "other" taxa
#'
#' @param colnames Character vector of taxonomic rank names
#' (e.g., c("Kingdom","Phylum","Class"))
#' @return A data.frame formatted as a taxonomy table row for "other" category
#' @examples
#' tax_table_other(c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
#' @export
tax_table_other <- function(colnames) {
  df <- data.frame(matrix('other', nrow = 1, ncol = length(colnames)))
  colnames(df) <- colnames
  rownames(df) <- 'other'
  return(df)
}

#' @title Hierarchical Composition Barplot with Dendrogram
#'
#' @description Creates a stacked barplot of taxonomic compositions
#' with hierarchical clustering dendrogram and sample metadata annotations.
#'
#' @param ps A phyloseq object containing OTU table and taxonomy table
#' @param taxrank Taxonomic rank to use for coloring barplot (e.g., "Genus")
#' @param feature Metadata feature to annotate in color strips
#' @param top Number of top taxa to show individually
#' @param colors_barplot Color palette for taxonomic groups
#' @param colors_annotation List of color palettes for each annotation feature
#' @param dist Distance method for clustering ("bray", "aitchinson")
#' @param size Relative heights for plot components
#' (dendrogram, annotation, barplot)
#' @param decreasing Logical for ordering "other" category position
#'
#' @return A patchwork composite ggplot object
#'
#' @examples
#' \dontrun{
#' data(GlobalPatterns)
#' tax_cols <- colnames(phyloseq::tax_table(GlobalPatterns))
#' other_tax <- tax_table_other(tax_cols)
#'
#' barplot_hierarchical(
#'   ps = GlobalPatterns,
#'   taxrank = "Genus",
#'   feature = "SampleType",
#'   top = 10,
#'   colors_barplot = c("steelblue", "darkorange"),
#'   colors_line = c("red", "blue"),
#'   dist = "bray",
#'   size = c(1, 0.2, 2),
#'   decreasing = TRUE
#' )
#' }
#'
#' @importFrom patchwork plot_layout
#' @export
barplot_hierarchical <- function(ps, taxrank,
                                 feature, top, colors_barplot, colors_annotation,
                                 dist, size,
                                 decreasing) {
  # Validate inputs
  if (!inherits(ps, "phyloseq")) stop("Input must be a phyloseq object")
  if (!taxrank %in% phyloseq::rank_names(ps)) stop("Specified taxonomic rank not found")

  # Extract metadata for annotation
  metadata <- ps@sam_data[, feature, drop = FALSE]

  # Prepare taxonomy matrix
  taxa_matrix <- as.data.frame(phyloseq::tax_table(ps))
  taxa_matrix <- cbind(ASV = rownames(taxa_matrix), taxa_matrix)

  # Create taxonomy table for 'other' category
  tax_cols <- colnames(phyloseq::tax_table(ps))
  tax_table_other <- tax_table_other(tax_cols)

  # Select top taxa and create separate objects
  topx <- microbiome::top_taxa(ps, n = top)
  not_topx <- !(taxa_matrix$ASV %in% topx)

  # Create phyloseq objects
  ps_top <- phyloseq::prune_taxa(topx, ps)
  ps_top@sam_data <- metadata

  ps_other <- phyloseq::prune_taxa(not_topx, ps)

  # Aggregate remaining taxa as "other"
  otu_table_other <- as.matrix(rowSums(phyloseq::otu_table(ps_other)))
  colnames(otu_table_other) <- "other"

  # Create phyloseq object for "other" category
  TAX <- phyloseq::tax_table(as.matrix(tax_table_other))
  OTU <- phyloseq::otu_table(otu_table_other, taxa_are_rows = FALSE)
  samples <- phyloseq::sample_data(metadata)
  ps_other <- phyloseq::phyloseq(OTU, TAX, samples)

  # Merge data for plotting
  data_plot <- phyloseq::psmelt(ps_top)
  data_plot_other <- phyloseq::psmelt(ps_other)
  colnames_to_keep <- colnames(data_plot)
  data_plot_other <- data_plot_other[, colnames_to_keep]
  df <- rbind(data_plot, data_plot_other)

  # Create dendrogram based on distance metric
  if (dist == 'bray') {
    dend <- stats::as.dendrogram(
      stats::hclust(
        vegan::vegdist(ps@otu_table, method = "bray")
      )
    )
  }
  if (dist == 'aitchinson') {
    dend <- stats::as.dendrogram(stats::hclust(stats::dist(ps@otu_table)))
  }

  # Prepare plot components
  hc <- ggdendro::dendro_data(stats::as.dendrogram(dend))
  df[, taxrank] <- factor(
    df[, taxrank],
    levels = unique(df[, taxrank][order(df[, taxrank] == 'other', decreasing = decreasing)]))

    df$Sample <- factor(df$Sample, levels = ggdendro::label(hc)$label)

  # Build plot components
  p1 <- ggplot2::ggplot(data = ggdendro::segment(hc)) +
    ggplot2::geom_segment(ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
    ggplot2::scale_x_discrete(labels = ggdendro::label(hc)$label) +
    ggplot2::ylab(dist) + ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank()
    )
  p2 <- ggplot2::ggplot(df, ggplot2::aes(fill = !!sym(taxrank), y = Abundance, x = Sample)) +
    ggplot2::geom_bar(stat = "identity", position = "fill", width = 1,
                      color = "black", size = 0.1) +
    ggplot2::scale_fill_manual("legend", values = colors_barplot, name = taxrank) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(size = 1, colour = "black"),
      axis.title.x = ggplot2::element_text(size = 10, hjust = 1),
      legend.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.position = 'right',
      legend.box = "vertical"
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))

  p3 <- lapply(features, function(feat) {
    ggplot2::ggplot(df, ggplot2::aes(x = Sample, y = 1, fill = !!sym(feat))) +
      ggplot2::geom_tile(color = 'black', size = 0.1) +
      ggplot2::scale_fill_manual(values = colors_annotation[[feat]]) +
      ggplot2::theme_void() +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 10),
        legend.title = ggplot2::element_text(size = 10, face = "bold"),
        legend.justification = "left"
      )+
      ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1,title.position = 'left'))
  })

  plot_list <- c(list(p1), p3,list(p2))
  combined_plot <- patchwork::wrap_plots(plot_list, ncol = 1, heights = size)
  return(combined_plot)
}
