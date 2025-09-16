#' @title Create Ordination Plots with Statistical Comparisons
#'
#' @description This function generates ordination plots (e.g., PCoA) from a phyloseq object
#' with statistical comparisons between groups. It creates a main ordination plot with
#' ellipses, group means, and ANOSIM results, along with marginal distribution plots
#' showing group differences along each axis.
#'
#' @param ps_obj A phyloseq object containing microbiome data
#' @param method Ordination method ("PCoA", "NMDS", etc.) (default: "PCoA")
#' @param distance_method Distance metric ("wunifrac", "bray", etc.) (default: "bray")
#' @param group Character string specifying the metadata column to group by
#' @param size Base font size for plot elements (default: 10)
#' @param palette Vector of colors for groups
#' @param taxa_plot Logical indicating whether to plot taxa vectors (default: TRUE)
#' @param taxrank Taxonomic level for taxa vector analysis ("Phylum", "Genus", etc.)
#'               (default: "Genus")
#' @param top_taxa Number of top significant taxa vectors to display
#'                (default: 10, set to NULL to show all significant taxa)
#' @param len_axis Scaling factor for taxa vectors (default: 1)
#' @param ratio Parameter indicating the ratio of figures
#' @param font_family Font family to use (default: "sans")
#'
#' @return A grid.arrange object containing the ordination plot with marginal distributions
#'
#' @details
#' The function performs the following analysis:
#' 1. Calculates distance matrix and performs ordination
#' 2. Computes ANOSIM to test group differences
#' 3. Plots ordination with ellipses and group means
#' 4. Adds significant taxa vectors (if taxa_plot = TRUE)
#' 5. Creates marginal violin plots with statistical comparisons
#'
#' When taxa_plot = TRUE, the function will:
#' - Aggregate taxa at the specified taxonomic rank (taxrank)
#' - Perform envfit analysis to find significant taxa vectors
#' - Display the top_taxa most significant vectors (by r value)
#' - Scale vectors by len_axis for better visualization
#' @examples
#' \dontrun{
#' # Example usage:
#' colors <- c("Group1" = "blue", "Group2" = "red")
#' p <- create_ordination_plots(
#'   ps_obj = my_phyloseq,
#'   method = "PCoA",
#'   distance_method = "bray",
#'   group = "Treatment",
#'   palette = colors
#' )
#' print(p)
#' }
#'
#' @importFrom phyloseq  sample_data tax_table otu_table sample_variables
#' @importFrom ggplot2 ggplot aes_string geom_point stat_ellipse geom_segment
#' @importFrom ggplot2 geom_label geom_text xlab ylab theme scale_color_manual
#' @importFrom ggplot2 element_text element_rect element_line unit ggplotGrob
#' @importFrom ggpubr compare_means stat_pvalue_manual
#' @importFrom rstatix add_significance add_x_position add_y_position
#' @importFrom gridExtra grid.arrange
#' @importFrom vegan anosim envfit
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats as.formula
#' @importFrom dplyr group_by summarise left_join
#' @importFrom grid unit.pmax arrow
#' @importFrom speedyseq tax_glom
#' @importFrom utils combn
#'
#' @export
create_ordination_plots <- function(ps_obj, method = "PCoA", distance_method = "bray",
                                    group, size = 10, palette, taxa_plot = TRUE,
                                    taxrank = "Genus", top_taxa = 10,
                                    len_axis = 1, ratio = c(1, 4),
                                    font_family = "sans") {

  # Validate input
  if (!inherits(ps_obj, "phyloseq")) stop("Input must be a phyloseq object")
  if (!group %in% colnames(sample_data(ps_obj))) {
    stop("Grouping variable not found in sample data")
  }


  # Automatically determine group levels (sorted)
  group_levels <- sort(unique(as.character(phyloseq::sample_data(ps_obj)[[group]])))
  phyloseq::sample_data(ps_obj)[[group]] <- factor(
    phyloseq::sample_data(ps_obj)[[group]],
    levels = group_levels
  )

  # Distance matrix and ordination
  dist <- phyloseq::distance(ps_obj, method = distance_method)
  ordination <- phyloseq::ordinate(ps_obj, method = method, distance = dist)

  # Prepare plot data
  P <- cbind(as.data.frame(ordination$vectors),
             as.data.frame(as.matrix(phyloseq::sample_data(ps_obj))))
  P[[group]] <- factor(P[[group]], levels = group_levels)

  # Calculate variance explained
  variance_explained <- 100 * ordination$values$Relative_eig[1:2]

  # Calculate group means
  means <- P %>%
    dplyr::group_by(!!rlang::sym(group)) %>%
    dplyr::summarise(mean_Axis1 = mean(Axis.1), mean_Axis2 = mean(Axis.2))

  P <- P %>% dplyr::left_join(means, by = group)

  # ANOSIM test
  anosim_res <- vegan::anosim(dist, P[[group]], parallel = 30)
  anosim_text <- paste("ANOSIM R:", round(anosim_res$statistic, 3),
                       "p-value:", round(anosim_res$signif, 3))

  # Taxa vectors (if enabled)
  if(taxa_plot) {
    ps_ph <- speedyseq::tax_glom(ps_obj, taxrank = taxrank, NArm = TRUE)
    env_fit <- vegan::envfit(ordination$vectors, as.data.frame(otu_table(ps_ph)),
                             permutations = 999, na.rm = TRUE)

    significant_vectors <- as.data.frame(env_fit$vectors$arrows) %>%
      dplyr::mutate(pvals = env_fit$vectors$pvals,
                    r = env_fit$vectors$r) %>%
      dplyr::filter(pvals < 0.05) %>%
      dplyr::arrange(desc(r)) %>%
      utils::head(top_taxa)

    tax_df <- as.data.frame(tax_table(ps_obj))[, c('Phylum','Genus'), drop = FALSE]
    top_significant_vectors <- merge(significant_vectors, tax_df, by = 0, all.x = TRUE)
    top_significant_vectors$taxa <- paste0(top_significant_vectors$Phylum, '_',
                                           top_significant_vectors$Genus)
  }

  # Main ordination plot
  pl <- ggplot(P, aes_string(x = "Axis.1", y = "Axis.2", color = group))  +
    geom_point(size = 1.5, alpha = 0.8) +
    stat_ellipse() +
    geom_segment(aes(xend = mean_Axis1, yend = mean_Axis2), alpha = 0.3) +
    geom_label(
      data = means,
      aes(x = mean_Axis1, y = mean_Axis2, label = !!rlang::sym(group)),
      inherit.aes = FALSE,  # Critical fix
      fill = "white",
      color = "black",
      fontface = "bold",
      size = 3
    ) +
    geom_text(x = Inf, y = Inf, label = anosim_text,
              hjust = 1.1, vjust = 1.1, size = size * 0.3) +
    scale_color_manual(values = palette) +
    xlab(paste0('Axis.1 [', round(variance_explained[1], 1), '%]')) +
    ylab(paste0('Axis.2 [', round(variance_explained[2], 1), '%]')) +
    theme_bw(base_size = size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, family = font_family, size = size),
      axis.text.y = element_text(family = font_family, size = size),
      axis.title = element_text(family = font_family, size = size, face = "bold"),
      legend.position = "none",
      panel.grid = element_line(size = 0.1, color = "gray"),
      text = element_text(family = font_family, size = size)
    )

  # Add taxa vectors if enabled
  if(taxa_plot && nrow(top_significant_vectors) > 0) {
    pl <- pl +
      geom_segment(
        data = top_significant_vectors,
        aes(x = 0, y = 0, xend = Axis.1 * len_axis, yend = Axis.2 * len_axis),
        arrow = arrow(length = unit(0.2, "cm")), color = "black", alpha = 0.8
      ) +
      ggrepel::geom_text_repel(
        data = top_significant_vectors,
        aes(x = Axis.1 * len_axis, y = Axis.2 * len_axis, label = taxa),
        color = "black", size = 3, box.padding = 0.5
      )+
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, family = font_family, size = size),
        axis.text.y = element_text(family = font_family, size = size),
        axis.title = element_text(family = font_family, size = size, face = "bold"),
        text = element_text(family = font_family, size = size)
      )
  }

  # Create comparison list for stats
  comparisons <- combn(group_levels, 2, simplify = FALSE)

  # Axis.1 statistics and plot
  stats_axis1 <- ggpubr::compare_means(
    as.formula(paste("Axis.1 ~", group)),
    data = P,
    method = "wilcox.test",
    p.adjust.method = "BH",
    comparisons = comparisons
  ) %>%
    rstatix::add_significance("p.adj") %>%
    rstatix::add_y_position(data = P,
                            formula = as.formula(paste("Axis.1 ~", group)),
                            step.increase = 0.2)

  # Convert group names to numeric positions
  stats_axis1$xmin <- as.numeric(factor(stats_axis1$group1, levels = group_levels))
  stats_axis1$xmax <- as.numeric(factor(stats_axis1$group2, levels = group_levels))
  print(stats_axis1)

  x_dens <- ggplot(P, aes_string(x = group, y = "Axis.1", color = group)) +
    geom_violin(trim = FALSE, alpha = 0.1) +
    geom_boxplot(width = 0.5, alpha = 0.75) +
    geom_jitter(size = 1.5, alpha = 0.5, width = 0.1) +
    scale_color_manual(values = palette) +
    ggpubr::stat_pvalue_manual(
      stats_axis1, label = "p.adj.signif",
      tip.length = 0.02, step.increase = 0.05,coord.flip = TRUE
    ) +
    coord_flip() +
    theme_bw(base_size = size) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, family = font_family, size = size),
          axis.text.y = element_text(family = font_family, size = size),
          axis.title = element_text(family = font_family, size = size, face = "bold"))

  # Axis.2 statistics and plot
  stats_axis2 <- ggpubr::compare_means(
    as.formula(paste("Axis.2 ~", group)),
    data = P,
    method = "wilcox.test",
    p.adjust.method = "BH",
    comparisons = comparisons
  ) %>%
    rstatix::add_significance("p.adj") %>%
    rstatix::add_y_position(data = P,
                            formula = as.formula(paste("Axis.2 ~", group)),
                            step.increase = 0.2)


  stats_axis2$xmin <- as.numeric(factor(stats_axis2$group1, levels = group_levels))
  stats_axis2$xmax <- as.numeric(factor(stats_axis2$group2, levels = group_levels))
  print(stats_axis2)

  y_dens <- ggplot(P, aes_string(x = group, y = "Axis.2", color = group)) +
    geom_violin(trim = FALSE, alpha = 0.1) +
    geom_boxplot(width = 0.5, alpha = 0.75) +
    geom_jitter(size = 1.5, alpha = 0.5, width = 0.1) +
    scale_color_manual(values = palette) +
    ggpubr::stat_pvalue_manual(
      stats_axis2, label = "p.adj.signif",
      tip.length = 0.02, step.increase = 0.05
    ) +
    theme_bw(base_size = size) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, family = font_family, size = size),
          axis.text.y = element_text(family = font_family, size = size),
          axis.title = element_text(family = font_family, size = size, face = "bold"))

  # Legend plot
  legend_plot <- ggplot(P, aes_string(x = "Axis.1", y = "Axis.2", color = group)) +
    geom_point(size = 2.5) +
    scale_color_manual(values = palette) +
    theme_bw() +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 45, hjust = 1, family = font_family, size = size),
          axis.text.y = element_text(family = font_family, size = size),
          axis.title = element_text(family = font_family, size = size, face = "bold"))
  legend <- cowplot::get_legend(legend_plot)

  # Arrange plots
  final_plot <- gridExtra::grid.arrange(
    y_dens, pl, legend, x_dens,
    ncol = 2, nrow = 2,
    widths = ratio, heights = rev(ratio)
  )
  return(final_plot)
}

