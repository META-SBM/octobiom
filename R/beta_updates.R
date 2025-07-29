#' @title Create Ordination Plots with Statistical Comparisons (Dual Group Version)
#'
#' @description This function generates ordination plots (e.g., PCoA) from a phyloseq object
#' with statistical comparisons between groups. It creates a main ordination plot with
#' ellipses, group means, and ANOSIM results, along with marginal distribution plots
#' showing group differences along each axis. Supports dual grouping variables.
#'
#' @param ps_obj A phyloseq object containing microbiome data
#' @param method Ordination method ("PCoA", "NMDS", etc.) (default: "PCoA")
#' @param distance_method Distance metric ("wunifrac", "bray", etc.) (default: "bray")
#' @param group Character string specifying the primary metadata column to group by
#' @param size Base font size for plot elements (default: 10)
#' @param palette Vector of colors for primary groups
#' @param taxa_plot Logical indicating whether to plot taxa vectors (default: FALSE)
#' @param taxrank Taxonomic level for taxa vector analysis ("Phylum", "Genus", etc.)
#'               (default: "Genus")
#' @param top_taxa Number of top significant taxa vectors to display
#'                (default: 10)
#' @param len_axis Scaling factor for taxa vectors (default: 1)
#' @param ratio Parameter indicating the ratio of figures (default: c(1, 3))
#' @param additional_group Secondary grouping variable (default: 'NONE')
#' @param shape_palette Vector of shapes for secondary grouping variable
#' @param additional_palette Vector of colors for secondary groups
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
#' When two grouping variables are provided:
#' - Primary variable represented by color
#' - Secondary variable represented by shape
#' - Separate violins shown for each variable
#'
#' @examples
#' \dontrun{
#' # Example with two grouping variables:
#' p <- create_ordination_plots_double(
#'   ps_obj = my_phyloseq,
#'   group = "BMI",
#'   additional_group = "Treatment",
#'   palette = c("Normal" = "green", "Obese" = "red"),
#'   additional_palette = c("Control" = "gray", "Treated" = "blue"),
#'   shape_palette = c(16, 17)
#' )
#' }
#'
#' @importFrom phyloseq sample_data tax_table otu_table distance ordinate
#' @importFrom ggplot2 ggplot aes geom_point stat_ellipse geom_segment geom_label
#' @importFrom ggplot2 geom_text xlab ylab theme scale_color_manual scale_shape_manual
#' @importFrom ggplot2 theme_bw element_text element_rect element_line unit
#' @importFrom ggpubr compare_means stat_pvalue_manual ggarrange
#' @importFrom rstatix add_significance add_y_position
#' @importFrom gridExtra grid.arrange
#' @importFrom vegan anosim envfit
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr group_by summarise left_join
#' @importFrom speedyseq tax_glom
#' @importFrom cowplot get_legend
#' @importFrom stats as.formula
#' @importFrom utils combn
#' @importFrom rlang sym !! .data
#' @export
create_ordination_plots_double <- function(ps_obj, method = "PCoA", distance_method = "bray",
                                    group = 'bmi_group', size = 10, palette, taxa_plot = FALSE,
                                    taxrank = "Genus", top_taxa = 10,
                                    len_axis = 1, ratio = c(1, 3),additional_group = NULL ,shape_palette= NULL,
                                    additional_palette= NULL) {

  # Validate input
  if (!inherits(ps_obj, "phyloseq")) stop("Input must be a phyloseq object")
  if (!group %in% colnames(sample_data(ps_obj))) {
    stop("Grouping variable not found in sample data")
  }


  # Automatically determine group levels (sorted)
  group_levels <- levels(phyloseq::sample_data(ps_obj)[[group]])
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
  pl <- ggplot(P, aes(x = `Axis.1`, y = `Axis.2`, color = !!rlang::sym(group)))  +
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
      alpha=0.4,
      size = 3
    ) +
    scale_color_manual(values = palette) +
    xlab(paste0('Axis.1 [', round(variance_explained[1], 1), '%]')) +
    ylab(paste0('Axis.2 [', round(variance_explained[2], 1), '%]')) +
    labs(subtitle = anosim_text)+
    #geom_text(x = Inf, y = Inf, label = anosim_text,
    #          hjust = 1.1, vjust = 1.1, size = size * 0.3) +
    theme_bw(base_size = size) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.position = "right",
      legend.background = element_rect(color = "black"),
      legend.direction = "horizontal",
      panel.grid = element_line(size = 0.1, color = "gray")
    )


  if(!is.null(additional_group)){
    if (!additional_group %in% colnames(sample_data(ps_obj))) {
      stop("Grouping variable not found in sample data")
    }
    means <- P %>%
      dplyr::group_by(!!rlang::sym(additional_group)) %>%
      dplyr::summarise(mean_Axis1 = mean(Axis.1), mean_Axis2 = mean(Axis.2))
    pl <- pl +
      geom_point(
        aes(shape = !!sym(additional_group)))+
      scale_shape_manual(values = shape_palette)+
      geom_label(
        data = means,
        aes(x = mean_Axis1, y = mean_Axis2, label = !!rlang::sym(additional_group)),
        inherit.aes = FALSE,  # Critical fix
        fill = "white",
        color = "black",
        fontface = "bold",
        alpha=0.4,
        size = 3
      )

  }
  legend <- cowplot::get_legend(pl)
  pl <- pl + theme(legend.position = "none")


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
      )
  }
  group_levels <- levels(phyloseq::sample_data(ps_obj)[[group]])
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
                            step.increase = 0.05)

  # Convert group names to numeric positions
  stats_axis1$xmin <- as.numeric(factor(stats_axis1$group1, levels = group_levels))
  stats_axis1$xmax <- as.numeric(factor(stats_axis1$group2, levels = group_levels))
  print(stats_axis1)

  x_dens <- ggplot(P, aes_string(x = group, y = "Axis.1", color = group)) +
    #geom_violin(trim = FALSE, alpha = 0.1) +
    geom_boxplot(width = 0.5, alpha = 0.75) +
    geom_jitter(size = 1.5, alpha = 0.5, width = 0.1) +
    scale_color_manual(values = palette) +
    ggpubr::stat_pvalue_manual(
      stats_axis1, label = "p.adj.signif",
      tip.length = 0.02, step.increase = 0.05,coord.flip = TRUE
    ) +
    coord_flip() +
    theme_bw(base_size = size) +
    theme(legend.position = "none")

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
                            step.increase = 0.05)


  stats_axis2$xmin <- as.numeric(factor(stats_axis2$group1, levels = group_levels))
  stats_axis2$xmax <- as.numeric(factor(stats_axis2$group2, levels = group_levels))
  print(stats_axis2)

  y_dens <- ggplot(P, aes_string(x = group, y = "Axis.2", color = group)) +
    #geom_violin(trim = FALSE, alpha = 0.1) +
    geom_boxplot(width = 0.5, alpha = 0.75) +
    geom_jitter(size = 1.5, alpha = 0.5, width = 0.1) +
    scale_color_manual(values = palette) +
    ggpubr::stat_pvalue_manual(
      stats_axis2, label = "p.adj.signif",
      tip.length = 0.02, step.increase = 0.05
    ) +
    theme_bw(base_size = size) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45))
  pl_y_range <- ggplot_build(y_dens)$layout$panel_params[[1]]$y.range
  y_range_pl <- ggplot_build(pl)$layout$panel_params[[1]]$y.range

  y_min <- min(pl_y_range[1], y_range_pl[1])
  y_max <- max(pl_y_range[2], y_range_pl[2])
  common_y_range <- c(y_min, y_max)

  y_dens <- y_dens + ggplot2::ylim(common_y_range)

  pl_x_range <- ggplot_build(x_dens)$layout$panel_params[[1]]$x.range
  x_range_pl <- ggplot_build(pl)$layout$panel_params[[1]]$x.range

  x_min <- min(pl_x_range[1], x_range_pl[1])
  x_max <- max(pl_x_range[2], x_range_pl[2])
  common_x_range <- c(x_min, x_max)

  x_dens <- x_dens + ggplot2::ylim(common_x_range)

  if(!is.null(additional_group)){
    group_levels <- levels(phyloseq::sample_data(ps_obj)[[additional_group]])
    comparisons <- combn(group_levels, 2, simplify = FALSE)

    # Axis.1 statistics and plot
    stats_axis1 <- ggpubr::compare_means(
      as.formula(paste("Axis.1 ~", additional_group)),
      data = P,
      method = "wilcox.test",
      p.adjust.method = "BH",
      comparisons = comparisons
    ) %>%
      rstatix::add_significance("p.adj") %>%
      rstatix::add_y_position(data = P,
                              formula = as.formula(paste("Axis.1 ~", additional_group)),
                              step.increase = 0)
    stats_axis1$xmin <- as.numeric(factor(stats_axis1$group1, levels = group_levels))
    stats_axis1$xmax <- as.numeric(factor(stats_axis1$group2, levels = group_levels))

    x_dens_add <- ggplot(P, aes_string(x = additional_group, y = "Axis.1", color = additional_group)) +
      #geom_violin(trim = FALSE, alpha = 0.1) +
      geom_boxplot(width = 0.5, alpha = 0.75) +
      geom_jitter(size = 1.5, alpha = 0.5, width = 0.1) +
      scale_color_manual(values = additional_palette) +
      ggpubr::stat_pvalue_manual(
        stats_axis1, label = "p.adj.signif",
        tip.length = 0.02, step.increase = 0.05,coord.flip = TRUE
      ) +
      coord_flip() +
      theme_bw(base_size = size) +
      theme(legend.position = "none")

    stats_axis2 <- ggpubr::compare_means(
      as.formula(paste("Axis.2 ~", additional_group)),
      data = P,
      method = "wilcox.test",
      p.adjust.method = "BH",
      comparisons = comparisons
    ) %>%
      rstatix::add_significance("p.adj") %>%
      rstatix::add_y_position(data = P,
                              formula = as.formula(paste("Axis.2 ~", additional_group)),
                              step.increase = 0.05)


    stats_axis2$xmin <- as.numeric(factor(stats_axis2$group1, levels = group_levels))
    stats_axis2$xmax <- as.numeric(factor(stats_axis2$group2, levels = group_levels))
    print(stats_axis2)

    y_dens_add <- ggplot(P, aes_string(x = additional_group, y = "Axis.2", color = additional_group)) +
      #geom_violin(trim = FALSE, alpha = 0.1) +
      geom_boxplot(width = 0.5, alpha = 0.75) +
      geom_jitter(size = 1.5, alpha = 0.5, width = 0.1) +
      scale_color_manual(values = additional_palette) +
      ggpubr::stat_pvalue_manual(
        stats_axis2, label = "p.adj.signif",
        tip.length = 0.02, step.increase = 0.05
      ) +
      theme_bw(base_size = size) +
      theme(legend.position = "right",
            legend.background = element_rect(color = "black"),
            legend.direction = "horizontal",
            axis.text.x = element_text(angle = 45))
    legend_add <- cowplot::get_legend(y_dens_add)
    y_dens_add <- y_dens_add + theme(legend.position = "none")

    pl_y_range <- ggplot_build(y_dens)$layout$panel_params[[1]]$y.range
    pl_y_range_add <- ggplot_build(y_dens_add)$layout$panel_params[[1]]$y.range
    y_range_pl <- ggplot_build(pl)$layout$panel_params[[1]]$y.range

    y_min <- min(pl_y_range[1], pl_y_range_add[1], y_range_pl[1])
    y_max <- max(pl_y_range[2], pl_y_range_add[2], y_range_pl[2])
    common_y_range <- c(y_min, y_max)

    y_dens <- y_dens + ggplot2::ylim(common_y_range)
    y_dens_add <- y_dens_add + ggplot2::ylim(common_y_range)
    y_dens <- ggarrange(y_dens_add,y_dens,ncol=2,nrow=1,align = "hv")

    pl_x_range <- ggplot_build(x_dens)$layout$panel_params[[1]]$x.range
    x_dens_x_range <- ggplot_build(x_dens_add)$layout$panel_params[[1]]$x.range
    x_range_pl <- ggplot_build(pl)$layout$panel_params[[1]]$x.range

    x_min <- min(pl_x_range[1], x_dens_x_range[1], x_range_pl[1])
    x_max <- max(pl_x_range[2], x_dens_x_range[2], x_range_pl[2])
    common_x_range <- c(x_min, x_max)

    x_dens <- x_dens + ggplot2::ylim(common_x_range)
    x_dens_add <- x_dens_add + ggplot2::ylim(common_x_range)
    x_dens <- ggarrange(x_dens,x_dens_add,ncol=1,nrow=2,align = "hv")

  }
  if(!is.null(additional_group)){
    legend_stacked <-
      cowplot::plot_grid(
        cowplot::ggdraw(legend$grobs[[1]]),NULL,
        cowplot::ggdraw(legend$grobs[[2]]),NULL,
        cowplot::ggdraw(legend_add$grobs[[1]]),
        ncol = 1,
        rel_heights = c(1,-0.5, 1,-0.5,1),
        align = "hv",
        axis = "tb",scale =1
      )
  } else {
    legend_stacked <-
      cowplot::plot_grid(
        cowplot::ggdraw(legend$grobs[[1]]),
        ncol = 1,
        align = "hv",
        axis = "tb",scale =1)
  }

  #pl <- pl +ggplot2::xlim(common_x_range) + ggplot2::ylim(common_y_range)

  # Arrange plots
  final_plot <- gridExtra::grid.arrange(
    y_dens, pl, legend_stacked, x_dens,
    ncol = 2, nrow = 2,
    widths = ratio, heights = rev(ratio)
  )

  return(final_plot)
}

