#' @title Plot Alpha Diversity (Basic)
#' @description Internal function for basic alpha diversity plot
#' @param ps_obj A phyloseq object containing microbiome data
#' @param col Character string specifying the metadata column to group by
#' @param measure Alpha diversity measure to plot (default: 'Shannon')
plot_alpha_div <- function(ps_obj, group, color, measure
) {
  # Generate boxplot of alpha diversity without transformation
  p <- phyloseq::plot_richness(ps_obj, x = group, color = color, measures = measure)
  return(p)
}
#' @title Create Alpha Diversity Plots with Statistical Comparisons
#'
#' @description This function generates alpha diversity plots (e.g., Shannon index) from a phyloseq object
#' with optional statistical comparisons between groups. It creates violin-boxplot combinations
#' with jittered points and adds significance annotations.
#'
#' @param ps_obj A phyloseq object containing microbiome data
#' @param col Character string specifying the metadata column to group by
#' @param measure Alpha diversity measure to plot (default: 'Shannon')
#' @param method Statistical test method (default: 'wilcox.test')
#' @param color Vector of colors for groups
#' @param size Base font size for plot elements (default: 10)
#'
#' @return A ggplot object showing alpha diversity with statistical annotations
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' colors <- c("Group1" = "blue", "Group2" = "red")
#' p <- create_alpha_plots(ps_obj = my_phyloseq,
#'                        col = "Treatment",
#'                        measure = "Shannon",
#'                        color = colors,
#'                        level_factors = c("Control", "Treatment"))
#' print(p)
#' }
#'
#' @importFrom ggplot2 ggplot geom_violin geom_boxplot geom_jitter position_jitterdodge
#' @importFrom ggplot2 scale_color_manual theme_bw theme element_text position_dodge
#' @importFrom rstatix add_significance add_x_position add_y_position
#'
#' @export
create_alpha_plots <- function(ps_obj, col, measure = "Shannon",
                               method = "wilcox.test",
                               colors, size = 10) {

  # Automatically get factor levels from metadata
  group_factor <- phyloseq::sample_data(ps_obj)[[col]]
  group_levels <- sort(unique(as.character(group_factor)))

  # Convert to factor with natural ordering
  phyloseq::sample_data(ps_obj)[[col]] <- factor(
    group_factor,
    levels = group_levels
  )

  # Create base plot
  p <- plot_alpha_div(ps_obj, group = col, color = col, measure = measure) +
    geom_violin(trim = FALSE, alpha = 0.1) +
    geom_boxplot(width = 0.5, alpha = 0.75, position = position_dodge(0.9)) +
    geom_jitter(size = 1.5, alpha = 0.5,
                position = position_jitterdodge(jitter.width = 0.1,
                                                dodge.width = 0.9)) +
    scale_color_manual(values = colors) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  # Prepare data
  df <- as.data.frame(p$data)
  group <- col

  # Statistical comparisons with automatic group detection
  comparisons <- combn(group_levels, 2, simplify = FALSE)

  anno_df <- ggpubr::compare_means(
    as.formula(paste("value ~", group)),
    data = df,
    method = method,
    p.adjust.method = "BH",
    comparisons = comparisons
  ) %>%
    add_significance("p.adj") %>%
    add_x_position() %>%
    add_y_position(
      data = df,
      formula = as.formula(paste("value ~", group)),
      step.increase = 0.2
    )

  # Auto-adjust x positions based on factor levels
  anno_df$xmin <- match(anno_df$group1, group_levels)
  anno_df$xmax <- match(anno_df$group2, group_levels)

  # Add annotations
  p <- p + ggpubr::stat_pvalue_manual(
    anno_df,
    label = "p.adj.signif",
    tip.length = 0.02,
    step.increase = 0.05
  )

  return(list(plot = p, stats = anno_df))
}


