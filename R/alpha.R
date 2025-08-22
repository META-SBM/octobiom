#' @title Plot Alpha Diversity (Basic)
#' @description Internal function for basic alpha diversity plot
#' @param ps_obj A phyloseq object containing microbiome data
#' @param col Character string specifying the metadata column to group by
#' @param measure Alpha diversity measure to plot (default: 'Shannon')
plot_alpha_div <- function(ps_obj, group, color, measure
) {
  # Generate boxplot of alpha diversity without transformation
  p <- phyloseq::plot_richness(ps_obj, x = group, color = color, measures = measure)
  p$layers <- p$layers[-1]
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
#' @param colors Vector of colors for groups
#' @param size Base font size for plot elements (default: 10)
#' @param ncol Number of columns in the combined plot (default: length(measures)).
#' @param text_size size of text in figure
#' @param bracket_nudge_y Vertical adjustment to nudge brackets (default 0)
#'
#' @return A list containing:
#' \itemize{
#'   \item `plots` - Combined ggplot object of all measures.
#'   \item `stats` - List of statistical results for each measure.
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' colors <- c("Group1" = "blue", "Group2" = "red")
#' result <- create_alpha_plots(
#'   ps_obj = my_phyloseq,
#'   col = "Treatment",
#'   measures = c("Shannon", "Simpson"),
#'   colors = colors,
#'   bracket_nudge_y = 0.1
#' )
#' print(p)
#' }
#'
#' @importFrom ggplot2 ggplot geom_violin geom_boxplot geom_jitter position_jitterdodge
#' @importFrom ggplot2 scale_color_manual theme_bw theme element_text position_dodge
#' @importFrom rstatix add_significance add_x_position add_y_position
#' @importFrom patchwork wrap_plots
#'
#' @export
create_alpha_plots <- function(ps_obj, col,
                               measures = c("Shannon"),
                               method = "wilcox.test",
                               colors, size = 10,ncol = NULL, text_size = 5,bracket_nudge_y = 0,
                               font_family = "Times New Roman") {
  if (!font_family %in% extrafont::fonts()) {
    warning(paste("Font", font_family, "not found. Using default sans font."))
    font_family <- "sans"  # Fallback to default
  }

  plot_list <- list()
  stat_list <- list()
  # Automatically get factor levels from metadata
  group_factor <- phyloseq::sample_data(ps_obj)[[col]]
  group_levels <- sort(unique(as.character(group_factor)))

  # Convert to factor with natural ordering
  phyloseq::sample_data(ps_obj)[[col]] <- factor(
    group_factor,
    levels = group_levels
  )
  sample_sizes <- table(sample_data(ps_obj)[[col]])

  for (measure in measures) {
    # Create base plot
    p <- plot_alpha_div(ps_obj, group = col, color = col, measure = measure) +
      ggplot2::geom_blank() +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.1, aes_string(fill = col)) +
      ggplot2::geom_boxplot(width = 0.2, alpha = 0.75, outlier.shape = NA) +
      ggplot2::geom_jitter(size = 0.5, alpha = 0.5, width = 0.1, aes_string(color = col)) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::labs(title = paste(measure),y='') +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = element_text(angle = 45, hjust = 1, family = font_family, size = size),
        axis.text.y = element_text(family = font_family, size = size),
        axis.title = element_text(family = font_family, size = size, face = "bold"),
        legend.text = element_text(family = font_family, size = size),
        legend.title = element_text(family = font_family, size = size, face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5, family = font_family),
        strip.text = element_blank(),
        text = element_text(family = font_family, size = size)
      )+
      ggplot2::scale_x_discrete(labels = function(x) paste0(x, "\n(n=", sample_sizes[x], ")"))
    # Prepare data
    df <- as.data.frame(p$data)
    # Statistical comparisons with automatic group detection
    comparisons <- combn(group_levels, 2, simplify = FALSE)

    stat_res <- ggpubr::compare_means(
      as.formula(paste("value ~", col)),
      data = df,
      method = method,
      p.adjust.method = "BH",
      comparisons = comparisons
    ) %>%
      rstatix::add_significance("p.adj") %>%
      rstatix::add_x_position() %>%
      rstatix::add_y_position(
        data = df,
        formula = as.formula(paste("value ~", col)),
        step.increase = 0.2
      )
    # Auto-adjust x positions based on factor levels
    #stat_res$xmin <- match(stat_res$group1, group_levels)
    #stat_res$xmax <- match(stat_res$group2, group_levels)

    stat_res$p.adj.formatted <- base::format.pval(
      stat_res$p.adj,
      digits = 2,
      eps = 0.001
    )


    stat_res$label <- paste0(
      stat_res$p.adj.signif
    )
    stat_res$label1 <-paste0(
     '(p=', stat_res$p.adj.formatted, ')'
    )

    # Добавляем на график
    p <- p +
      ggpubr::stat_pvalue_manual(
        stat_res,
        label = "label",
        tip.length = 0.01,
        step.increase = 0.15,
        size = text_size,
        bracket.size = 0.1,
        vjust = 0.5,
        parse = TRUE,
        bracket.nudge.y = bracket_nudge_y
      )+
      ggpubr::stat_pvalue_manual(
        stat_res,
        label = "label1",
        tip.length = 0.01,
        step.increase = 0.15,
        size = text_size,
        bracket.size = 0.1,
        vjust = 2.3,
        parse = TRUE,
        bracket.nudge.y = bracket_nudge_y
      )+
      theme(text = element_text(family = font_family, lineheight = 0.5))

    # Save results
    plot_list[[measure]] <- p
    stat_list[[measure]] <- stat_res
  }


  if (is.null(ncol)) ncol <- length(measures)
  combined_plot <- patchwork::wrap_plots(plot_list, ncol = ncol) &
    theme(text = element_text(family = font_family))

  # Возвращаем результаты
  return(list(
    plots = combined_plot,
    stats = stat_list
  ))
}


