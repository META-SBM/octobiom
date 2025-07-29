#' PERMANOVA Analysis and Visualization
#'
#' Performs PERMANOVA analysis and creates visualization plots showing the variance
#' explained by different metadata variables.
#'
#' @param ps_obj A phyloseq object containing OTU table and sample data
#' @param formula Character string specifying the PERMANOVA formula (e.g., "Group + Condition")
#' @param method Distance method ("bray", "euclidean", "unifrac", "wunifrac")
#' bray and wunifrac take into account the relative abundance, when euclidean
#' takes into account clr-data (aitchison distance)
#' @param show_plot Logical, whether to display the plot (default: TRUE)
#' @param size Base size for plot elements
#' @param cols_meta Dataframe mapping metadata columns to categories
#' @param palette Color palette for categories
#' @param threshold Significance threshold line position (default:0)
#'
#' @return A list containing:
#' \itemize{
#'   \item PERMANOVA results (adonis2 object)
#'   \item Combined ggplot object (main plot + circular plot)
#'   \item Processed results dataframe
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' result <- plot_permanova(
#'   ps_obj = ps,
#'   formula = "Group + Timepoint",
#'   method = "bray",
#'   cols_meta = metadata_mapping,
#'   palette = color_list
#' )
#' }
#'
#' @export
#' @importFrom phyloseq distance
#' @importFrom vegan adonis2
#' @importFrom ggplot2 ggplot geom_col labs theme_classic coord_flip geom_label
#' @importFrom ggplot2 theme element_rect element_blank element_text
#' @importFrom ggplot2 aes element_line facet_grid guides guide_legend geom_rect coord_polar xlim theme_void
#' @importFrom ggtext element_markdown
#' @importFrom ggplot2 scale_fill_manual scale_color_manual geom_hline
#' @importFrom patchwork plot_layout
#' @importFrom glue glue
#' @importFrom dplyr %>% mutate filter group_by arrange desc ungroup left_join
#' @importFrom stats as.formula
#' @importFrom utils head
plot_permanova <- function(ps_obj, formula,
                           method = "bray", show_plot = TRUE,
                           size = 10,
                           cols_meta,
                           palette, threshold = 0) {
  Feature <- R2 <- ymin <- ymax <- label <- labelPosition <- signif_disp <- category <- new_name <- NULL

  # Convert sample data to a data frame
  metadf <- as.data.frame(as.matrix(ps_obj@sam_data))

  if(method == 'euclidean'){
    distance <- 'aitchison'
    bd1 <- ps_obj %>%
      microViz::dist_calc(distance,dist= distance) %>%
      microViz::tax_transform(trans = "clr")%>%
      microViz::ord_calc(method = "PCA") %>%
      microViz::dist_bdisp(variables = cols_meta$new_name) %>%
      microViz::bdisp_get()
  }else{
    distance <- method
    bd1 <- ps_obj %>%
      microViz::dist_calc(distance,dist= distance) %>%
      microViz::dist_calc(distance,dist= distance) %>%
      microViz::dist_bdisp(variables = cols_meta$new_name) %>%
      microViz::bdisp_get()
  }

  # Calculate distance matrix based on the specified method
  dist <- phyloseq::distance(ps_obj, method = method)

  # Create formula from the string
  formula1 <- as.formula(paste("dist ~", formula))

  # Perform adonis analysis with the specified formula
  permanova_result <- adonis2(formula1, data = metadf,by='term',parallel = 30)

  # Process adonis results for plotting
  res <- as.data.frame(permanova_result)
  res$meta <- row.names(res)
  res <- subset(res, res$meta != 'Total')
  res$meta <- gsub('Residuals', 'other', res$meta)
  res[which(res$`Pr(>F)` > 0.05), 'signif'] <- ' '
  res[which(res$`Pr(>F)` <= 0.05), 'signif'] <- '*'
  res[which(res$`Pr(>F)` <= 0.01), 'signif'] <- '**'
  res[which(res$`Pr(>F)` <= 0.001), 'signif'] <- '***'
  res$Feature <- as.character(res$meta)
  res <- res %>%
    mutate(  Feature = ifelse((signif == '*') | (signif == '**') | (signif == '***') ,
                              glue("**{Feature}**"),  # Use glue to format as bold
                              Feature))
  res$Feature <- factor(res$Feature, levels = res$Feature[order(res$R2, decreasing = FALSE)])

  res$Feature <- stringr::str_wrap(res$Feature, width = 1)
  res <- tibble::rownames_to_column(res, var = "new_name")
  res <- res %>%
    left_join(cols_meta, by = 'new_name')
  res <- res %>%
    filter(new_name != "Residual")
  permdisp_df <- data.frame(new_name = character(), p_val = numeric(), stringsAsFactors = FALSE)

  for (elem in cols_meta$new_name) {
    # Извлекаем p-value из bd1 для текущего элемента
    p_value <- bd1[[elem]]$anova$`Pr(>F)`[1]  # Используем двойные квадратные скобки для доступа к элементу

    # Добавляем новую строку в results
    permdisp_df <- rbind(permdisp_df, data.frame(new_name = elem, p_val = p_value))
  }
  permdisp_df[which(permdisp_df$p_val > 0.05), 'signif_disp'] <- ' '
  permdisp_df[which(permdisp_df$p_val <= 0.05), 'signif_disp'] <- '*'
  permdisp_df[which(permdisp_df$p_val <= 0.01), 'signif_disp'] <- '**'
  permdisp_df[which(permdisp_df$p_val <= 0.001), 'signif_disp'] <- '***'

  res <- merge(res,permdisp_df,by.x='meta',by.y='new_name',all.x=T)
  #res <- res %>%
  #  dplyr::left_join(permdisp_df,by='new_name')
  #print(res[1,])
  # Create the ggplot object for plotting
  res <- res %>%
    group_by(category) %>%          # Group by 'category'
    arrange(desc(R2), .by_group = TRUE) %>%  # Sort within each group by 'R2' in descending order
    ungroup()
  res$Feature <- factor(res$Feature, levels = res$Feature[order(-res$R2)])
  #print(res)
  main_plot <- ggplot(res, aes(x = Feature, y = R2, fill = category, order =R2)) +
    geom_col() +  # Use geom_col() to plot actual values
    labs(title = "",
         x = "",
         y = "R2") +
    theme_classic() +  # Optional: use a minimal theme
    coord_flip()+
    geom_label(data = res %>% filter(signif != ' '),
               aes(label = signif,color='PERMANOVA \n significance'),  # Set fontface based on significance
               size = 3, hjust = -0.2,fill='white')+
    geom_label(data = res %>% filter(signif_disp != ' '),
               aes(label = signif_disp,color = 'PERMDISP \n significance'),  # Set fontface based on significance
               size = 3,fill='white', hjust = 1)+
    theme(
      panel.background = element_rect(fill="white"),
      legend.position = 'top',               # Position legend at top left corner
      legend.title = element_blank(),
      legend.text = element_text(size=12,color='black'),
      legend.direction = 'horizontal',
      legend.justification = "left",
      axis.title.y = element_blank(),
      axis.line.x = element_line(color="black"),
      axis.line.y = element_blank(),
      axis.text.y = ggtext::element_markdown(size=12, hjust=1, color='black', linewidth = 1.2),
      axis.text.x = element_text(size=12, hjust=1, color='black'),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      strip.text.x = element_blank(),  # Remove x strip labels
      strip.text.y = element_blank()) +
    facet_grid(category ~ ., scales = "free", space = "free")+
    guides(fill = guide_legend(nrow = 1))+
    scale_fill_manual(values = palette)+
    scale_color_manual(values=c("PERMANOVA \n significance"="red", "PERMDISP \n significance"="black"))+
    #scale_fill_manual(values = c("Significance"="white", "Dispersion"="black"))+
    geom_hline(yintercept = threshold, col = "grey30", lty = "dashed")

  total_effect_size <- round(sum(res$R2)*100,3)
  unexplained_variance <- 100 - total_effect_size

  # Prepare data for circular barplot
  circular_data <- data.frame(
    category = c("Explained Variance", "Unexplained Variance"),
    value = c(total_effect_size, unexplained_variance)
  )


  circular_data$fraction <- circular_data$value / sum(circular_data$value)

  # Compute the cumulative percentages (top of each rectangle)
  circular_data$ymax <- cumsum(circular_data$fraction)

  # Compute the bottom of each rectangle
  circular_data$ymin <- c(0, head(circular_data$ymax, n=-1))

  # Compute label position
  circular_data$labelPosition <- (circular_data$ymax + circular_data$ymin) / 2

  # Compute a good label
  circular_data$label <- paste0(circular_data$category, "\n value: ", circular_data$value)

  # Make the plot
  circular_plot <- ggplot(circular_data, aes(ymax=ymax, ymin=ymin, xmax=3.5, xmin=3, fill=category)) +
    geom_rect() +
    geom_label( x=4, aes(y=labelPosition, label=label), size=5) +
    scale_fill_manual(values = palette)+
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    theme_void() +
    theme(legend.position = "none")
  #plot(circular_plot)
  combined_plot <- main_plot + circular_plot + plot_layout(ncol=2,widths = c(2,1))

  # Check if plot should be generated
  if (show_plot) {
    print(combined_plot)
  }
  # Return both adonis results and the ggplot object
  return(list(permanova_result, combined_plot,res))
}


