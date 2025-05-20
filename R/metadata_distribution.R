#' Plot Distributions of Metadata Variables in a Phyloseq Object
#'
#' Creates distribution plots for specified metadata variables, handling both
#' categorical and numeric variables appropriately.
#'
#' @param physeq A phyloseq object
#' @param features Character vector of metadata variables to plot
#' @param color_palette Optional color palette for categorical variables
#' @param text_size Base text size (default: 12)
#' @param ncol Number of columns for combined plot (NULL for automatic)
#' @param nrow Number of rows for combined plot (NULL for automatic)
#' @param combined_plot Return combined plot (TRUE) or list of plots (FALSE)
#' @return Either a combined plot or list of ggplot objects
#' @importFrom rlang sym
#' @export
plot_metadata_distributions <- function(physeq,
                                        features,
                                        color_palette = NULL,
                                        text_size = 12,
                                        ncol = NULL,
                                        nrow = NULL,
                                        combined_plot = TRUE) {

  # Validate input features
  if (!all(features %in% colnames(phyloseq::sample_data(physeq)))) {
    missing_features <- setdiff(features, colnames(phyloseq::sample_data(physeq)))
    stop(paste("The following features are missing from metadata:",
               paste(missing_features, collapse = ", ")))
  }

  # Extract and prepare metadata
  meta_df <- phyloseq::sample_data(physeq) %>%
    tibble::as_tibble() %>%
    dplyr::select(all_of(features)) %>%
    dplyr::mutate(across(where(is.character), as.factor))

  # Create default color palette if not provided
  if (is.null(color_palette)) {
    factor_levels <- meta_df %>%
      select(where(is.factor)) %>%        # select only factor columns
      map(levels) %>%                     # get levels of each factor column
      flatten_chr()                      # flatten to a single character vector

    # Count unique levels
    n_colors <- length(unique(factor_levels))

    # Create color palette
    color_palette <- scales::hue_pal()(n_colors)
  }

  # Internal plotting function for each variable
  plot_variable <- function(var_name, data) {
    var <- data[[var_name]]

    if (is.factor(var) || is.character(var)) {
      # Handle categorical variables
      if (is.character(var)) var <- as.factor(var)

      levels_count <- length(levels(var))

      # Prepare plot data
      plot_data <- data %>%
        dplyr::count(!!sym(var_name))

      # Create plot
      p <- ggplot2::ggplot(plot_data,
                           ggplot2::aes(x = !!sym(var_name), y = n, fill = !!sym(var_name))) +
        ggplot2::geom_col() +
        ggplot2::scale_fill_manual(values = color_palette) +
        ggplot2::labs(x = var_name, y = "Count") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          text = ggplot2::element_text(size = text_size),
          axis.text.x = ggplot2::element_text(angle = ifelse(levels_count > 5, 45, 0), hjust = 1),
          legend.position = "none"
        )

      if (levels_count > 5) {
        p <- p + ggplot2::coord_flip()
      }

    } else if (is.numeric(var)) {
      # Handle numeric variables
      p <- ggplot2::ggplot(data, ggplot2::aes(x = !!sym(var_name))) +
        ggplot2::geom_histogram(bins = 30, fill = color_palette[1], color = "white") +
        ggplot2::labs(x = var_name, y = "Count") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          text = ggplot2::element_text(size = text_size)
        )
    } else {
      return(NULL)
    }

    return(p)
  }

  # Create plots for all variables
  plots <- purrr::map(features, ~ plot_variable(., meta_df)) %>%
    purrr::set_names(features) %>%
    purrr::compact()

  if (length(plots) == 0) {
    message("No plottable variables found in metadata")
    return(invisible(NULL))
  }

  # Return either combined plot or list of plots
  if (combined_plot) {
    if (is.null(ncol)) ncol <- min(2, length(plots))
    if (is.null(nrow)) nrow <- ceiling(length(plots) / ncol)
    return(patchwork::wrap_plots(plots, ncol = ncol, nrow = nrow))
  } else {
    return(plots)
  }
}
