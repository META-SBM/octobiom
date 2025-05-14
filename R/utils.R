#' Get color palette by name
#'
#' @param name Palette name (e.g., "set_2_1")
#' @return Vector of HEX color codes
#' @export
#'
#' @examples
#' get_palette("set_2_1")
get_palette <- function(name) {
  if (!name %in% names(my_palettes)) {
    stop("Palette '", name, "' not found. Available palettes: ",
         paste(names(my_palettes), collapse = ", "))
  }
  my_palettes[[name]]
}

#' Display color palette
#'
#' @param name Palette name
#' @export
#'
#' @examples
#' show_palette("set_10_1")
show_palette <- function(name) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required to display palettes")
  }

  pal <- get_palette(name)
  df <- data.frame(
    colors = factor(pal, levels = pal),
    y = 1
  )

  ggplot2::ggplot(df, ggplot2::aes(x = colors, y = y, fill = colors)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste("Palette:", name)) +
    ggplot2::theme(legend.position = "none",
                   axis.text.y = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank())
}
