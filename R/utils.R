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
#' @title Setup Fonts for Plotting
#' @description Initialize and load Times New Roman font for ggplot2
#' @param font_family Font family name to load (default: "Times New Roman")
#' @importFrom extrafont font_import loadfonts fonts
#' @export
setup_plot_fonts <- function(font_family = "Times New Roman") {
  # Check if font is already available
  if (!font_family %in% extrafont::fonts()) {
    message("Importing fonts...")
    extrafont::font_import(prompt = FALSE)
    extrafont::loadfonts(device = "win", quiet = TRUE)
  }

  if (font_family %in% extrafont::fonts()) {
    message(paste("Font", font_family, "successfully loaded"))
  } else {
    warning(paste("Font", font_family, "not found. Using default sans font."))
  }
}
