#' Color palettes for MyPackage
#'
#' A collection of curated color palettes for data visualization.
#'
#' @format A named list of color vectors in HEX format.
#' @export
#'
#' @examples
#' # View all palette names
#' names(my_palettes)
#'
#' # Access a specific palette
#' my_palettes$set_2_1
my_palettes <- list(
  set_2_1 = c("#FF1B6B", "#45CAFF"),
  set_2_2 = c("#D5B3FF", "#FBD960"),
  set_2_3 = c("#08415C", "#D62828"),
  set_2_4 = c("#696EFF", "#F8ACFF"),

  set_10_1 = c("#001219", "#005f73", "#0a9396", "#94d2bd", "#e9d8a6",
               "#ee9b00", "#ca6702", "#bb3e03", "#ae2012", "#9b2226"),

  set_10_2 = c("#F55B66", "#ACACAC", "#5CD5DA", "#0089DF", "#19427D",
               "#FF8300", "#FFF53B", "#E7AD99", "#EA7A67", "#7D4F50"),

  set_10_3 = c("#CDB4DB", "#FFC0D8", "#FF76A8", "#8ECBFF", "#647F98",
               "#FFA69E", "#FFE699", "#80F4DC", "#6ACFE0", "#4A4E59"),

  set_10_4 = c("#FF8400", "#3C096C", "#FF0A54", "#2D00F7", "#029191",
               "#08415C", "#CC2936", "#7E82F0", "#F8ACFF", "#B400FF"),

  set_20_1 = c("#001219", "#005f73", "#0a9396", "#94d2bd", "#e9d8a6",
               "#ee9b00", "#ca6702", "#bb3e03", "#ae2012", "#9b2226",
               "#F55B66", "#ACACAC", "#5CD5DA", "#0089DF", "#19427D",
               "#FF8300", "#FFF53B", "#E7AD99", "#EA7A67", "#7D4F50"),

  set_20_2 = c("#CDB4DB", "#FFC0D8", "#FF76A8", "#8ECBFF", "#647F98",
               "#FFA69E", "#FFE699", "#80F4DC", "#6ACFE0", "#4A4E59",
               "#FF8400", "#3C096C", "#FF0A54", "#2D00F7", "#029191",
               "#08415C", "#CC2936", "#7E82F0", "#F8ACFF", "#B400FF")
)
# Save to data/ directory
usethis::use_data(my_palettes, overwrite = TRUE)
