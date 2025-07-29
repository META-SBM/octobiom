#' Generate customized microbiome heatmap with taxonomic and sample annotations
#'
#' @param ps_obj Phyloseq object containing microbiome data
#' @param tax_level_annotation Taxonomic level for row annotation (e.g., "Phylum"). Default: NULL.
#' @param features Vector of sample metadata features for column annotation. Default: NULL.
#' @param tax_colors List of color palettes for taxonomic levels (named by level)
#' @param feature_colors List of color palettes for sample features (named by feature)
#' @param row_name_size Numeric parameter for the font size of line names
#' @param col_name_size Numeric parameter for the font size of column names

#' @return Invisibly returns the ComplexHeatmap object
#' @export
#'
#' @examples
#' \dontrun{
#' data(GlobalPatterns)
#'
#' # Define color palettes
#' phylum_colors <- c(
#'   "Proteobacteria" = "#1B9E77",
#'   "Actinobacteria" = "#D95F02",
#'   "Bacteroidetes" = "#7570B3"
#' )
#'
#' sampletype_colors <- c(
#'   "Feces" = "#E41A1C",
#'   "Skin" = "#377EB8",
#'   "Ocean" = "#4DAF4A"
#' )
#'
#' generate_microbiome_heatmap(
#'   ps_obj = GlobalPatterns,
#'   tax_level_annotation = "Phylum",
#'   features = "SampleType",
#'   tax_colors = list(Phylum = phylum_colors),
#'   feature_colors = list(SampleType = sampletype_colors)
#' )
#' }
generate_microbiome_heatmap <- function(ps_obj,
                                       tax_level_annotation,
                                       features,
                                       tax_colors,
                                       feature_colors,
                                       row_name_size = 0,
                                       col_name_size = 0) {

  # Input validation
  if (!inherits(ps_obj, "phyloseq")) {
    stop("Input must be a phyloseq object")
  }

  # Extract data from phyloseq object
  otu <- as(phyloseq::otu_table(ps_obj), "matrix")
  tax <- as(phyloseq::tax_table(ps_obj), "data.frame")
  meta <- as(phyloseq::sample_data(ps_obj), "data.frame")

  # Row annotation (taxonomic)
  row_annot <- NULL
  if (!is.null(tax_level_annotation)) {
    if (!tax_level_annotation %in% colnames(tax)) {
      stop("Taxonomic level not found in tax table: ", tax_level_annotation)
    }

    tax_vec <- tax[[tax_level_annotation]]
    tax_unique <- sort(unique(na.omit(tax_vec)))  # Удаляем NA значения

    # Получаем цвета для таксономии
    if (!is.null(tax_colors) && !is.null(tax_colors[[tax_level_annotation]])) {
      tax_palette <- tax_colors[[tax_level_annotation]]
      # Проверяем, все ли таксоны имеют цвета
      missing_taxa <- setdiff(tax_unique, names(tax_palette))
      if (length(missing_taxa) > 0) {
        warning("Colors missing for some taxa in ", tax_level_annotation,
                ": ", paste(missing_taxa, collapse = ", "))
        default_cols <- scales::hue_pal()(length(missing_taxa))
        names(default_cols) <- missing_taxa
        tax_palette <- c(tax_palette, default_cols)
      }
    } else {
      message("Using default palette for ", tax_level_annotation)
      tax_palette <- scales::hue_pal()(length(tax_unique))
      names(tax_palette) <- tax_unique
    }

    row_annot <- ComplexHeatmap::rowAnnotation(
      df = data.frame(Tax = tax_vec),
      col = list(Tax = tax_palette),
      annotation_legend_param = list(
        Tax = list(title = tax_level_annotation)
      ),
      show_annotation_name = FALSE
    )
  }

  # Top annotation (sample features)
  top_annot <- NULL
  if (!is.null(features)) {
    if (!all(features %in% colnames(meta))) {
      stop("Features not found in sample_data: ",
           paste(setdiff(features, colnames(meta)), collapse = ", "))
    }

    ann_df <- meta[, features, drop = FALSE]
    ann_col <- list()

    for (f in features) {
      uniq_vals <- sort(unique(na.omit(ann_df[[f]])))  # Удаляем NA значения

      # Получаем цвета для признаков
      if (!is.null(feature_colors) && !is.null(feature_colors[[f]])) {
        feature_palette <- feature_colors[[f]]
        # Проверяем, все ли уровни имеют цвета
        missing_levels <- setdiff(uniq_vals, names(feature_palette))
        if (length(missing_levels) > 0) {
          warning("Colors missing for some levels in feature ", f,
                  ": ", paste(missing_levels, collapse = ", "))
          default_cols <- scales::hue_pal()(length(missing_levels))
          names(default_cols) <- missing_levels
          feature_palette <- c(feature_palette, default_cols)
        }
      } else {
        message("Using default palette for feature: ", f)
        feature_palette <- scales::hue_pal()(length(uniq_vals))
names(feature_palette) <- uniq_vals
      }

      ann_col[[f]] <- feature_palette
    }

    top_annot <- ComplexHeatmap::HeatmapAnnotation(
      df = ann_df,
      col = ann_col,
      show_annotation_name = TRUE
    )
  }

  hm <- ComplexHeatmap::Heatmap(
    matrix = t(otu),
    col = circlize::colorRamp2(
      c(min(otu), 0, max(otu)),
      c("#45CAFF", "white", "#FF1B6B")
    ),
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D2",
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "ward.D2",
    row_names_gp = grid::gpar(fontsize = row_name_size),
    column_names_gp = grid::gpar(fontsize = col_name_size),
    #show_row_names = FALSE,
    #show_column_names = FALSE,
    left_annotation = row_annot,
    top_annotation = top_annot,
    heatmap_legend_param = list(
      title = "Abundance",
      direction = "horizontal",
      legend_width = unit(4, "cm"),
      labels_gp = grid::gpar(fontsize = 8),
      title_gp = grid::gpar(fontsize = 8)
    )
  )
  heatmap_grob <- grid::grid.grabExpr({
    ComplexHeatmap::draw(hm,heatmap_legend_side = "right",
         annotation_legend_side = "right",merge_legend = T)
  })
  heatmap_plot <- cowplot::plot_grid(heatmap_grob)
  return(heatmap_plot)

  }

