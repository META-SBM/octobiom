#' Prepare Phyloseq Object with Transformations and Core Taxa Filtering
#'
#' Applies taxonomic aggregation, transformation
#'  and core taxa filtering to a phyloseq object.
#'
#' @param phyloseq A phyloseq object to process.
#' @param transformation Transformation method:
#' 'compositional' (relative abundance), 'clr' (centered log-ratio) or NONE.
#' @param detection Abundance threshold for considering
#' taxa presence
#' @param prevalence Prevalence threshold
#' (proportion of samples where taxa must be present; 0-1).
#' @param taxonomic_level Taxonomic level for aggregation (e.g., "Genus").
#' @param first Parameter indicating that the data
#' is filtered or transformed the first.
#' @param filter_column Column name for sample filtering (optional)
#' @param filter_value Value to filter by in filter_column (optional)
#' @param filter_unclassified is flag to filter Kingdom == UNCLASSIFIED from MetaPhlAn
#' @return A processed phyloseq object with core taxa meeting criteria.
#' @examples
#' # ps <- prepare_table(phyloseq_obj, transformation = "compositional",
#' #                    detection = 0.001, prevalence = 0.1,
#' # taxonomic_level = "Genus", first = 'filter')
#' @export
prepare_phyloseq <- function(phyloseq,
                          transformation = "compositional",
                          detection = 0,
                          prevalence = 0.1,
                          taxonomic_level = "Species",
                          first = "transform",
                          filter_column = NULL,
                          filter_value = NULL,
                          filter_unclassified = T) {
  # Input validation
  if (!inherits(phyloseq, "phyloseq")) {
    stop("Input must be a phyloseq object.")
  }
  if (!transformation %in% c("compositional", "clr", "NONE")) {
    stop("Transformation must be 'compositional', 'clr' or NONE.")
  }
  if ((detection < 0) & (transformation == "compositional")) {
    stop("Detection threshold cannot be negative in compositional data.")
  }
  if (prevalence < 0 || prevalence > 1) {
    stop("Prevalence must be between 0 and 1.")
  }
  # Sample filtering if requested
  if (!is.null(filter_column)) {
    if (!filter_column %in% sample_variables(phyloseq)) {
      stop("Filter column '", filter_column, "' not found in sample data.")
    }
    if (is.null(filter_value)) {
      stop("Filter value must be specified when filter_column is provided.")
    }
    # Subset phyloseq for this time point
    metadata <- as(phyloseq@sam_data, "data.frame")
    samples_tp <- rownames(metadata[metadata[[filter_column]] == filter_value,])
    # Subsample phyloseq
    phyloseq <- phyloseq::prune_samples(samples_tp, phyloseq)
    if (phyloseq::nsamples(phyloseq) == 0) {
      stop("No samples remaining after filtering for ",
           filter_column, " == ", filter_value)
    }
    message("Filtered samples by ", filter_column, " == ", filter_value,
            ". Remaining samples: ", phyloseq::nsamples(phyloseq))
  }
  # Taxonomic aggregation
  phy_agg <- speedyseq::tax_glom(
    physeq = phyloseq,
    taxrank = taxonomic_level,
    NArm = TRUE
  )
  if (filter_unclassified == T){
    phy_agg <- phyloseq::subset_taxa(phy_agg,Kingdom != "UNCLASSIFIED")
  }
  if ((first == "transform") & (transformation != 'NONE')) {
    phy_trans <- microbiome::transform(phy_agg, transformation)

    core_taxa <- tryCatch(
      expr = {
        microbiome::core_members(
          phy_trans,
          detection = detection,
          prevalence = prevalence
        )
      },
      error = function(e) {
        message("No taxa met the detection/prevalence criteria.")
        return(character(0))  # Return empty vector instead of error
      }
    )
    # Subset phyloseq object to core taxa
    if (length(core_taxa) > 0) {
      phy_filtered <- phyloseq::prune_taxa(core_taxa, phy_trans)
    } else {
      warning("Returning empty phyloseq object - no core taxa found.")
      phy_filtered <- phyloseq::prune_taxa(character(0), phy_trans)
    }
    print(paste("Number of taxa after filtering:",
                phyloseq::ntaxa(phy_filtered)))
    return(phy_filtered)
  } else if ((first == "filter") & (transformation != 'NONE')) {
    # Core taxa identification
    core_taxa <- tryCatch(
      expr = {
        microbiome::core_members(
          phy_agg,
          detection = detection,
          prevalence = prevalence
        )
      },
      error = function(e) {
        message("No taxa met the detection/prevalence criteria.")
        return(character(0))  # Return empty vector instead of error
      }
    )
    # Subset phyloseq object to core taxa
    if (length(core_taxa) > 0) {
      phy_filtered <- phyloseq::prune_taxa(core_taxa, phy_agg)
    } else {
      warning("Returning empty phyloseq object - no core taxa found.")
      phy_filtered <- phyloseq::prune_taxa(character(0), phy_agg)
    }
    print(paste("Number of taxa after filtering:",
                phyloseq::ntaxa(phy_filtered)))

    phy_trans <- microbiome::transform(phy_filtered, transformation)

    return(phy_trans)
  } else {
    core_taxa <- tryCatch(
      expr = {
        microbiome::core_members(
          phy_agg,
          detection = detection,
          prevalence = prevalence
        )
      },
      error = function(e) {
        message("No taxa met the detection/prevalence criteria.")
        return(character(0))  # Return empty vector instead of error
      }
    )
    # Subset phyloseq object to core taxa
    if (length(core_taxa) > 0) {
      phy_filtered <- phyloseq::prune_taxa(core_taxa, phy_agg)
    } else {
      warning("Returning empty phyloseq object - no core taxa found.")
      phy_filtered <- phyloseq::prune_taxa(character(0), phy_agg)
    }
    print(paste("Number of taxa after filtering:",
                phyloseq::ntaxa(phy_filtered)))

    return(phy_filtered)
  }
}
