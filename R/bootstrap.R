#' Bootstrap OTU Table
#'
#' Perform bootstrap resampling on an OTU table with options for sampling with or without replacement.
#'
#' @param otu_table A matrix or data.frame with samples as columns and OTUs as rows
#' @param sample_size Integer specifying the size of the bootstrap sample (default: NULL = same as original)
#' @param with_replacement Logical indicating whether to sample with replacement (TRUE) or without (FALSE)
#' @param random_seed Integer for reproducible random sampling (default: NULL)
#'
#' @return A bootstrapped OTU table with the same structure as input
#'
#' @examples
#' # Create example OTU table
#' otu_data <- data.frame(
#'   Sample1 = c(10, 5, 3, 0),
#'   Sample2 = c(8, 2, 6, 1),
#'   Sample3 = c(15, 3, 0, 2),
#'   Sample4 = c(12, 4, 2, 0),
#'   Sample5 = c(9, 1, 5, 3),
#'   row.names = c("OTU1", "OTU2", "OTU3", "OTU4")
#' )
#'
#' # Bootstrap with replacement
#' boot_with_repl <- bootstrap_otu_table(otu_data, sample_size = 5, random_seed = 42)
#'
#' # Bootstrap without replacement
#' boot_without_repl <- bootstrap_otu_table(otu_data, sample_size = 3,
#'
#'                                        with_replacement = FALSE, random_seed = 42)
#'@export
bootstrap_otu_table <- function(otu_table, sample_size = NULL, with_replacement = TRUE, random_seed = NULL) {

  # Set random seed if provided
  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }

  # Get number of samples (columns)
  n_samples <- ncol(otu_table)

  # Set sample size if not specified
  if (is.null(sample_size)) {
    sample_size <- n_samples
  }

  # Check for valid sample size when sampling without replacement
  if (!with_replacement && sample_size > n_samples) {
    stop("Sample size cannot be larger than number of samples when with_replacement = FALSE")
  }

  # Generate bootstrap indices
  if (with_replacement) {
    # Sample with replacement
    boot_indices <- sample(n_samples, size = sample_size, replace = TRUE)
  } else {
    # Sample without replacement (permutation)
    boot_indices <- sample(n_samples, size = sample_size, replace = FALSE)
  }

  # Create bootstrapped table (preserve row names)
  boot_table <- otu_table[, boot_indices, drop = FALSE]

  # Copy row names if they exist
  if (!is.null(rownames(otu_table))) {
    rownames(boot_table) <- rownames(otu_table)
  }

  # Copy column names if they exist
  if (!is.null(colnames(otu_table))) {
    # For with_replacement = TRUE, we might have duplicates
    if (with_replacement && any(duplicated(boot_indices))) {
      # Append _repN to duplicate names
      col_counts <- table(boot_indices)
      new_names <- character(sample_size)
      counter <- rep(1, n_samples)

      for (i in seq_along(boot_indices)) {
        idx <- boot_indices[i]
        if (col_counts[as.character(idx)] > 1) {
          new_names[i] <- paste0(colnames(otu_table)[idx], "_rep", counter[idx])
          counter[idx] <- counter[idx] + 1
        } else {
          new_names[i] <- colnames(otu_table)[idx]
        }
      }
      colnames(boot_table) <- new_names
    } else {
      colnames(boot_table) <- colnames(otu_table)[boot_indices]
    }
  }

  return(boot_table)
}
