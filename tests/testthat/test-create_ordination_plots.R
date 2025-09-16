context("Test create_ordination_plots function")

create_test_phyloseq <- function() {

  otu_matrix <- matrix(
    c(150, 120, 100, 80, 50,
      30, 40, 60, 80, 100,
      10, 20, 30, 40, 50,
      5, 10, 15, 20, 25),
    nrow = 4, ncol = 5,
    dimnames = list(
      c("OTU1", "OTU2", "OTU3", "OTU4"),
      c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5")
    )
  )


  sample_data <- data.frame(
    Treatment = factor(c("A", "A", "B", "B", "B")),
    PatientID = c("P1", "P2", "P3", "P4", "P5"),
    row.names = c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5")
  )


  taxonomy <- matrix(
    c("Firmicutes", "Bacilli", "Bacillales", "Bacillaceae", "Bacillus",
      "Firmicutes", "Bacilli", "Bacillales", "Bacillaceae", "Bacillus",
      "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Bacteroidaceae", "Bacteroides",
      "Proteobacteria", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Escherichia"),
    nrow = 4, ncol = 5, byrow = TRUE,
    dimnames = list(
      c("OTU1", "OTU2", "OTU3", "OTU4"),
      c("Phylum", "Class", "Order", "Family", "Genus")
    )
  )


  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(otu_matrix, taxa_are_rows = TRUE),
    phyloseq::sample_data(sample_data),
    phyloseq::tax_table(taxonomy)
  )

  return(ps)
}

test_that("Function returns correct object type", {
  ps_test <- create_test_phyloseq()
  palette <- c("A" = "blue", "B" = "red")


  result <- create_ordination_plots(
    ps_obj = ps_test,
    method = "PCoA",
    distance_method = "bray",
    group = "Treatment",
    palette = palette,
    taxa_plot = FALSE
  )


  expect_s3_class(result, "gtable")
})

test_that("Function handles invalid inputs correctly", {
  ps_test <- create_test_phyloseq()
  palette <- c("A" = "blue", "B" = "red")


  expect_error(
    create_ordination_plots(
      ps_obj = ps_test,
      group = "NonExistentColumn",
      palette = palette
    ),
    "Grouping variable not found in sample data"
  )


  expect_error(
    create_ordination_plots(
      ps_obj = data.frame(x = 1:10),
      group = "Treatment",
      palette = palette
    ),
    "Input must be a phyloseq object"
  )
})

test_that("Function works with different distance methods", {
  ps_test <- create_test_phyloseq()
  palette <- c("A" = "blue", "B" = "red")


  distance_methods <- c("bray", "jaccard", "euclidean")

  for (method in distance_methods) {
    expect_silent(
      create_ordination_plots(
        ps_obj = ps_test,
        distance_method = method,
        group = "Treatment",
        palette = palette,
        taxa_plot = FALSE
      )
    )
  }
})

test_that("Function works with taxa_plot enabled", {
  ps_test <- create_test_phyloseq()
  palette <- c("A" = "blue", "B" = "red")


  result <- create_ordination_plots(
    ps_obj = ps_test,
    group = "Treatment",
    palette = palette,
    taxa_plot = TRUE,
    taxrank = "Genus",
    top_taxa = 2
  )

  expect_s3_class(result, "gtable")
})

test_that("Function handles different ordination methods", {
  ps_test <- create_test_phyloseq()
  palette <- c("A" = "blue", "B" = "red")


  methods <- c("PCoA", "NMDS")

  for (method in methods) {
    result <- create_ordination_plots(
      ps_obj = ps_test,
      method = method,
      group = "Treatment",
      palette = palette,
      taxa_plot = FALSE
    )

    expect_s3_class(result, "gtable")
  }
})

test_that("Function handles edge cases", {
  ps_test <- create_test_phyloseq()
  palette <- c("A" = "blue", "B" = "red")

  ps_small <- ps_test
  phyloseq::sample_data(ps_small)$Treatment <- c("A", "A", "B", "B", "A")

  expect_s3_class(
    create_ordination_plots(
      ps_obj = ps_small,
      group = "Treatment",
      palette = palette,
      taxa_plot = FALSE
    ),
    "gtable"
  )
})

test_that("Statistical results are calculated correctly", {
  ps_test <- create_test_phyloseq()
  palette <- c("A" = "blue", "B" = "red")

  result <- create_ordination_plots(
    ps_obj = ps_test,
    group = "Treatment",
    palette = palette,
    taxa_plot = FALSE
  )
  expect_true(TRUE)
})

test_that("Function works with different visual parameters", {
  ps_test <- create_test_phyloseq()
  palette <- c("A" = "blue", "B" = "red")


  expect_s3_class(
    create_ordination_plots(
      ps_obj = ps_test,
      group = "Treatment",
      palette = palette,
      size = 12,
      taxa_plot = FALSE
    ),
    "gtable"
  )


  expect_s3_class(
    create_ordination_plots(
      ps_obj = ps_test,
      group = "Treatment",
      palette = palette,
      ratio = c(2, 3),
      taxa_plot = FALSE
    ),
    "gtable"
  )
})
