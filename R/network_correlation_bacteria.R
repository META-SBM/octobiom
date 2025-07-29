#' Run FastSpar Network Analysis Pipeline
#'
#' This function performs a complete microbial co-occurrence network analysis using FastSpar,
#' including correlation calculation, bootstrap analysis, p-value estimation, and network visualization.
#'
#' @param conda_env Character. Name of the conda environment containing FastSpar (default: "fastspar").
#' @param input_dir Character. Path to directory containing input OTU table (default: ".").
#' @param output_dir Character. Path to directory for output files (default: "results").
#' @param otu_table Character. Filename of OTU table (default: "otu_table.tsv").
#' @param iterations Integer. Number of iterations for FastSpar (default: 50).
#' @param exclusion_iterations Integer. Number of exclusion iterations (default: 10).
#' @param threshold Numeric. Correlation strength exclusion threshold (default: 0.1).
#' @param threads Integer. Number of threads for parallel processing (default: 1).
#' @param seed Integer. Random number generator seed (default: 1).
#' @param permutations Integer. Number of bootstrap permutations (default: 1000).
#' @param cor_threshold Numeric. Minimum absolute correlation for edge inclusion (default: 0.3).
#' @param p_adjust_method Character. Method for p-value adjustment (default: "BH").
#' @param remove_temp Logical. Remove temporary bootstrap files? (default: TRUE).
#' @param ps_obj phyloseq object. Optional phyloseq object for CLR abundance calculation.
#'
#' @return A list containing:
#' \itemize{
#'   \item graph - igraph network object
#'   \item adjacency_matrix - Filtered adjacency matrix
#'   \item correlation_matrix - Raw correlation matrix
#'   \item pvalues - Adjusted p-value matrix
#'   \item plot - ggraph visualization of the network
#' }
#'
#' @details The function performs the following steps:
#' \enumerate{
#'   \item Calculates correlations using FastSpar
#'   \item Performs bootstrap analysis
#'   \item Estimates p-values
#'   \item Filters correlations based on significance and threshold
#'   \item Constructs and visualizes the network
#'   \item Calculates network metrics (modularity, connectivity, etc.)
#' }
#'
#' @examples
#' \dontrun{
#' result <- run_fastspar_analysis(
#'   conda_env = "fastspar",
#'   input_dir = "path/to/data",
#'   ps_obj = phyloseq_object
#' )
#' print(result$plot)
#' }
#'
#' @export
#' @importFrom igraph graph_from_data_frame set_edge_attr edge_density degree cluster_louvain
#' @importFrom ggraph ggraph geom_edge_arc geom_node_point geom_node_text scale_edge_width
#' @importFrom phyloseq otu_table
#' @importFrom microbiome transform
#' @importFrom dplyr rename filter mutate
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot_build theme_void scale_size_continuous
#' @importFrom stats p.adjust
run_fastspar_analysis <- function(
    conda_env = "fastspar",
    input_dir = ".",
    output_dir = "results",
    otu_table = "otu_table.tsv",
    iterations = 50,
    exclusion_iterations = 10,
    threshold = 0.1,
    threads = 1,
    seed = 1,
    permutations = 1000,
    cor_threshold = 0.3,
    p_adjust_method = "BH",
    remove_temp = TRUE,
    ps_obj = ps_obj) {

  # Load required packages
  require(igraph)
  require(ggraph)
  require(tidyverse)
  require(phyloseq)
  require(microbiome)

  # Create output directory
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Activate conda environment
  #system(paste("conda init"))
  #system(paste("conda activate", conda_env))
  system(paste('conda run -n', conda_env))
  tryCatch({
    # Убедимся, что все пути абсолютные
    input_dir <- normalizePath(input_dir)
    output_dir <- normalizePath(output_dir)
    otu_table_path <- file.path(input_dir, otu_table)

    correlation_file <- file.path(output_dir, "median_correlation.tsv")
    covariance_file <- file.path(output_dir, "covariance.tsv")
    pvalues_file <- file.path(output_dir, "pvalues.tsv")

    # Проверка существования OTU таблицы
    if (!file.exists(otu_table_path)) {
      stop("OTU table not found: ", otu_table_path)
    }
    if (file.exists(correlation_file) && file.exists(covariance_file) && file.exists(pvalues_file)) {
      message("Output files already exist. Skipping FastSpar run and proceeding to network construction.")
    } else {
      message("Output files not found. Running FastSpar analysis...")

      # Step 1: Run FastSpar
      fastspar_cmd <- sprintf(
        "conda run -n %s fastspar --otu_table %s --correlation %s --covariance %s -i %d -x %d -e %.2f -t %d -s %d",
        conda_env, otu_table_path, correlation_file, covariance_file,
        iterations, exclusion_iterations, threshold, threads, seed
      )
      message("Running: ", fastspar_cmd)
      system(fastspar_cmd)

      # Step 2: Bootstrap
      bootstrap_counts_dir <- file.path(output_dir, "bootstrap_counts")
      if (!dir.exists(bootstrap_counts_dir)) dir.create(bootstrap_counts_dir, recursive = TRUE)

      fastspar_bootstrap_cmd <- sprintf(
        "conda run -n %s fastspar_bootstrap --otu_table %s --number %d --prefix %s/data",
        conda_env, otu_table_path, permutations, bootstrap_counts_dir
      )
      message("Running: ", fastspar_bootstrap_cmd)
      system(fastspar_bootstrap_cmd)

      # Step 3: Parallel
      bootstrap_correlation_dir <- file.path(output_dir, "bootstrap_correlation")
      if (!dir.exists(bootstrap_correlation_dir)) dir.create(bootstrap_correlation_dir)

      fastspar_parallel_cmd <- sprintf(
        "conda run -n %s parallel fastspar --otu_table {} --correlation %s/cor_{/} --covariance %s/cov_{/} -i 5 ::: %s/*",
        conda_env, bootstrap_correlation_dir, bootstrap_correlation_dir, bootstrap_counts_dir
      )
      message("Running: ", fastspar_parallel_cmd)
      system(fastspar_parallel_cmd)

      # Step 4: P-values
      fastspar_pvalues_cmd <- sprintf(
        "conda run -n %s fastspar_pvalues --otu_table %s --correlation %s --prefix %s/cor_data_ --permutations %d --outfile %s",
        conda_env, otu_table_path, correlation_file, bootstrap_correlation_dir, permutations, pvalues_file
      )
      message("Running: ", fastspar_pvalues_cmd)
      system(fastspar_pvalues_cmd)

      # Clean up temporary files
      if (remove_temp) {
        if (dir.exists(bootstrap_counts_dir)) unlink(bootstrap_counts_dir, recursive = TRUE, force = TRUE)
        if (dir.exists(bootstrap_correlation_dir)) unlink(bootstrap_correlation_dir, recursive = TRUE, force = TRUE)
      }
    }})

# Read and process results
tryCatch({
  # Read correlation matrix
  NPS.corr <- as.matrix(read.table(file.path(output_dir, "median_correlation.tsv"),
                                   row.names = 1))
  colnames(NPS.corr) <- rownames(NPS.corr)

  # Read p-values
  NPS.oneside <- as.matrix(read.table(file.path(output_dir, "pvalues.tsv"),
                                                            row.names = 1))
  colnames(NPS.oneside) <- rownames(NPS.oneside)

  # Adjust p-values
  NPS.oneside.adjusted <- matrix(
    p.adjust(as.vector(NPS.oneside), method = p_adjust_method),
    nrow = nrow(NPS.oneside),
    dimnames = dimnames(NPS.oneside)
  )

  # Create adjacency matrix
  adj_matrix <- NPS.corr
  adj_matrix[NPS.oneside.adjusted >= 0.05 | abs(NPS.corr) < cor_threshold] <- 0
  if (all(adj_matrix == 0)) {
    stop("The adjacency matrix is completely zero after applying the threshold conditions. Function terminated.")
  }
  adj_matrix <- adj_matrix[rowSums(abs(adj_matrix)) > 0, colSums(abs(adj_matrix)) > 0]


  clr_data <- microbiome::transform(ps_obj, "clr")
  clr_matrix <- as(t(otu_table(clr_data)), "matrix")
  mean_abundance <- rowMeans(clr_matrix)
  mean_abundance <- mean_abundance %>%
    as.data.frame()
  mean_abundance <- mean_abundance %>%
    filter(rownames(mean_abundance) %in% rownames(adj_matrix))%>%
    rename(abundance='.')

  # Add edge attributes
  results_df <- adj_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column("from") %>%
    pivot_longer(cols = -from, names_to = "to", values_to = "cor") %>%
    #filter(cor != 0) %>%
    mutate(
      direction = ifelse(cor > 0, "positive", "negative"),
      cor_abs = abs(cor)
    )
  # Create graph
  g <- igraph::graph_from_data_frame(
    d = results_df,
    directed = FALSE,
    vertices = data.frame(
      name = row.names(mean_abundance),
      abundance = mean_abundance
    )
  )

  g <- set_edge_attr(g, "edge_direction", value = results_df$direction)
  g <- set_edge_attr(g, "edge_cor_abs", value = results_df$cor_abs)
  g <- set_edge_attr(g, "edge_cor_value", value = results_df$cor)

  # Network metrics
  cl <- cluster_louvain(g)
  modularity_value <- round(modularity(cl),2)
  natural_connectivity <- round(pulsar::natural.connectivity(adj_matrix),4)
  edge_density_value <- round(edge_density(g),2)
  edge_cor_values <- E(g)$edge_cor_value
  mean_edge_cor <- round(mean(edge_cor_values),3)
  positive_edge_percentage <- round(sum(edge_cor_values > 0) / length(edge_cor_values) * 100,2)
  num_components <- igraph::components(g)$no
  clustering_coeff <- round(transitivity(g, type = "average"),2)
  total_edges <- ecount(g)
  total_nodes <- vcount(g)
  edges_to_nodes_ratio <- round(total_edges / total_nodes,2)
  network_annotation <- paste("Modularity:", modularity_value,"\n", "Natural connectivity:", natural_connectivity,"\n",
                              "Number of components:", num_components,"\n", "Positive edge percentage:", positive_edge_percentage,"\n",
                              "Edge density:", edge_density_value,"\n","Clustering coefficient:", clustering_coeff,"\n",
                              "Average degree:", edges_to_nodes_ratio,"\n", "Mean correlation:", mean_edge_cor)
  network_stats <- data.frame(
    Metric = c("Modularity",
               "Natural connectivity",
               "Number of components",
               "Positive edge percentage",
               "Edge density",
               "Clustering coefficient",
               "Average degree (Edges/Nodes ratio)",
               "Mean correlation"),
    Value = c(modularity_value,
              natural_connectivity,
              num_components,
              positive_edge_percentage,
              edge_density_value,
              clustering_coeff,
              edges_to_nodes_ratio,
              mean_edge_cor),
    stringsAsFactors = FALSE
  )


  # Visualization
  p <- ggraph(g, layout = "fr") +
    geom_edge_arc(
      aes(width = edge_cor_abs, color = edge_direction),
      strength = 0.1, alpha = 0.7) +
    geom_node_point(aes(size = abundance), alpha = 0.3, color = "black") +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_edge_color_manual(values = c(positive = "red", negative = "blue")) +
    scale_edge_width(range = c(0.5, 2)) +
    scale_size_continuous(
      range = c(2, 10),
      name = "Mean CLR Abundance",
      breaks = seq(round(min(mean_abundance),1), round(max(mean_abundance),1), length.out = 5)
    ) +
    theme_void() +
    labs(title = "Microbial Co-occurrence Network",
         subtitle = 'SparCC algorithm')+
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  x_range <- ggplot_build(p)$layout$panel_params[[1]]$x.range
  y_range <- ggplot_build(p)$layout$panel_params[[1]]$y.range


  x_pos <- x_range[2]
  y_pos <- y_range[2]


  p <- p +
    geom_text(
      aes(x = x_pos, y = y_pos, label = network_annotation),
      hjust = 1, vjust = 1,alpha = 0.1,
      size = 4
    )

  return(list(
    graph = g,
    adjacency_matrix = adj_matrix,
    correlation_matrix = NPS.corr,
    pvalues = NPS.oneside.adjusted,
    plot = p,
    network_stats = network_stats
  ))

}, error = function(e) {
  stop("Error processing results: ", e$message)
})
}
