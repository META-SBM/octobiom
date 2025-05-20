library(devtools)
load_all()
devtools::document()
devtools::check()
devtools::install()


library(octobiom)
ps <- readRDS('~/core_microbiome/data/metaphlan.rds')
cols_meta <- read.csv('~/snakemake_project/DATA/cols_meta.csv',row.names = 1)
usethis::use_data(cols_meta)
usethis::use_data(my_palettes, overwrite = TRUE)

load(file='./data/ps.rda')
ps_obj <- prepare_phyloseq(ps, transformation = 'compositional',detection = 0, prevalence= 0,
                           taxonomic_level = 'Species', first = 'filter')
load(file='./data/cols_meta.rda')
permanova_result <- plot_permanova(ps_obj = ps_obj,formula = 'YEAR + bmi_group + age_group',
                                   method = 'bray',cols_meta = cols_meta,palette = get_palette("set_10_2"))
permanova_result[[2]]
my_comparison <- c('2017','2022')
alpha_results <- octobiom::create_alpha_plots(ps_obj = ps_obj,col= 'YEAR',measure = c("Shannon",'Simpson'),method = "wilcox.test",
                   colors = get_palette("set_10_2"))
alpha_results
beta_results <- octobiom::create_ordination_plots(ps_obj,group = 'YEAR',palette = get_palette("set_2_1"))
barplot_hierarchical_results <- barplot_hierarchical(ps_obj,taxrank = 'Species',feature = 'YEAR',top = 15,colors_barplot = get_palette("set_20_2"),
                                                     colors_line = get_palette('set_2_1'),dist = 'bray',size = c(1,1,5),
                                                     decreasing = T)
ps_obj_clr <- prepare_phyloseq(ps, transformation = 'clr',detection = 10, prevalence= 0.05,
                               taxonomic_level = 'Species', first = 'filter')
unique(ps_obj_clr@tax_table[,'Phylum'])
phylum_colors <- c(
  'Euryarchaeota'              = "#001219" ,
  'Thaumarchaeota'             = "#005f73" ,
  'Actinobacteria'             = "#0a9396" ,
  'Bacteria_unclassified'      = "#94d2bd" ,
  'Bacteroidetes'              = "#e9d8a6" ,
  'Candidatus_Melainabacteria' = "#ee9b00" ,
  'Candidatus_Saccharibacteria'= "#ca6702" ,
  'Firmicutes'                 = "#bb3e03" ,
  'Fusobacteria'               = "#ae2012" ,
  'Lentisphaerae'              = "#9b2226" ,
  'Proteobacteria'             = '#F55B66' ,
  'Synergistetes'              = '#ACACAC' ,
  'Tenericutes'                = '#5CD5DA' ,
  'Verrucomicrobia'            = '#0089DF' ,
  'Ascomycota'                 = '#19427D'
)
YEAR_colors <- c(
  '2017' = "#D5B3FF",
  '2022' = "#FBD960"
)
bmi_colors <- c(
  'Lean' = "#FF8400",
  'Overweight' = "#3C096C",
  'Obesity' = "#FF0A54"
)
age_colors <- c(
  '37and45'='#48cae4',#blue
  '46and55' = '#a24f7c' ,
  '56and65'='#ff6b35',#orange
  '65more' = '#ef233c' #red
)
heatmap_results <- octobiom::generate_microbiome_heatmap(ps_obj_clr,tax_level_annotation = 'Phylum',features = 'YEAR',
                                               tax_colors = list('Phylum' = phylum_colors),
                                               feature_colors = list('YEAR' = YEAR_colors))
heatmap_results
heatmap_results <- octobiom::generate_microbiome_heatmap(ps_obj_clr,tax_level_annotation = 'Phylum',features = c('YEAR','bmi_group'),
                                               tax_colors = list('Phylum' = phylum_colors),
                                               feature_colors = list('YEAR' = YEAR_colors,'bmi_group' = bmi_colors))
ps_obj_clr <- octobiom::prepare_phyloseq(ps, transformation = 'clr',detection = 10, prevalence= 0.15,
                               taxonomic_level = 'Species', first = 'filter')
feature_colors <- c(
  "age_at_baseline" = "#ff9999",
  "bmi" = "#99ccff"
)
network_correlation_results <- octobiom::build_multi_correlation_network(
  ps_obj_clr,
  features = c("age_at_baseline", "bmi"),
  tax_rank = "Species",
  cor_method = 'pearson',
  feature_colors = feature_colors
)
octobiom::plot_metadata_distributions(ps_obj_clr,features = c('bmi','age_at_baseline','bmi_group','age_group'),color_palette = octobiom::get_palette("set_10_2"))
network_correlation_results
heatmap_results
?prepare_phyloseq
?plot_permanova
?octobiom::create_alpha_plots
?create_ordination_plots
?barplot_hierarchical
?generate_microbiome_heatmap
?build_multi_correlation_network
?plot_metadata_distributions
get_palette("set_2_2")

ps_obj <- octobiom::prepare_phyloseq(ps,transformation = 'NONE',detection = 10, prevalence = 0.15,taxonomic_level = 'Species',first = 'filter')
ps_obj@sam_data$bmi_group <- factor(ps_obj@sam_data$bmi_group,levels = c('Lean','Overweight','Obesity'))
ps_obj@sam_data$age_group <- ifelse(ps_obj@sam_data$age_at_baseline <= 36, 'less_36',
                                    ifelse(ps_obj@sam_data$age_at_baseline <= 45, '37and45',
                                           ifelse(ps_obj@sam_data$age_at_baseline <= 55, '46and55',
                                                  ifelse(ps_obj@sam_data$age_at_baseline <= 65, '56and65',
                                                         '65more'))))
diagdds <- phyloseq_to_deseq2(ps_obj ,~ age_group)
diagdds <- DESeq(diagdds, test="Wald", fitType="parametric", sfType = "poscounts")

contrasts <- list(
  c('bmi_group', 'Lean', 'Overweight'),
  c('bmi_group', 'Lean', 'Obesity'),
  c('bmi_group', 'Overweight', 'Obesity')
)
contrasts <- list(
  c('age_group', '46and55', '37and45'),
  c('age_group', '56and65', '37and45'),
  c('age_group', '65more', '37and45')
)
resultsNames(diagdds)
# Применяем функцию с обработкой NULL результатов
sigtabs_list <- lapply(contrasts, function(contrast) {
  tryCatch(
    octobiom::create_sigtab(ps_obj, diagdds, contrast),
    error = function(e) {
      message(paste("Error processing contrast:", paste(contrast, collapse = " vs "), " - ", e$message))
      return(NULL)
    }
  )
})
?octobiom::create_sigtab
combined_df <- do.call(rbind, sigtabs_list)
?octobiom::create_comparison_plots
combined_df
comparison_names <- unique(combined_df$comparison)
diff_plots <- octobiom::create_comparison_plots(data = combined_df, comparison_name = comparison_names[1],group_var_abund_prev = "bmi_group",
                                  group_var_fold_change = "Phylum",
                                  ps_obj = ps_obj,
                                  group_colors_abund_prev =  bmi_colors,
                                  group_colors_fold_change = phylum_colors)
diff_plots <- octobiom::create_comparison_plots(data = combined_df, comparison_name = comparison_names[1],group_var_abund_prev = "age_group",
                                                group_var_fold_change = "Phylum",
                                                ps_obj = ps_obj,
                                                group_colors_abund_prev =  age_colors,
                                                group_colors_fold_change = phylum_colors)
usethis::use_testthat()
usethis::use_test('permanova.R')
devtools::test()
usethis::use_vignette("my-vignette")
?octobiom::retention_plot
det <- c(0, 0.01,0.05,0.1)/100
prevalences <- seq(.5, 1, .1)
ps_obj <- octobiom::prepare_phyloseq(ps, transformation = 'compositional',detection = 0, prevalence= 0,
                                     taxonomic_level = 'Species', first = 'filter')
octobiom::retention_plot(phyloseq = ps_obj,detection_values = det, prevalence_values = prevalences,
                         prevalence_colors = get_palette('set_10_4'))
?octobiom::prevalence_abundance_plot
octobiom::prevalence_abundance_plot(ps,top_n=10)
?octobiom::group_prevalence

