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
# Install released version from CRAN
install.packages("pkgdown")
usethis::use_github()
usethis::use_pkgdown_github_pages()
usethis::use_pkgdown()
pkgdown::build_site()
usethis::use_git_config(
  user.name = "Kuzmichenko Polina",
  user.email = "kuzmichenko.p.a@gmail.com"
)
#install.packages("hexSticker")
library(hexSticker)
knitr::include_graphics('figures/octobiom.png',dpi=300)
imgurl <- "./figures/octobiom.png"   # Use your local image path here
sysfonts::font_add_google("Montserrat", "montserrat")
?hexSticker::sticker
sticker(imgurl, package="octobiom", p_size=20, s_x=1, s_y=.75, s_width=.6,h_fill = 'white',p_color = 'black',
        filename="./figures/imgfile.png",h_color="#f39c12")
sticker(
  imgurl,
  package = "octobiom",       # Package name (appears on sticker)
  p_size = 20,                # Text size for package name
  s_x = 1,                    # Logo X position (centered)
  s_y = 0.8,                 # Logo Y position (slightly higher for balance)
  s_width = 0.6,              # Logo width (fill hexagon effectively)
  h_fill = "#F2EEE9",         # Hex fill color (white or leave default)
  h_color = "#A6171C",        # Border color (choose a contrast, e.g., blue or black)
  p_color = "#A6171C",        # Text color (black)
  #p_family = "sans",          # Font family (sans for modern look)
  url = "github.com/META-SBM/octobiom", # Add a URL if desired
  u_size = 4,                 # URL text size
  filename = "./figures/imgfile.png", # Output file location
  p_y=1.5,
  p_family =  "montserrat"
)
usethis::git_sitrep()

