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

library(hexSticker)
library(showtext)
library(sysfonts)


font_add_google("Montserrat", "montserrat")
font_add_google("Roboto", "roboto")
showtext_auto()
?hexSticker::sticker
sticker(imgurl, package="octobiom", p_size=20, s_x=1, s_y=.75, s_width=.6,h_fill = 'white',p_color = 'black',
        filename="./figures/imgfile.png",h_color="#f39c12")
sticker(
  imgurl,
  package = "octobiom",
  p_size = 18,
  s_x = 1,
  s_y = 0.78,
  s_width = 0.62,
  h_fill = "#FFFFFF",
  h_color = "#3498DB",
  h_size = 1.5,
  filename = "./figures/imgfile.png",
  p_y=1.5,
  p_color = "#2C3E50",
  p_family = "montserrat",
  p_fontface = "bold",
  url = "github.com/META-SBM/octobiom",
  u_size = 4.5,
  u_color = "#7F8C8D",
  u_y = 0.08,
  u_family = "roboto",
)
usethis::git_sitrep()

