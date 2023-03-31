# save figures
library("targets")
library("tarchetypes")
library("dataDownloader")
library("tidyverse")
library("readxl")
library("traitstrap")
library("vegan")
library("ggvegan")
library("kableExtra")
library("viridis")
library("patchwork")
library("broom")
library("ape")
#library("nlme")
library("Hmisc")
library("lme4")
library("broom.mixed")
library("MuMIn")
library(glue)

# Main MS
tar_load(community_trait_plot)
ggsave("output/community_trait_plot.jpg", community_trait_plot, dpi = 300, height = 6, width = 8, bg = "white")

tar_load(community_trait_variance_plot)
ggsave("output/community_trait_variance_plot.jpg", community_trait_variance_plot, dpi = 300, height = 6, width = 8, bg = "white")
tar_load(single_trait_variance_plot)
ggsave("output/single_trait_variance_plot.jpg", single_trait_variance_plot, dpi = 300, height = 6, width = 8, bg = "white")
tar_load(trait_variance_plot)
ggsave("output/trait_variance_plot.jpg", trait_variance_plot, dpi = 300, height = 6, width = 8, bg = "white")


tar_load(trait_ordination_plot)
ggsave("output/trait_ordination_plot.jpg", trait_ordination_plot, dpi = 300, height = 6, width = 6, bg = "white")
tar_load(trait_oridination_PC3)
ggsave("output/trait_oridination_PC3.jpg", trait_oridination_PC3, dpi = 300, height = 6, width = 6, bg = "white")

tar_load(ITV_plot)
ggsave("output/ITV_plot.jpg", ITV_plot, dpi = 300, height = 8, width = 10, bg = "white")

tar_load(bryo_plot)
ggsave("output/bryo_plot.jpg", bryo_plot, dpi = 300, height = 6, width = 6, bg = "white")


# SI figures

tar_load(environment_plot)
ggsave("output/environment_plot.jpg", environment_plot, dpi = 300, height = 6, width = 8, bg = "white")

tar_load(community_pca_plot)
ggsave("output/community_pca_plot.jpg", community_pca_plot, dpi = 300, height = 6, width = 7, bg = "white")

tar_load(diversity_plot)
ggsave("output/diversity_plot.jpg", diversity_plot, dpi = 300, height = 4, width = 6, bg = "white")

tar_load(imputation_plot)
ggsave("output/imputation_plot.jpg", imputation_plot, dpi = 300, height = 6, width = 8, bg = "white")

tar_load(bryophyte_density_plot)
ggsave("output/bryophyte_density_plot.jpg", bryophyte_density_plot, dpi = 300, height = 6, width = 8, bg = "white")


# tar_load(correlation_plot)
# ggsave("output/correlation_plot.jpg", correlation_plot, dpi = 300, height = 10, width = 7, bg = "white")

tar_load(trait_microclimate_figure)
ggsave("output/trait_microclimate_figure.jpg", trait_microclimate_figure, dpi = 300, height = 12, width = 8, bg = "white")
