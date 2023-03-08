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
#library("Hmisc")
library("lme4")
library("broom.mixed")
library("MuMIn")
library(glue)
library(egg)

# Main MS
tar_load(community_trait_plot)
ggsave("output/community_trait_plot.jpg", community_trait_plot, dpi = 300, height = 6, width = 8, bg = "white")

tar_load(community_trait_variance_plot)
ggsave("output/community_trait_variance_plot.jpg", community_trait_variance_plot, dpi = 300, height = 6, width = 8, bg = "white")
tar_load(single_trait_variance_plot)
ggsave("output/single_trait_variance_plot.jpg", single_trait_variance_plot, dpi = 300, height = 6, width = 8, bg = "white")
tar_load(trait_variance_plot)
ggsave("output/trait_variance_plot.jpg", trait_variance_plot, dpi = 300, height = 6, width = 8, bg = "white")

tar_load(full_trait_ordination_plot)
ggsave("output/trait_ordination_plot1.jpg", full_trait_ordination_plot, dpi = 300, height = 6, width = 6, bg = "white")
tar_load(separate_trait_ordination_plot)
ggsave("output/trait_ordination_plot2.jpg", separate_trait_ordination_plot, dpi = 300, height = 6, width = 6, bg = "white")
tar_load(trait_ordination_plot)
ggsave("output/trait_ordination_plot3.jpg", trait_ordination_plot, dpi = 300, height = 6, width = 6, bg = "white")

tar_load(ITV_plot)
ggsave("output/ITV_plot.jpg", ITV_plot, dpi = 300, height = 8, width = 10, bg = "white")

tar_load(bryo_plot)
ggsave("output/bryo_plot.jpg", bryo_plot, dpi = 300, height = 6, width = 6, bg = "white")


# SI figures

tar_load(climate_plot)
ggsave("output/climate_plot.jpg", climate_plot, dpi = 300, height = 4, width = 6, bg = "white")

tar_load(community_pca_plot)
ggsave("output/community_pca_plot.jpg", community_pca_plot, dpi = 300, height = 6, width = 7, bg = "white")

tar_load(ordination_plot)
ggsave("output/ordination_plot.jpg", ordination_plot, dpi = 300, height = 4, width = 6, bg = "white")

tar_load(diversity_plot)
ggsave("output/diversity_plot.jpg", diversity_plot, dpi = 300, height = 4, width = 6, bg = "white")

tar_load(imputation_plot)
ggsave("output/imputation_plot.jpg", imputation_plot, dpi = 300, height = 6, width = 8, bg = "white")

# tar_load(correlation_plot)
# ggsave("output/correlation_plot.jpg", correlation_plot, dpi = 300, height = 10, width = 7, bg = "white")

tar_load(trait_soil_temp_figure)
ggsave("output/trait_soil_temp_figure.jpg", trait_soil_temp_figure, dpi = 300, height = 6, width = 8, bg = "white")


tar_load(trait_soil_moisture_figure)
ggsave("output/trait_soil_moisture_figure.jpg", trait_soil_moisture_figure, dpi = 300, height = 6, width = 8, bg = "white")
