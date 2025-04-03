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
library("Hmisc")
library("lme4")
library("broom.mixed")
library("MuMIn")
library(glue)

# Main MS
tar_load(community_trait_plot)
ggsave("output/community_trait_plot.pdf", community_trait_plot, dpi = 300, height = 6, width = 8, bg = "white", device = cairo_pdf)

tar_load(dN15_plot)
ggsave("output/dn15_plot.pdf", dN15_plot, dpi = 300, height = 6, width = 8, bg = "white", device = cairo_pdf)


tar_load(trait_ordination_plot)
#ggsave("output/trait_ordination_plot2.jpg", trait_ordination_plot, dpi = 300, height = 8, width = 8.5, bg = "white")
ggsave("output/trait_ordination_plot4.pdf", trait_ordination_plot, dpi = 300, height = 8, width = 6, bg = "white", device = cairo_pdf)
tar_load(trait_oridination_PC3)
ggsave("output/trait_oridination_PC3.pdf", trait_oridination_PC3, dpi = 300, height = 8, width = 8.5, bg = "white", device = cairo_pdf)

tar_load(ITV_plot)
ggsave("output/ITV_plot_new.pdf", ITV_plot, dpi = 300, height = 8, width = 10, bg = "white", device = cairo_pdf)

tar_load(vegetation_plot)
ggsave("output/vegetation_plot.pdf", vegetation_plot, dpi = 300, height = 8, width = 10, bg = "white", device = cairo_pdf)

# tar_load(bryo_plot)
# ggsave("output/bryo_plot.jpg", bryo_plot, dpi = 300, height = 6, width = 6, bg = "white")


# SI figures

tar_load(environment_plot)
ggsave("output/environment_plot.pdf", environment_plot, dpi = 300, height = 6, width = 8, bg = "white", device = cairo_pdf)

tar_load(community_pca_plot)
ggsave("output/community_pca_plot2.pdf", community_pca_plot, dpi = 300, height = 5, width = 7, bg = "white", device = cairo_pdf)

# tar_load(diversity_plot)
# ggsave("output/diversity_plot.jpg", diversity_plot, dpi = 300, height = 4, width = 6, bg = "white")

tar_load(imputation_plot)
ggsave("output/imputation_plot.pdf", imputation_plot, dpi = 300, height = 6, width = 8, bg = "white", device = cairo_pdf)

# tar_load(bryophyte_density_plot)
# ggsave("output/bryophyte_density_plot.jpg", bryophyte_density_plot, dpi = 300, height = 6, width = 8, bg = "white")


# tar_load(correlation_plot)
# ggsave("output/correlation_plot.jpg", correlation_plot, dpi = 300, height = 10, width = 7, bg = "white")

tar_load(trait_microclimate_figure)
ggsave("output/trait_microclimate_figure.pdf", trait_microclimate_figure, dpi = 300, height = 12, width = 8, bg = "white", device = cairo_pdf)
