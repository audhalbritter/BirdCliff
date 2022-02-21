# save figures

# Main MS
tar_load(trait_plot)
ggsave("output/trait_plot.jpg", trait_plot, dpi = 300, height = 6, width = 8, bg = "white")

tar_load(trait_ordination_plot)
ggsave("output/trait_ordination_plot.jpg", trait_ordination_plot, dpi = 300, height = 7, width = 10, bg = "white")

tar_load(ind_species_figure)
ggsave("output/ind_species_figure.jpg", ind_species_figure, dpi = 300, height = 6, width = 10, bg = "white")



# SI figures

tar_load(climate_plot)
ggsave("output/climate_plot.jpg", climate_plot, dpi = 300, height = 4, width = 6, bg = "white")

tar_load(ordination_plot)
ggsave("output/ordination_plot.jpg", ordination_plot, dpi = 300, height = 5, width = 10, bg = "white")

tar_load(diversity_plot)
ggsave("output/diversity_plot.jpg", diversity_plot, dpi = 300, height = 4, width = 6, bg = "white")

tar_load(correlation_plot)
ggsave("output/correlation_plot.jpg", correlation_plot, dpi = 300, height = 10, width = 7, bg = "white")

