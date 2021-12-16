#
analysis_plan <- list(
  # do something interesting

  # diversity
  tar_target(
    name = diversity_grad,
    command = comm_raw %>%
    group_by(Gradient, Site, PlotID, Elevation_m) %>%
    summarise(Richness = n(),
              Diversity = diversity(Cover),
              Evenness = Diversity/log(Richness),
              sumAbundance = sum(Cover)) %>%
    pivot_longer(cols = Richness:sumAbundance, names_to = "DiversityIndex", values_to = "Value")
  ),

  # make species ordination
  tar_target(
    name = fNMDS,
    command = make_ordination(comm_raw)
  ),

  # make trait ordination
  tar_target(
    name = trait_pca,
    command = make_trait_pca(trait_mean)
  )
)
