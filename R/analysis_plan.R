#
analysis_plan <- list(
  # do something interesting
  # calculate mean trait value
  tar_target(
    name = trait_mean,
    command = traits_gradient %>%
      group_by(Gradient, Trait) %>%
      summarise(Mean = mean(Value)) %>%
      pivot_wider(names_from = Gradient, values_from = Mean)
    ),

  # diversity
  tar_target(
    name = diversity_grad,
    command = comm_gradient %>%
    group_by(Gradient, Site, PlotID, Elevation_m) %>%
    summarise(Richness = n(),
              Diversity = diversity(Cover),
              Evenness = Diversity/log(Richness),
              sumAbundance = sum(Cover)) %>%
    pivot_longer(cols = Richness:sumAbundance, names_to = "DiversityIndex", values_to = "Value")
  )
)
