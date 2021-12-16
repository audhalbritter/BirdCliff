figures_plan <- list(

  # climate data
  tar_target(
    name = climate_plot,
    command = {
      climate_data %>%
        mutate(Variable = recode(Variable, "SoilMoisture" = "soil moisture in %", "SoilTemperature" = "soil temperature in °C")) %>%
        left_join(coordinates, by = c("Gradient", "Site")) %>%
        ggplot(aes(x = Elevation_m, y = Value, colour = Gradient)) +
        geom_point(alpha = 0.5) +
        geom_smooth(method = "lm") +
        scale_colour_manual(values = c("green4", "grey"), labels = c("Birdcliff", "Control")) +
        labs(x = "Elevation in m a.s.l.", y = "") +
        facet_wrap(~ Variable, scales = "free_y") +
        theme_minimal()
      }),

  # trait change along gradients
  tar_target(
    name = trait_plot,
    command = {
      fancy_trait_name_dictionary(trait_mean) %>%
        ggplot(aes(x = Elevation_m, y = mean, colour = Gradient)) +
        geom_point(alpha = 0.3) +
        geom_smooth(method = "lm") +
        scale_colour_manual(values = c("green4", "grey"), labels = c("Birdcliff", "Control")) +
        labs(x = "Elevation in m a.s.l.", y = "Bootstrapped trait mean") +
        facet_wrap(~ trait_fancy, scales = "free_y") +
        theme_minimal() +
        theme(legend.position = c(0.9, 0.05))

    }),

  # tar_target(
  #   name = trait_histogram,
  #   command = {
  #
  #     trait_plot <- traits_gradient %>%
  #       mutate(Value = if_else(Trait %in% c("Dry_Mass_g", "Wet_Mass_g", "Leaf_Area_cm2", "Plant_Height_cm"), log(Value), Value),
  #              Trait = factor(Trait, levels = c("Plant_Height_cm", "Wet_Mass_g", "Dry_Mass_g", "Leaf_Area_cm2", "Leaf_Thickness_mm", "SLA_cm2_g", "LDMC", "C_percent", "N_percent", "P_percent", "CN_ratio", "dC13_permil", "dN15_permil"))) %>%
  #         filter(!is.na(Trait)) %>%
  #       ggplot(aes(x = Value, fill = Gradient, colour = Gradient)) +
  #       geom_density(alpha = 0.5) +
  #       scale_fill_manual(values = c("grey", "green4"), labels = c("Control", "Birdcliff")) +
  #       scale_colour_manual(values = c("grey", "green4"), labels = c("Control", "Birdcliff")) +
  #       labs(x = "Trait values", y = "Density") +
  #       facet_wrap(~Trait, scales = "free")
  #
  #   }),

  # trait histogram
  tar_target(
    name = diversity_plot,
    command = {
      ggplot(diversity_grad, aes(x = Elevation_m, y = Value, colour = Gradient, fill = Gradient)) +
        geom_point() +
        geom_smooth(method = "lm") +
        scale_colour_manual(values = c("green4", "grey"), labels = c("Birdcliff", "Control")) +
        scale_fill_manual(values = c("green4", "grey"), labels = c("Birdcliff", "Control")) +
        labs(x = "Elelvation in m a.s.l.", y = "") +
        facet_wrap(~ DiversityIndex, scales = "free")
    })

)

