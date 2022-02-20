figures_plan <- list(

  # climate data
  tar_target(
    name = climate_plot,
    command = {
      climate_data %>%
        mutate(Variable = recode(Variable, "SoilMoisture" = "soil moisture in %", "SoilTemperature" = "soil temperature in °C")) %>%
        left_join(coordinates, by = c("Gradient", "Site")) %>%
        ggplot(aes(x = Elevation_m, y = Value, colour = Gradient, fill = Gradient)) +
        geom_point(alpha = 0.5) +
        geom_smooth(method = "lm") +
        scale_colour_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
        scale_fill_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
        labs(x = "Elevation in m a.s.l.", y = "") +
        facet_wrap(~ Variable, scales = "free_y") +
        theme_minimal() +
        theme(legend.position = "bottom")
    }),


  # species ordination
  tar_target(
    name = ordination_plot,
    command = make_ordination_plot(comm_raw,
                                   NMDS = sp_ordination[[1]],
                                   fNMDS = sp_ordination[[2]])),


  # correlation plot
  tar_target(
    name = correlation_plot,
    command = {
      # trait correlation both gradients
      trait_corr_all = get_trait_correlations(trait_mean) %>%
        mutate(r = if_else(r == 1, NA_real_, r)) %>%
        ggplot(aes(x = trait1, y = trait2, fill = r, label = round(r_sig, 2))) +
        geom_tile() +
        labs(x = NULL, y = NULL,
             title = "Both gradients",
             fill = "Pearson's\nCorrelation") +
        scale_fill_gradient2(mid = "#f7f7f7", low = "#f1a340",
                             high = "#998ec3",
                             na.value = "grey80",
                             limits = c(-1, 1)) +
        geom_text() +
        theme_minimal() +
        theme(legend.position = "top") +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0))

      # Bird cliff
      trait_corr_B = get_trait_correlations(trait_mean %>%
                                              filter(Gradient == "B")) %>%
        mutate(r = if_else(r == 1, NA_real_, r)) %>%
        ggplot(aes(x = trait1, y = trait2, fill = r, label = round(r_sig, 2))) +
        geom_tile() +
        labs(x = NULL, y = NULL,
             title = "Bird cliff",
             fill = "Pearson's\nCorrelation") +
        scale_fill_gradient2(mid = "#f7f7f7", low = "#f1a340",
                             high = "#998ec3",
                             na.value = "grey80",
                             limits = c(-1, 1)) +
        geom_text() +
        theme_minimal() +
        theme(legend.position = "none") +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0))

      # Control
      trait_corr_C = get_trait_correlations(trait_mean %>%
                                              filter(Gradient == "C")) %>%
        mutate(r = if_else(r == 1, NA_real_, r)) %>%
        ggplot(aes(x = trait1, y = trait2, fill = r, label = round(r_sig, 2))) +
        geom_tile() +
        labs(x = NULL, y = NULL,
             title = "Reference",
             fill = "Pearson's\nCorrelation") +
        scale_fill_gradient2(mid = "#f7f7f7", low = "#f1a340",
                             high = "#998ec3",
                             na.value = "grey80",
                             limits = c(-1, 1)) +
        geom_text() +
        theme_minimal() +
        theme(legend.position = "none") +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0))

      trait_corr = trait_corr_all / trait_corr_B / trait_corr_C
    })


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


)
