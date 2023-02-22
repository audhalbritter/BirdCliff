figures_plan <- list(

  ### TRAITS

  # FIGURE 2: trait change along gradients
  tar_target(
    name = community_trait_plot,
    command = community_model_output %>%
        # remove ratios
        #filter(!trait_trans %in% c("CN_ratio", "NP_ratio")) %>%
        make_trait_figure(.)
  ),

  # variance
  tar_target(
    name = community_trait_variance_plot,
    command = make_trait_variance_figure(community_variance_output)
  ),

  # density curves
  tar_target(
    name = single_trait_variance_plot,
    command = trait_mean |>
      ggplot(aes(x = mean, fill = Gradient)) +
      geom_density(alpha = 0.4) +
      scale_fill_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
      labs(x = "Bootstrapped trait mean",
           y = "Density") +
      facet_wrap( ~ trait_trans, scales = "free") +
      theme_minimal()
  ),

  # variation in PC1 scores
  tar_target(
    name = trait_variance_plot,
    command = trait_pca[[1]] |>
      ggplot(aes(x = factor(Site), y = PC1, fill = Gradient, colour = Gradient)) +
      geom_violin(alpha = 0.3) +
      geom_point(position=position_jitterdodge()) +
      scale_fill_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
      scale_colour_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
      scale_x_discrete(labels = c("1" = "10.8", "2" = "41.8", "3" = "84.0", "4" = "123.0", "5" = "174.0", "6" = "224.0", "7" = "238.0")) +
      labs(x = "Elevation m a.s.l.",
           y = "PC1 score") +
      theme_minimal()
  ),

  # FIGURE 3: trait ordination
  tar_target(
    name = full_trait_ordination_plot,
    command = make_full_trait_pca_plot(trait_pca)
  ),

  tar_target(
    name = separate_trait_ordination_plot,
    command = make_separate_pca_plot(trait_pca)
  ),

  # separate ordinations (appendix!!!)
  tar_target(
    name = trait_ordination_plot,
    command = make_trait_pca_plot(trait_pca_B, trait_pca_C)
  ),

  # FIGURE 5: VASCULAR PLANTS
  # tar_target(
  #   name = vascular_plot,
  #   command = make_vascular_figure(vascular_model_output)),


  # FIGURE 6: BRYOPHYTES
  tar_target(
    name = bryo_plot,
    command = make_bryo_figure(bryophyte_model_output)),

  # FIGURE 4: ITV PLOT
  tar_target(
    name = ITV_plot,
    command = make_ITV_plot(itv_output)
  )

)



# # soilM and T figures
# tar_target(
#   name = soilM_plot,
#   command = {
#
#     bm = fancy_trait_name_dictionary(trait_mean) |>
#       distinct(trait_fancy) |>
#       mutate(best_model = c("", "G", "", "G", rep("", 3), "G", rep("", 5)),
#              y.pos = c(45, 25, -30, 16, -5, 0.5, 1, 3.8, 55, 0.3, 2, 200, -1))
#
#     soilM <- fancy_trait_name_dictionary(trait_mean) |>
#       left_join(bm, by = "trait_fancy") |>
#       ggplot(aes(x = SoilMoisture, y = mean, colour = Gradient)) +
#       geom_point(alpha = 0.5) +
#       geom_smooth(method = "lm") +
#       scale_colour_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
#       geom_text(aes(label = best_model, x = 10, y = y.pos), colour = "black") +
#       labs(x = "Soil moisture in %", y = "Trait mean") +
#       facet_wrap(~ trait_fancy, scales = "free_y")
#
#     ggsave(soilM, filename = "output/soilM.jpeg", dpi = 300, width = 8, height = 6)
#
#
#     soilT <- fancy_trait_name_dictionary(trait_mean) |>
#       left_join(bm, by = "trait_fancy") |>
#       ggplot(aes(x = SoilTemperature, y = mean, colour = Gradient)) +
#       geom_point(alpha = 0.5) +
#       geom_smooth(method = "lm") +
#       scale_colour_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
#       geom_text(aes(label = best_model, x = 10, y = y.pos), colour = "black") +
#       labs(x = "Soil temperature in°C", y = "Trait mean") +
#       facet_wrap(~ trait_fancy, scales = "free_y")
#
#     ggsave(soilT, filename = "output/soilT.jpeg", dpi = 300, width = 8, height = 6)
#
#   }),
