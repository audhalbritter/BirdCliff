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
  # axis 1 and 2
  tar_target(
    name = trait_ordination_plot,
    command = make_pca_plot(trait_pca)
  ),

  # axis 1 and 3
  tar_target(
    name = trait_oridination_PC3,
    command = make_pca_3_plot(trait_pca)
  ),

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
