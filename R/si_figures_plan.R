si_figures_plan <- list(

  # CLIMATE DATA
  # analysis
  tar_target(
    name = run_climate_model,
    command = run_trait_model(dat = climate_data,
                              group = "Variable",
                              response = Value,
                              continous_predictor = Elevation_m) |>
      pivot_longer(cols = -c(Variable, data),
                   names_sep = "_",
                   names_to = c(".value", "names"))
  ),

  # select best model
  tar_target(
    name = climate_model,
    command = run_climate_model |>
      # remove models that have singular fit
      filter(singular == FALSE) |>
      filter(aic == min(aic))
  ),

  # lrt
  tar_target(
    name = climate_lrt_output,
    command = climate_lrt(climate_data)
  ),

  # Produce model output and prediction
  tar_target(
    name = climate_model_output,
    command = model_output_prediction(climate_model)
  ),


  # model output
  tar_target(
    name = climate_table,
    command = make_climate_table(climate_model_output)
  ),

  # climate figure
  tar_target(
    name = climate_plot,
    command = make_climate_figure(climate_model_output)
  ),

  # COMMUNITY DATA
  # PCA
  tar_target(
    name = community_pca_plot,
    command = make_sp_pca_figure(comm_pca_B, comm_pca_C)
    ),


  # TRAIT DATA
  # trait imputation plot
  tar_target(
    name = imputation_plot,
    command = {

      gradient <- c(
        B = "Nutrient input",
        C = "Reference")

      trait_names <- c(
        "Plant_Height_cm_log" = "Height cm",
        "Dry_Mass_g_log" = "Dry mass g",
        "Thickness_mm_log" = "Thickness mm",
        "Leaf_Area_cm2_log" = "Area cm2",
        "SLA_cm2_g" = "SLA cm2/g",
        "LDMC" = "LDMC",
        "C_percent" = "C %",
        "N_percent" = "N % ",
        "CN_ratio" = "CN",
        "dN15_permil" = "δN15 ‰",
        "dC13_permil" = "δC13 ‰",
        "P_percent" = "P %",
        "NP_ratio" = "NP")

      #check trait coverage
      imputation_plot <- trait_impute %>%
        autoplot(.) +
        scale_fill_manual(labels = c("Plot", "Site", "Locality"),
                          values = c("#56B4E9", "#009E73", "#E69F00")) +
        scale_y_continuous(breaks = c(0, 0.5, 1)) +
        scale_x_discrete(breaks = c("B_1_D", "B_2_D", "B_3_D", "B_4_D", "B_5_D", "C_1_D", "C_2_D", "C_3_D", "C_4_D", "C_5_D", "C_6_D", "C_7_D"),
                         labels = c("1", "2", "3", "4", "5", "1", "2", "3", "4", "5", "6", "7")) +

        facet_grid(trait_trans ~ Gradient, scales = "free_x",
                   labeller = labeller(Gradient = gradient, trait_trans = trait_names)) +
        labs(x = "Site") +
        theme_minimal() +
        theme(strip.text.y = element_text(angle = 0))

      return(imputation_plot)

    }
  ),

  # # trait ratios
  # tar_target(
  #   name = trait_plot_appendix,
  #   command = community_model_output |>
  #       # remove ratios
  #       filter(trait_trans %in% c("CN_ratio", "NP_ratio")) %>%
  #       make_trait_figure(.)
  #
  # ),

  # trait mean and variance output
  tar_target(
    name = trait_mean_var,
    command = {

      fancy_trait_name_dictionary(trait_mean) %>%
        group_by(Gradient, trait_fancy) %>%
        summarise(se = round(sd(mean)/sqrt(n()), 2),
                  mean = round(mean(mean), 2),
                  se_var = round(sd(var)/sqrt(n()), 2),
                  var = round(mean(var), 2)) %>%
        select(Trait = trait_fancy, Gradient, Mean = mean, "SE Mean" = se, Variance = var, "SE Variance" = se_var) %>%
        write_csv(file = "output/Mean_var.csv")

    }),


  # make species list
  tar_target(
    name = species_list,
    command = {

      bind_rows(vascular = comm_raw |>
                  group_by(Gradient, Taxon) |>
                  summarise(Cover = round(mean(Cover), 1),
                            Cover = as.factor(Cover)),
                bryophyte = bryo_traits_raw |>
                  distinct(Gradient, Taxon) |>
                  mutate(Cover = "x"),
                .id = "Functionalgroup") |>
        pivot_wider(names_from = Gradient, values_from = Cover) |>
        arrange(Functionalgroup, Taxon) |>
        rename("Bird cliff" = "B", "Reference" = "C", "Functional group" = "Functionalgroup") |>
        filter(Taxon != "unknown sp") |>
        write_csv(file = "output/species_list.csv")

    })



#     comm_raw |>
#       distinct(Gradient, Taxon) |>
#       mutate(presence = "x") |>
#       pivot_wider(names_from = Gradient, values_from = presence) |>
#       arrange(Taxon) |>
#       filter(is.na(B) | is.na(C)) |>
#       print(n = Inf)
#
#       comm_raw |> group_by(Gradient, Taxon) |> summarise(m_c = mean(Cover)) |> pivot_wider(names_from = Gradient, values_from = m_c) |>
#         arrange(Taxon) |>View()
#
# traits_raw |>
#   distinct(Gradient, Taxon) |>
#   mutate(presence = "x") |>
#   pivot_wider(names_from = Gradient, values_from = presence) |>
#   arrange(Taxon) |>
#   #filter(is.na(B) | is.na(C)) |>
#   print(n = Inf)



  # correlation plot
  # tar_target(
  #   name = correlation_plot,
  #   command = {
  #     # trait correlation both gradients
  #     trait_corr_all = get_trait_correlations(trait_mean) %>%
  #       mutate(r = if_else(r == 1, NA_real_, r)) %>%
  #       ggplot(aes(x = trait1, y = trait2, fill = r, label = round(r_sig, 2))) +
  #       geom_tile() +
  #       labs(x = NULL, y = NULL,
  #            title = "Both gradients",
  #            fill = "Pearson's\nCorrelation") +
  #       scale_fill_gradient2(mid = "#f7f7f7", low = "#f1a340",
  #                            high = "#998ec3",
  #                            na.value = "grey80",
  #                            limits = c(-1, 1)) +
  #       geom_text() +
  #       theme_minimal() +
  #       theme(legend.position = "top") +
  #       scale_x_discrete(expand=c(0,0)) +
  #       scale_y_discrete(expand=c(0,0))
  #
  #     # Bird cliff
  #     trait_corr_B = get_trait_correlations(trait_mean %>%
  #                                             filter(Gradient == "B")) %>%
  #       mutate(r = if_else(r == 1, NA_real_, r)) %>%
  #       ggplot(aes(x = trait1, y = trait2, fill = r, label = round(r_sig, 2))) +
  #       geom_tile() +
  #       labs(x = NULL, y = NULL,
  #            title = "Bird cliff",
  #            fill = "Pearson's\nCorrelation") +
  #       scale_fill_gradient2(mid = "#f7f7f7", low = "#f1a340",
  #                            high = "#998ec3",
  #                            na.value = "grey80",
  #                            limits = c(-1, 1)) +
  #       geom_text() +
  #       theme_minimal() +
  #       theme(legend.position = "none") +
  #       scale_x_discrete(expand=c(0,0)) +
  #       scale_y_discrete(expand=c(0,0))
  #
  #     # Control
  #     trait_corr_C = get_trait_correlations(trait_mean %>%
  #                                             filter(Gradient == "C")) %>%
  #       mutate(r = if_else(r == 1, NA_real_, r)) %>%
  #       ggplot(aes(x = trait1, y = trait2, fill = r, label = round(r_sig, 2))) +
  #       geom_tile() +
  #       labs(x = NULL, y = NULL,
  #            title = "Reference",
  #            fill = "Pearson's\nCorrelation") +
  #       scale_fill_gradient2(mid = "#f7f7f7", low = "#f1a340",
  #                            high = "#998ec3",
  #                            na.value = "grey80",
  #                            limits = c(-1, 1)) +
  #       geom_text() +
  #       theme_minimal() +
  #       theme(legend.position = "none") +
  #       scale_x_discrete(expand=c(0,0)) +
  #       scale_y_discrete(expand=c(0,0))
  #
  #     trait_corr = trait_corr_all / trait_corr_B / trait_corr_C
  #   })


)
