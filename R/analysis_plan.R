#
analysis_plan <- list(

  # CLIMATE
  tar_target(
    name = climate_result,
    command = climate_data |>
      group_by(Variable, Gradient) |>
      summarise(mean = mean(Value),
                se = sd(Value)/sqrt(n()))
  ),


  # COMMUNITY
  # Species community PCA
  # make species ordination (separate by Gradient)
  tar_target(
    name = comm_pca_B,
    command = make_community_pca(comm_raw |>
                                   filter(Gradient == "B"))
  ),

  tar_target(
    name = comm_pca_C,
    command = make_community_pca(comm_raw |>
                                   filter(Gradient == "C"))
  ),


  # FUNCTIONAL TRAITS
  # Community trait mean

  # Elevation
  # run linear and quadratic model
  tar_target(
    name = trait_community_model,
    command = run_trait_model(dat = trait_mean,
                               group = "trait_trans",
                               response = mean,
                               continous_predictor = Elevation_m) |>
      pivot_longer(cols = -c(trait_trans, data),
                   names_sep = "_",
                   names_to = c(".value", "names"))
  ),

  # select best model
  tar_target(
    name = community_model,
    command = trait_community_model |>
      # remove models that have singular fit
      filter(singular == FALSE) |>
      filter(aic == min(aic))
  ),

  # Likelihood ratio test
  tar_target(
    name = community_lrt,
    command = likelihood_ratio_test(trait_mean)
  ),

  # Produce model output and prediction
  tar_target(
    name = community_model_output,
    command = model_output_prediction(community_model) |>
      # add LRT text
      mutate(text = case_when(trait_trans %in% c("NP_ratio", "dC13_permil") ~ "Null",
                              trait_trans %in% c("CN_ratio", "N_percent") ~ "N+E",
                              TRUE ~ "NxE"))
    ),

  # trait table
  tar_target(
    name = trait_model_table,
    command = make_trait_table(community_model_output)
    ),

  # MICROCLIMATE
  # run linear and quadratic model
  tar_target(
    name = trait_soil_temp_model,
    command = run_trait_model(dat = trait_mean,
                              group = "trait_trans",
                              response = mean,
                              continous_predictor = SoilTemperature) |>
      pivot_longer(cols = -c(trait_trans, data),
                   names_sep = "_",
                   names_to = c(".value", "names")) |>
      filter(singular == "FALSE") |>
      filter(aic == min(aic))
  ),

  # Likelihood ratio test
  tar_target(
    name = soiltemp_lrt,
    command = likelihood_ratio_test_ST(trait_mean)  #|>
    #filter(`Pr(>Chisq)` <= 0.05)
  ),

  # Produce model output and prediction
  # others: NULL
  # PH, dC13: NxST
  # N%, CN, dN15: ST
  tar_target(
    name = soil_temp_model_output,
    command = model_output_prediction(trait_soil_temp_model) |>
      # add LRT text
      mutate(text = case_when(trait_trans %in% c("Plant_Height_cm_log", "dC13_permil") ~ "Temperature",
                              trait_trans %in% c("CN_ratio", "N_percent", "dN15_permil") ~ "NxTemperature",
                              TRUE ~ "Null"))
  ),


  # run linear and quadratic model
  tar_target(
    name = trait_soil_moisture_model,
    command = run_trait_model(dat = trait_mean,
                              group = "trait_trans",
                              response = mean,
                              continous_predictor = SoilMoisture) |>
      pivot_longer(cols = -c(trait_trans, data),
                   names_sep = "_",
                   names_to = c(".value", "names")) |>
      filter(singular == "FALSE") |>
      filter(aic == min(aic))
  ),

  # Likelihood ratio test
  tar_target(
    name = soilmoisture_lrt,
    command = likelihood_ratio_test_SM(trait_mean)  #|>
    #filter(`Pr(>Chisq)` <= 0.05)
  ),

  # Produce model output and prediction
  # others: NULL
  # N%, CN, dN15: SM
  tar_target(
    name = soil_moisture_model_output,
    command = model_output_prediction(trait_soil_moisture_model) |>
      # add LRT text
      mutate(text = case_when(trait_trans %in% c("CN_ratio", "N_percent", "dN15_permil") ~ "Moisture",
                              TRUE ~ "Null"))
  ),


  # Community trait variance
  # run linear and quadratic model
  # NP ratio, both models are singular fit
  tar_target(
    name = trait_variance_model,
    command = run_trait_model(dat = trait_mean,
                              group = "trait_trans",
                              response = var,
                              continous_predictor = Elevation_m) |>
      pivot_longer(cols = -c(trait_trans, data),
                   names_sep = "_",
                   names_to = c(".value", "names"))
  ),

  # select best model
  tar_target(
    name = community_variance_model,
    command = trait_variance_model |>
      # remove models that have singular fit
      filter(singular == FALSE) |>
      filter(aic == min(aic))
  ),

  # LRT for variance
  tar_target(
    name = ltr_variance,
    command = likelihood_ratio_test_variance(trait_mean |>
                                               filter(trait_trans != "NP_ratio"))
  ),

  # Produce model output and prediction
  tar_target(
    name = community_variance_output,
    command = model_output_prediction(community_variance_model) |>
      # add LRT text
      mutate(text = case_when(trait_trans %in% c("SLA_cm2_g", "dC13_permil") ~ "Null",
                              trait_trans %in% c("Plant_Height_cm_log") ~ "E",
                              TRUE ~ "NxE"))
  ),


  # Trait ordination (PCA)
  tar_target(
    name = trait_pca_B,
    command = make_trait_pca(trait_mean %>% filter(Gradient == "B"))
  ),

  tar_target(
    name = trait_pca_C,
    command = make_trait_pca(trait_mean %>% filter(Gradient == "C"))
  ),

  # pca output
  tar_target(
    name = trait_pca_output,
    command = {

      bind_rows(Birdcliff = trait_pca_B[[2]],
                Reference = trait_pca_C[[2]], #%>%
                  # mutate(PC1 = PC1 * -1,
                  #        PC2 = PC2 * -1,
                  #        PC3 = PC3 * -1,
                  #        PC4 = PC4 * -1),
                .id = "Gradient") %>%
        select(Gradient, Trait = trait_fancy, PC1:PC4) %>%
        mutate(PC1 = round(PC1, digits = 2),
               PC2 = round(PC2, digits = 2),
               PC3 = round(PC3, digits = 2),
               PC4 = round(PC4, digits = 2)) %>%
        arrange(Gradient, -PC1) %>%
        # order traits
        write_csv(., file = "output/Loadings_trait_PCA.csv")

    }
  ),

  tar_target(
    name = trait_ord_expl_var,
    command = {

      bind_cols(
        bird = vegan::eigenvals(trait_pca_B[[3]])/sum(vegan::eigenvals(trait_pca_B[[3]])) * 100,
        reference = vegan::eigenvals(trait_pca_C[[3]])/sum(vegan::eigenvals(trait_pca_C[[3]])) * 100)

    }),


  ### ITV
  tar_target(
    name = itv_output,
    command = make_ITV_analysis(trait_mean)),


  # INDIVIDUAL LEVEL TRAITS

  # Vascular plants

  # combine data
  tar_target(
    name = ind_traits,
    command = combine_traits(traits_raw, bryo_traits_raw)),

  # run linear and quadratic model
  tar_target(
    name = trait_vascular_model,
    command = run_trait_model(dat = ind_traits |>
                                filter(Functional_group == "vascular"),
                              group = c("trait_trans", "Taxon"),
                              response = Value,
                              continous_predictor = Elevation_m) |>
      pivot_longer(cols = -c(Taxon, trait_trans, data),
                   names_sep = "_",
                   names_to = c(".value", "names"))
  ),

  # select best model
  tar_target(
    name = vascular_model,
    command = trait_vascular_model |>
      # remove models that have singular fit
      filter(singular == FALSE) |>
      filter(aic == min(aic))
  ),

  # Likelihood ratio test
  tar_target(
    name = vascular_lrt,
    command = likelihood_ratio_test_vasc(ind_traits)
  ),

  # Produce model output and prediction
  tar_target(
    name = vascular_model_output,
    command = model_output_prediction(vascular_model) |>
      # add LRT text
      mutate(text = case_when(
        Taxon == "salix polaris" & trait_trans %in% c("Leaf_Area_cm2_log", "N_percent") ~ "Null",
        trait_trans == "dN15_permil" ~ "NxE",
        TRUE ~ "N+E"))
  ),

  # trait table
  tar_target(
    name = vascular_model_table,
    command = make_vascular_table(vascular_model_output)
  ),


  # BRYOPHYTES
  # run linear and quadratic model
  tar_target(
    name = trait_bryophyte_model,
    command = run_bryophyte_model(dat = ind_traits |>
                                filter(Functional_group == "bryophyte"),
                              group = c("trait_trans", "Taxon"),
                              response = Value,
                              continous_predictor = Elevation_m) |>
      pivot_longer(cols = -c(Taxon, trait_trans, data),
                   names_sep = "_",
                   names_to = c(".value", "names"))
  ),

  # select best model
  tar_target(
    name = bryophyte_model,
    command = trait_bryophyte_model |>
      filter(aic == min(aic))
  ),

  # Produce model output and prediction
  tar_target(
    name = bryophyte_model_output,
    command = bryophyte_model |>
      # make model output and prediction
      mutate(model_output = map(mod, tidy),
             prediction = map2(.x = mod, .y = data,
                               .f = ~augment(.x, interval = "confidence") |>
                                 select(.fitted, .lower, .upper) |>
                                 bind_cols(.y))) #|>
      #unnest(prediction)
  ),

  # Produce model output and prediction
  tar_target(
    name = bryo_table,
    command = make_bryo_table(bryophyte_model_output)
  )


)

