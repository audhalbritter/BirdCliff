#
analysis_plan <- list(

  # ENVIRONMENT DATA (Soil temperature, moisture, C and N)
  tar_target(
    name = environment_result,
    command = environment |>
      group_by(Variable, Gradient) |>
      summarise(mean = mean(Value),
                se = sd(Value)/sqrt(n()))
  ),


  # VASCULAR PLANT COMMUNITY
  # Taxonomic PCA
  tar_target(
    name = comm_pca,
    command = make_community_pca(comm_raw)
  ),

  tar_target(
    name = adonis_comm_output,
    command = tidy(comm_pca[[4]]$aov.tab)
  ),


  # PLANT FUNCTIONAL TRAIT REGRESSIONS
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

  # TRAIT VS MICROCLIMATE
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


  # dN15 trait analysis

  # prep data
  tar_target(
    name = trait_mean_dN15,
    command = trait_mean |>
      select(Gradient:mean, GS) |>
      pivot_wider(names_from = trait_trans, values_from = mean) |>
      pivot_longer(cols = c(Plant_Height_cm_log:dC13_permil), names_to = "trait_trans", values_to = "mean")
  ),

  # run model
  tar_target(
    name = dN15_model,
    command = run_dN15_model(dat = trait_mean_dN15,
                              group = "trait_trans",
                              response = mean,
                              continous_predictor = dN15_permil) |>
        pivot_longer(cols = -c(trait_trans, data),
                     names_sep = "_",
                     names_to = c(".value", "names")) |>
        # remove models that have singular fit
        filter(singular == FALSE)
  ),

  # Likelihood ratio test
  tar_target(
    name = dN15_lrt,
    command = trait_mean_dN15 |>
      group_by(trait_trans) |>
      nest() |>
      mutate(LTR = map(data, ~{
        ExN = lmer(mean ~  Gradient * dN15_permil + (1|GS), REML=FALSE, data = .x)
        EplusN = lmer(mean ~  Gradient + dN15_permil + (1|GS), REML=FALSE, data = .x)
        N = lmer(mean ~  Gradient + (1|GS), REML=FALSE, data = .x)
        E = lmer(mean ~  dN15_permil + (1|GS), REML=FALSE, data = .x)
        Null = lmer(mean ~  1 + (1|GS), REML=FALSE, data = .x)
        # lr test
        test = anova(ExN, EplusN, N, E, Null)
      })) |>
      unnest(LTR)

  ),

  # Produce model output and prediction
  tar_target(
    name = dN15_model_output,
    command = model_output_prediction(dN15_model) |>
      # add LRT text
      mutate(text = case_when(trait_trans %in% c("Plant_Height_cm_log", "C_percent") ~ "NxδN15",
                              trait_trans %in% c("Dry_Mass_g_log", "Leaf_Area_cm2_log", "Thickness_mm_log", "LDMC") ~ "N+δN15",
                              trait_trans %in% c("N_percent", "CN_ratio") ~ "δN15",
                              TRUE ~ "Null"))
  ),

  # regression model output
  tar_target(
    name = dN15_model_table,
    command = dN15_model_output |>
      select(-data, -mod, -singular, -prediction, -r) |>
      unnest(model_output) %>%
      fancy_trait_name_dictionary(.) %>%
      filter(effect == "fixed") %>%
      ungroup() |>
      select(Trait = trait_fancy, text, term, estimate:statistic, class) %>%
      mutate(term = recode(term,
                           "(Intercept)" = "Intercept",
                           "GradientB" = "N",
                           ".continous_predictor" = "δN15",
                           "GradientB:.continous_predictor" = "NxδN15")) %>%
      # join R squared
      left_join(dN15_model_output |>
                  select(-data, -mod, -singular, -prediction, -model_output) |>
                  unnest_wider(col = r, names_sep = "_") |>
                  rename(Rm = r_1, Rc = r_2) |>
                  ungroup() %>%
                  fancy_trait_name_dictionary(.) %>%
                  select(trait_fancy, Rm, Rc),
                by = c("Trait" = "trait_fancy")) %>%
      mutate(estimate = round(estimate, digits = 2),
             std.error = round(std.error, digits = 2),
             statistic = round(statistic, digits = 2),
             Rm = round(Rm, digits = 2),
             Rc = round(Rc, digits = 2)) %>%
      select("Trait category" = class, Trait, Term = term, Estimate = estimate, "Std error" = std.error, "t-value" = "statistic", "Marginal R2" = Rm, "Conditional R2" = Rc) %>%
      write_csv(., file = "output/dN15_trait_regression_output.csv")
  ),



  # # Community trait variance
  # # run linear and quadratic model
  # # NP ratio, both models are singular fit
  # tar_target(
  #   name = trait_variance_model,
  #   command = run_trait_model(dat = trait_mean,
  #                             group = "trait_trans",
  #                             response = var,
  #                             continous_predictor = Elevation_m) |>
  #     pivot_longer(cols = -c(trait_trans, data),
  #                  names_sep = "_",
  #                  names_to = c(".value", "names"))
  # ),
  #
  # # select best model
  # tar_target(
  #   name = community_variance_model,
  #   command = trait_variance_model |>
  #     # remove models that have singular fit
  #     filter(singular == FALSE) |>
  #     filter(aic == min(aic))
  # ),
  #
  # # LRT for variance
  # tar_target(
  #   name = ltr_variance,
  #   command = likelihood_ratio_test_variance(trait_mean |>
  #                                              filter(trait_trans != "NP_ratio"))
  # ),
  #
  # # Produce model output and prediction
  # tar_target(
  #   name = community_variance_output,
  #   command = model_output_prediction(community_variance_model) |>
  #     # add LRT text
  #     mutate(text = case_when(trait_trans %in% c("SLA_cm2_g", "dC13_permil") ~ "Null",
  #                             trait_trans %in% c("Plant_Height_cm_log") ~ "E",
  #                             TRUE ~ "NxE"))
  # ),
  #
  # # trait table
  # tar_target(
  #   name = trait_variance_table,
  #   command = make_trait_variance_table(community_variance_output)
  # ),


  # MULTIVARIATE TRAIT ANALYSIS
  # TRAIT ORDINATION (PCA)
  # both localities
  tar_target(
    name = trait_pca,
    command = make_trait_pca(trait_mean)
  ),

  # pca output
  tar_target(
    name = trait_pca_output,
    command = {

      trait_pca[[2]] %>%
        filter(class != "Environment") |>
        select(Trait = trait_fancy, PC1:PC4) %>%
        mutate(PC1 = round(PC1, digits = 2),
               PC2 = round(PC2, digits = 2),
               PC3 = round(PC3, digits = 2),
               PC4 = round(PC4, digits = 2)) %>%
        arrange(-PC1) %>%
        # order traits
        write_csv(., file = "output/Loadings_trait_PCA.csv")

    }
  ),

  tar_target(
    name = trait_ord_expl_var,
    command = {

      vegan::eigenvals(trait_pca[[3]])/sum(vegan::eigenvals(trait_pca[[3]])) * 100

    }),

  tar_target(
    name = adonis_trait_output,
    command = tidy(trait_pca[[4]]$aov.tab)
    ),


  # test PCA1 regression
  # run linear and quadratic model
  tar_target(
    name = PCA1_community_model,
    command = {

      dat <- trait_pca[[1]] |>
        rename(Elevation_m = Mean_elevation)

      mod_linear <-  lmer(PC1 ~  Gradient * Elevation_m + (1|GS), data = dat)
      mod_quadratic <-  lmer(PC1 ~  Gradient * poly(Elevation_m, 2) + (1|GS), data = dat)
      AIC(mod_linear, mod_quadratic)

      ExN = lmer(PC1 ~  Gradient * poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = dat)
      EplusN = lmer(PC1 ~  Gradient + poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = dat)
      N = lmer(PC1 ~  Gradient + (1|GS), REML=FALSE, data = dat)
      E = lmer(PC1 ~  poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = dat)
      Null = lmer(PC1 ~  1 + (1|GS), REML=FALSE, data = dat)
      # lr test
      test = anova(ExN, EplusN, N, E, Null)


      r_squared <- r.squaredGLMM(mod_quadratic) |>
        as_tibble()

      # best model is NxE
      tidy(mod_quadratic) |>
        filter(effect == "fixed") |>
        mutate(term = recode(term,
                             "(Intercept)" = "Intercept",
                             "GradientB" = "N",
                             "Elevation_m" = "E",
                             "GradientB:.continous_predictor" = "NxE" ,
                             "poly(Elevation_m, 2)1" = "E",
                             "poly(Elevation_m, 2)2" = "E2",
                             "GradientB:poly(Elevation_m, 2)1" = "NxE",
                             "GradientB:poly(Elevation_m, 2)2" = "NxE2")) |>
        bind_cols(r_squared) |>
        mutate(estimate = round(estimate, digits = 2),
               std.error = round(std.error, digits = 2),
               statistic = round(statistic, digits = 2),
               R2m = round(R2m, digits = 2),
               R2c = round(R2c, digits = 2)) %>%
        ungroup() |>
        select(Term = term, Estimate = estimate, "Std error" = std.error, "t-value" = "statistic", "Marginal R2" = R2m, "Conditional R2" = R2c) %>%
        write_csv(., file = "output/PCA_regression_output.csv")

      }
  ),


  ### ITV ANALYSIS
  tar_target(
    name = itv_output,
    command = make_ITV_analysis(trait_mean))


  # # BRYOPHYTES
  # # combine vascular and bryophyte data
  # tar_target(
  #   name = ind_traits,
  #   command = combine_traits(traits_raw, bryo_traits_raw)),
  #
  #
  # # BRYOPHYTES (maybe remove?)
  # # run linear and quadratic model
  # tar_target(
  #   name = trait_bryophyte_model,
  #   command = run_bryophyte_model(dat = ind_traits |>
  #                               filter(Functional_group == "bryophyte"),
  #                             group = c("trait_trans", "Taxon"),
  #                             response = Value,
  #                             continous_predictor = Elevation_m) |>
  #     pivot_longer(cols = -c(Taxon, trait_trans, data),
  #                  names_sep = "_",
  #                  names_to = c(".value", "names"))
  # ),
  #
  # # select best model
  # tar_target(
  #   name = bryophyte_model,
  #   command = trait_bryophyte_model |>
  #     filter(aic == min(aic))
  # ),
  #
  # # Produce model output and prediction
  # tar_target(
  #   name = bryophyte_model_output,
  #   command = bryophyte_model |>
  #     # make model output and prediction
  #     mutate(model_output = map(mod, tidy),
  #            prediction = map2(.x = mod, .y = data,
  #                              .f = ~augment(.x, interval = "confidence") |>
  #                                select(.fitted, .lower, .upper) |>
  #                                bind_cols(.y))) #|>
  #     #unnest(prediction)
  # ),
  #
  # # Produce model output and prediction
  # tar_target(
  #   name = bryo_table,
  #   command = make_bryo_table(bryophyte_model_output)
  # )


)

