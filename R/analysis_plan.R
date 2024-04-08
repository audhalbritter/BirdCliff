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
    command = tidy(comm_pca[[4]])
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


  ### SOIL dN15 TRAIT ANALYSIS

  # run model
  tar_target(
    name = dN15_trait_model,
    command = run_trait_model(dat = trait_mean |>
                               filter(!is.na(d15n)),
                              group = "trait_trans",
                              response = mean,
                              continous_predictor = d15n) |>
        pivot_longer(cols = -c(trait_trans, data),
                     names_sep = "_",
                     names_to = c(".value", "names"))
  ),

  # select best model
  tar_target(
    name = dN15_model,
    command = dN15_trait_model |>
      # remove models that have singular fit
      filter(singular == FALSE) |>
      filter(aic == min(aic))
  ),

  # Likelihood ratio test
  tar_target(
    name = dN15_lrt,
    command = trait_mean |>
      filter(!is.na(d15n)) |>
      filter(trait_trans == "LDMC") |>
      group_by(trait_trans) |>
      nest() |>
      mutate(LTR = map(data, ~{
        ExN = lmer(mean ~  Gradient * poly(d15n, 2) + (1|GS), REML=FALSE, data = .x)
        EplusN = lmer(mean ~  Gradient + poly(d15n, 2) + (1|GS), REML=FALSE, data = .x)
        N = lmer(mean ~  Gradient + (1|GS), REML=FALSE, data = .x)
        E = lmer(mean ~  poly(d15n, 2) + (1|GS), REML=FALSE, data = .x)
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
      mutate(text = case_when(trait_trans %in% c("NP_ratio") ~ "NxδN15",
                              trait_trans %in% c("dN15_permil") ~ "N+δN15",
                              trait_trans %in% c("N_percent", "CN_ratio", "dC13_permil") ~ "N",
                              trait_trans %in% c("Plant_Height_cm_log", "Dry_Mass_g_log") ~ "Null",
                              TRUE ~ "δN15"))
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
                           ".continous_predictor" = "E",
                           "GradientB:.continous_predictor" = "NxE" ,
                           "poly(.continous_predictor, 2)1" = "E",
                           "poly(.continous_predictor, 2)2" = "E2",
                           "GradientB:poly(.continous_predictor, 2)1" = "NxE",
                           "GradientB:poly(.continous_predictor, 2)2" = "NxE2")) %>%
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

  # correlation between leaf and soil dN15
  tar_target(
    name = dN15_cor,
    command = trait_mean |>
      filter(trait_trans == "dN15_permil",
             !is.na(d15n)) |>
      group_by(Gradient) |>
      summarise(cor_dn15 = cor(mean, d15n))
  ),


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
        filter(class != "Soil abiotics") |>
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
    command = tidy(trait_pca[[4]])
    ),


  # test PCA1 regression
  # run linear and quadratic model
  # tar_target(
  #   name = PCA1_community_model,
  #   command = {
  #
  #     dat <- trait_pca[[1]] |>
  #       rename(Elevation_m = Mean_elevation)
  #
  #     mod_linear <-  lmer(PC1 ~  Gradient * Elevation_m + (1|GS), data = dat)
  #     mod_quadratic <-  lmer(PC1 ~  Gradient * poly(Elevation_m, 2) + (1|GS), data = dat)
  #     AIC(mod_linear, mod_quadratic)
  #
  #     ExN = lmer(PC1 ~  Gradient * poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = dat)
  #     EplusN = lmer(PC1 ~  Gradient + poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = dat)
  #     N = lmer(PC1 ~  Gradient + (1|GS), REML=FALSE, data = dat)
  #     E = lmer(PC1 ~  poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = dat)
  #     Null = lmer(PC1 ~  1 + (1|GS), REML=FALSE, data = dat)
  #     # lr test
  #     test = anova(ExN, EplusN, N, E, Null)
  #
  #
  #     r_squared <- r.squaredGLMM(mod_quadratic) |>
  #       as_tibble()
  #
  #     # best model is NxE
  #     tidy(mod_quadratic) |>
  #       filter(effect == "fixed") |>
  #       mutate(term = recode(term,
  #                            "(Intercept)" = "Intercept",
  #                            "GradientB" = "N",
  #                            "Elevation_m" = "E",
  #                            "GradientB:.continous_predictor" = "NxE" ,
  #                            "poly(Elevation_m, 2)1" = "E",
  #                            "poly(Elevation_m, 2)2" = "E2",
  #                            "GradientB:poly(Elevation_m, 2)1" = "NxE",
  #                            "GradientB:poly(Elevation_m, 2)2" = "NxE2")) |>
  #       bind_cols(r_squared) |>
  #       mutate(estimate = round(estimate, digits = 2),
  #              std.error = round(std.error, digits = 2),
  #              statistic = round(statistic, digits = 2),
  #              R2m = round(R2m, digits = 2),
  #              R2c = round(R2c, digits = 2)) %>%
  #       ungroup() |>
  #       select(Term = term, Estimate = estimate, "Std error" = std.error, "t-value" = "statistic", "Marginal R2" = R2m, "Conditional R2" = R2c) %>%
  #       write_csv(., file = "output/PCA_regression_output.csv")
  #
  #     }
  # ),


  ### ITV ANALYSIS
  tar_target(
    name = itv_output,
    command = make_ITV_analysis(trait_mean))

)

