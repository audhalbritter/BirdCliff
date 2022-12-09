#
analysis_plan <- list(

  # COMMUNITY
  # make species ordination
  # tar_target(
  #   name = sp_ordination,
  #   command = make_ordination(comm_raw)
  # ),
  #
  # # test species ordination
  # tar_target(
  #   name = output_sp_ordination,
  #   command = test_ordination(comm_raw)
  # ),

  # community pca
  # make species ordination
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

  # Community traits
  # run linear and quadratic model
  tar_target(
    name = trait_community_model,
    command = run_trait_model(dat = trait_mean,
                               group = "trait_trans",
                               response = mean,
                               continous_predictor = Elevation_m) |>
      pivot_longer(cols = -c(trait_trans, data),
                   #names_pattern = "(.*)(linear|quadratic)$",
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

  # PROBABLY NOT NEEDED ANYMORE!!!
  # test top site
  # tar_target(
  #   name = top_site,
  #   command = test_top_site(trait_mean)
  # ),




  # make trait ordination
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
        mutate(Gradient = recode(Gradient, Birdcliff = "Bird cliff"),
               Gradient = factor(Gradient, levels = c("Bird cliff", "Reference"))) %>%
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

  ### NEEDS TO BE MADE
  # trait table
  # tar_target(
  #   name = vascular_model_table,
  #   command = make_trait_table(vascular_model_output)
  # ),

  # vascular: run model
  # tar_target(
  #   name = ind_trait_output,
  #   command = run_vascular_model(ind_traits)
  # ),

  # vascular: run LRT
  # tar_target(
  #   name = vascular_lrt,
  #   command = likelihood_ratio_test_ind(ind_traits)
  # ),

  # run model selection
  # does not work yet!

  # vascular: model output
  # tar_target(
  #   name = ind_vascular_traits_output,
  #   command = run_vascular_plant_models(ind_traits)),


  # BRYOPHYTES
  tar_target(
    name = bryo_trait_output,
    command = make_bryo_trait_model(ind_traits)),


  ### ITV
  tar_target(
    name = itv_output,
    command = make_ITV_analysis(trait_mean))

)

