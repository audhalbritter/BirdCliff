#
analysis_plan <- list(

  # COMMUNITY
  #test diversity along gradients
  tar_target(
    name = diversity_analysis,
    command = {

      best <- diversity_grad %>%
        ungroup() %>%
        distinct(DiversityIndex) %>%
        mutate(best_model = c("G", "G", "Null", "Null"))

      # GRADIENT MODEL only richness and diverity
      nest <- diversity_grad %>%
        mutate(GS = paste0(Gradient, Site)) %>%
        filter(DiversityIndex %in% c("Richness", "Diversity")) %>%
        group_by(DiversityIndex) %>%
        nest(data = -c(DiversityIndex))

      r_square <- nest %>%
        mutate(r = map(data, ~{
          mod <- lmer(Value ~ Gradient + (1|GS), data = .x)
          r = as.numeric(r.squaredGLMM(mod))
        })) %>%
        unnest_wider(col = r) %>%
        select(DiversityIndex, "Rm" = "...1", "Rc" = "...2")

      estimate <- nest %>%
        mutate(mod = map(data, ~lmer(Value ~ Gradient + (1|GS), data = .x)),
               result = map(mod, tidy)) %>%
          unnest(result)

      # NULL MODEL evenness and sumAbundance
      nest_2 <- diversity_grad %>%
        mutate(GS = paste0(Gradient, Site)) %>%
        filter(DiversityIndex %in% c("Evenness", "sumAbundance")) %>%
        group_by(DiversityIndex) %>%
        nest(data = -c(DiversityIndex))

      r_square_2 <- nest_2 %>%
        mutate(r = map(data, ~{
          mod <- lmer(Value ~ 1 + (1|GS), data = .x)
          r = as.numeric(r.squaredGLMM(mod))
        })) %>%
        unnest_wider(col = r) %>%
        select(DiversityIndex, "Rm" = "...1", "Rc" = "...2")

      estimate_2 <- nest_2 %>%
        mutate(mod = map(data, ~lmer(Value ~ 1 + (1|GS), data = .x)),
               result = map(mod, tidy)) %>%
        unnest(result)

      diversity_output <- bind_rows(estimate, estimate_2) %>%
        select(-data, -mod) %>%
        filter(effect == "fixed") %>%
        left_join(best, by = "DiversityIndex") %>%
        left_join(bind_rows(r_square, r_square_2), by = "DiversityIndex") %>%
        select(Index = DiversityIndex, "Best model" = best_model, term, Estimate = estimate, "Standard error" = std.error, "t-value" = statistic, "Marginal R2" = Rm, "Conditional R2" = Rc)

      return(diversity_output)

    }),

  # test best diversity model
  tar_target(
    name = div_best_model,
    command = {
      diversity_grad %>%
        mutate(GS = paste0(Gradient, Site)) %>%
        group_by(DiversityIndex) %>%
        nest(data = -c(DiversityIndex)) %>%
        mutate(model.set = map(data, ~{
          mod <- lmer(Value ~  Gradient * Elevation_m + (1|GS), REML = FALSE, na.action = "na.fail", data = .x)
          model.set = dredge(mod, rank = "AICc", extra = "R^2")
          })) %>%
        unnest(model.set)
    }),

  # check the nr of dimensions for NMDS
  tar_target(
    name = stress_plot,
    command = check_dimensions_NMDS(comm_raw)
  ),

  # make species ordination
  tar_target(
    name = sp_ordination,
    command = make_ordination(comm_raw)
  ),

  # test species ordination
  tar_target(
    name = output_sp_ordination,
    command = test_ordination(comm_raw)
  ),

  # test vascular plants along gradients
  tar_target(
    name = trait_analysis,
    command = {
      trait_mean %>%
        filter(trait_trans != "dC13_permil") %>%
        group_by(trait_trans) %>%
        nest(data = -c(trait_trans)) %>%
        mutate(mod = map(data, ~lmer(mean ~ Gradient * Elevation_m + (1|Site), data = .x)),
               result = map(mod, tidy)) %>%
        unnest(result)
    }),


  # TRAITS
  # model selection
  tar_target(
    name = model_selection,
    command = make_trait_model_selection(trait_mean)
  ),

  # Best model: likelihoot ratio test result
  tar_target(
    name = best_model,
    command = {
      trait_mean %>%
        distinct(trait_trans) %>%
        filter(trait_trans != "dC13_permil") %>%
        mutate(model = c("", "G", "GxE", "", "G+E", "G+E", "G", "", "", "E", "G+E", ""))
    }),

  tar_target(
    name = model_output,
    command = make_trait_output(trait_mean)
  ),


  # FIGURE 2a: trait change along gradients
  tar_target(
    name = trait_plot,
    command = {
      make_trait_figure(trait_mean)
    }),


  # make trait ordination
  tar_target(
    name = trait_pca_B,
    command = make_trait_pca(trait_mean %>% filter(Gradient == "B"))
  ),

  tar_target(
    name = trait_pca_C,
    command = make_trait_pca(trait_mean %>% filter(Gradient == "C"))
  ),

  # tar_target(
  #   name = trait_pca,
  #   command = make_trait_pca(trait_mean)
  # ),

  # FIGURE 2B: trait ordination
  tar_target(
    name = trait_ordination_plot,
    command = make_trait_pca_plot(trait_pca_B, trait_pca_C)
  ),

  tar_target(
    name = trait_ord_expl_var,
    command = {

      bind_cols(
        bird = vegan::eigenvals(trait_pca_B[[3]])/sum(vegan::eigenvals(trait_pca_B[[3]])) * 100,
        reference = vegan::eigenvals(trait_pca_C[[3]])/sum(vegan::eigenvals(trait_pca_C[[3]])) * 100)

    }),

  # INDIVIDUAL TRAITS (vascular plants)

  # combine data
  tar_target(
    name = ind_traits,
    command = combine_traits(traits_raw, bryo_traits_raw)),

  # run model selection
  # does not work yet!

  # vascular: model output
  tar_target(
    name = ind_traits_output,
    command = run_ind_models(ind_traits)),

  tar_target(
    name = ind_species_figure,
    command = make_ind_sp_plot(ind_traits)),

  # bryophytes
  tar_target(
    name = bryo_trait_output,
    command = make_ind_sp_plot(ind_traits)),


  ### ITV
  tar_target(
    name = itv_output,
    command = make_ITV_analysis(trait_mean)),

  tar_target(
    name = ITV_plot,
    command = make_ITV_plot(itv_output)
  )

)

