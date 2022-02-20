#
analysis_plan <- list(

  # COMMUNITY
  # test diversity along gradients
  tar_target(
    name = diversity_analysis,
    command = {

      best <- diversity_grad %>%
        ungroup() %>%
        distinct(DiversityIndex) %>%
        mutate(best_model = c("G", "G", "Null", "E"))

      # only richness and diverity
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

      diversity <- nest %>%
        mutate(mod = map(data, ~lmer(Value ~ Gradient + (1|GS), data = .x)),
               result = map(mod, tidy)) %>%
          unnest(result) %>%
        filter(effect == "fixed") %>%
        select(DiversityIndex, term:statistic) %>%
        left_join(best, by = "DiversityIndex") %>%
        left_join(r_square, by = "DiversityIndex") %>%
        select(Index = DiversityIndex, "Best model" = best_model, Estimate = estimate, "Standard error" = std.error, "t-value" = statistic, "Marginal R2" = Rm, "Conditional R2" = Rc)

      return(diversity)

    }),

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

  # make species ordination
  tar_target(
    name = sp_ordination,
    command = make_ordination(comm_raw)
  ),

  # test vascular plants along gradients
  tar_target(
    name = trait_analysis,
    command = {
      trait_mean %>%
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
        mutate(model = c("", "G", "GxE", "", "GxE", "G+E", "G", "", "", "E", "GxE", ""))
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

  # model output
  tar_target(
    name = ind_traits_output,
    command = run_ind_models(ind_traits)),

  tar_target(
    name = ind_species_figure,
    command = make_ind_sp_plot(ind_traits)),



  ### ITV
  tar_target(
    name = variation_split_exp,
    command = Intra_vs_Inter(traits_raw, trait_mean)
  ),

  tar_target(
    name = variation_split,
    command = Intra_vs_Inter_var_split(variation_split_exp)
  )

  #make_intra_vs_inter_figure(var_split_exp, var_split)


)

#variation_split_exp <- Intra_vs_Inter(traits_raw, trait_mean)
#variation_split <- Intra_vs_Inter_var_split(variation_split_exp)

