#
analysis_plan <- list(

  # COMMUNITY
  # test diversity along gradients
  tar_target(
    name = diversity_analysis,
    command = {
      diversity_grad %>%
        mutate(GS = paste0(Gradient, Site)) %>%
        group_by(DiversityIndex) %>%
        nest(data = -c(DiversityIndex)) %>%
        mutate(mod = map(data, ~lmer(Value ~ Gradient * Elevation_m + (1|GS), data = .x)),
               result = map(mod, tidy)) %>%
          unnest(result)
    }),


  # make species ordination
  tar_target(
    name = fNMDS,
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

  # INDIVIDUAL TRAITS

  tar_target(
    name = all_traits,
    command = combine_traits(traits_raw, bryo_traits_raw)),

  tar_target(
    name = ind_species_figure,
    command = make_ind_sp_plot(all_traits)),

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

