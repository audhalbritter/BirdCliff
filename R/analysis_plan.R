#
analysis_plan <- list(
  # do something interesting

  # test diversity along gradients
  tar_target(
    name = diversity_analysis,
    command = {
      diversity_grad %>%
        group_by(DiversityIndex) %>%
        nest(data = -c(DiversityIndex)) %>%
        mutate(mod = map(data, ~lm(Value ~ Gradient * Elevation_m, data = .x)),
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

  # make trait ordination
  tar_target(
    name = trait_pca_B,
    command = make_trait_pca(trait_mean %>% filter(Gradient == "B"))
  ),

  tar_target(
    name = trait_pca_C,
    command = make_trait_pca(trait_mean %>% filter(Gradient == "C"))
  ),

  ## Intra vs. Inter
  tar_target(
    name = variation_split_exp,
    command = Intra_vs_Inter(trait_mean)
  ),

  tar_target(
    name = variation_split,
    command = Intra_vs_Inter_var_split(variation_split_exp)
  )


  # some preliminary rough analysis
  # trait_mean %>%
  #   group_by(trait_trans) %>%
  #   nest(data = -c(trait_trans)) %>%
  #   mutate(mod = map(data, ~lm(mean ~ Gradient * Site, data = .x)),
  #          result = map(mod, tidy)) %>%
  #   unnest(result) %>%
  #   filter(p.value < 0.05) %>%
  #   select(trait_trans, term:p.value) %>% View()
  #
  #
  # diversity_grad %>%
  #   group_by(DiversityIndex) %>%
  #   nest(data = -c(DiversityIndex)) %>%
  #   mutate(mod = map(data, ~lm(Value ~ Gradient * Site, data = .x)),
  #          result = map(mod, tidy)) %>%
  #   unnest(result) %>%
  #   filter(p.value < 0.05)

)


# Check model fit for diversity Indices
# library(performance)
# dd <- diversity_grad %>%
#   filter(DiversityIndex == "sumAbundance")
# fit <- lm(Value ~ Gradient * Elevation_m, data = dd)
# check_model(fit)

