### TRAIT ANALYSIS

# test linear vs polynomial model
lin_vs_poly_model <- function(trait_mean){

  trait_mean |>
    group_by(trait_trans) |>
    nest() |>
    mutate(lin_mod = map(data, ~ safely(lmer)(mean ~  Gradient * Elevation_m + (1|GS), data = .x)$result),
           poly_mod = map(data, ~ safely(lmer)(mean ~  Gradient * poly(Elevation_m, 2) + (1|GS), data = .x)$result),
           mod_aic = map2(.x = lin_mod, .y = poly_mod, .f = AIC)) |>
    unnest(col = mod_aic)
}

likelihood_ratio_test <- function(trait_mean){

  trait_mean |>
    filter(trait_trans != "SLA_cm2_g") |>
    group_by(trait_trans) |>
    nest() |>
    mutate(LTR = map(data, ~{
      # models
      ExN = lmer(mean ~  Gradient * poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = .x)
      EplusN = lmer(mean ~  Gradient + poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = .x)
      N = lmer(mean ~  Gradient + (1|GS), REML=FALSE, data = .x)
      E = lmer(mean ~  poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = .x)
      Null = lmer(mean ~  1 + (1|GS), REML=FALSE, data = .x)
      # ltr
      test = anova(ExN, EplusN, N, E, Null)
      })) |>
    unnest(LTR) |> print(n = Inf)

}

#singular fit issue: LDMC, Dry mass, LA, dC13,
# dat <- trait_mean |>
#   filter(trait_trans == "dN15_permil")
# ExN = lmer(mean ~  Gradient * poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = dat)
# EplusN = lmer(mean ~  Gradient + poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = dat)
# N = lmer(mean ~  Gradient + (1|GS), REML=FALSE, data = dat)
# E = lmer(mean ~  poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = dat)
# Null = lmer(mean ~  1 + (1|GS), REML=FALSE, data = dat)
#test = anova(ExN, EplusN, N, E, Null)


# run full model
# trait_mean %>%
#   group_by(trait_trans) %>%
#   nest(data = -c(trait_trans)) %>%
#   mutate(model = map(data, ~ safely(lmer)(mean ~  Gradient * Elevation_m + (1|GS), data = .)$result),
#          # add performance
#          tidy_result = map(model, tidy)) %>%
#   unnest(tidy_result)

# make model selection to find
make_trait_model_selection <- function(trait_mean){

  model.sel <- trait_mean %>%
    group_by(trait_trans) %>%
    nest(data = -c(trait_trans)) %>%
    mutate(model.set = map(data, ~{
      mod <- lmer(mean ~  Gradient * Elevation_m + (1|GS), REML = FALSE, na.action = "na.fail", data = .x)
      model.set = dredge(mod, rank = "AICc", extra = "R^2")
    })) %>%
    unnest(model.set)

  return(model.sel)

}


# make prediction for lmer
lmer_prediction <- function(dat, fit){

  newdat <- dat %>%
    select(Gradient, Elevation_m) %>%
    mutate(mean = 0)
  newdat$mean <- predict(fit, newdat, re.form = NA)

  mm <- model.matrix(terms(fit), newdat)

  prediction <- newdat %>%
    mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
           tvar1 = pvar1 + VarCorr(fit)$GS[1],  ## must be adapted for more complex models
           cmult = 1.96) %>%
    mutate(plo = mean - cmult*sqrt(pvar1),
           phi = mean + cmult*sqrt(pvar1),
           tlo = mean - cmult*sqrt(tvar1),
           thi = mean + cmult*sqrt(tvar1))

  return(prediction)
}


run_best_model_for_prediction <- function(trait_mean){

  #best_model
  out <- fancy_trait_name_dictionary(trait_mean) |>
    # add best model (info from model selection)
    mutate(best_model = case_when(trait_trans == "Plant_Height_cm_log" ~ "mean ~ Elevation_m + (1|GS)",
                                  trait_trans %in% c("CN_ratio", "N_percent", "dC13_permil") ~ "mean ~ Gradient + (1|GS)",
                                  trait_trans %in% c("Leaf_Area_cm2_log", "LDMC") ~ "mean ~ Gradient + Elevation_m + (1|GS)",
                                  trait_trans %in% c("SLA_cm2_g", "dN15_permil") ~ "mean ~ Gradient * Elevation_m + (1|GS)",
                                  TRUE ~ "mean ~ 1 + (1|GS)"),
           # text for model selection
           # text = case_when(trait_trans == "Plant_Height_cm_log" ~ "E",
           #                        trait_trans %in% c("CN_ratio", "N_percent", "dC13_permil") ~ "N",
           #                        trait_trans %in% c("Leaf_Area_cm2_log", "LDMC") ~ "N+E",
           #                        trait_trans %in% c("SLA_cm2_g", "dN15_permil") ~ "NxE",
           #                        TRUE ~ "NULL")
           # text for full LRT
           text = case_when(trait_trans %in% c("NP_ratio", "dC13_permil") ~ "Null",
                            trait_trans %in% c("CN_ratio", "N_percent") ~ "N",
                            TRUE ~ "NxE")
           ) |>
    # run model
    group_by(trait_trans, best_model, text) |>
    nest() |>
    mutate(
      # run best model
      #model = map(data, ~ safely(lmer)(formula = best_model, data = .)$result),
      # run full model
      model = map(data, ~ safely(lmer)(mean ~  Gradient * poly(Elevation_m, 2) + (1|GS), data = .x)$result),
           # model output
           model_output = map(model, tidy),
           r = map(model, r.squaredGLMM),
           r = map(r, as.numeric),
           # make model prediction
           prediction = map2(.x = data, .y = model, ~ safely(lmer_prediction)(.x, .y)$result))

  return(out)

}



# trait output: model output and r2
make_trait_output <- function(community_trait_model){

  community_trait_model |>
    select(-data, -model, -prediction, -r) |>
    unnest(model_output) %>%
    fancy_trait_name_dictionary(.) %>%
    filter(effect == "fixed") %>%
    ungroup() |>
    select(Trait = trait_fancy, text, term, estimate:statistic, class) %>%
    mutate(term = recode(term, "(Intercept)" = "Intercept", "Elevation_m" = "E", "GradientC" = "G", "GradientC:Elevation_m" = "GxE")) %>%
    # join R squared
    left_join(community_trait_model |>
                select(-data, -model, -prediction, -model_output) |>
                unnest_wider(col = r) |>
                rename(Rm = ...1, Rc = ...2) |>
                ungroup() %>%
                fancy_trait_name_dictionary(.) %>%
                select(trait_fancy, Rm, Rc),
              by = c("Trait" = "trait_fancy")) %>%
    mutate(estimate = round(estimate, digits = 2),
           std.error = round(std.error, digits = 2),
           statistic = round(statistic, digits = 2),
           Rm = round(Rm, digits = 2),
           Rc = round(Rc, digits = 2)) %>%
    select(Class = class, Trait, Term = term, Estimate = estimate, "Std error" = std.error, "t-value" = "statistic", "Marginal R2" = Rm, "Conditional R2" = Rc) %>%
    arrange(Trait) %>%
    write_csv(., file = "output/Community_trait_regression_output.csv")

}
