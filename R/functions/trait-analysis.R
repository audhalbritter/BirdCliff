### TRAIT ANALYSIS

# run linear and quadratic model
run_trait_model <- function(dat, group, response, continous_predictor){

  dat |>
    rename(.response = {{response}},
           .continous_predictor = {{continous_predictor}}) |>
    group_by(across(all_of({{group}})))|>
    nest() |>
    mutate(mod_linear = map(data, ~ safely(lmer)(.response ~  Gradient * .continous_predictor + (1|GS), data = .x)$result),
           singular_linear = map_lgl(mod_linear, isSingular),
           mod_quadratic = map(data, ~ safely(lmer)(.response ~  Gradient * poly(.continous_predictor, 2) + (1|GS), data = .x)$result),
           singular_quadratic = map_lgl(mod_quadratic, isSingular),
           aic_linear = map(.x = mod_linear, .f = ~ safely(AIC)(.x)$result),
           aic_linear = as.numeric(aic_linear),
           aic_quadratic = map(.x = mod_quadratic, .f = ~ safely(AIC)(.x)$result),
           aic_quadratic = as.numeric(aic_quadratic))
}


# likelihood ratio test
likelihood_ratio_test <- function(trait_mean){

  trait_mean |>
    group_by(trait_trans) |>
    nest() |>
    mutate(LTR = map(data, ~{

      if (trait_trans != "SLA_cm2_g") {
        # quadratic model
        ExN = lmer(mean ~  Gradient * poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = .x)
        EplusN = lmer(mean ~  Gradient + poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = .x)
        N = lmer(mean ~  Gradient + (1|GS), REML=FALSE, data = .x)
        E = lmer(mean ~  poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = .x)
        Null = lmer(mean ~  1 + (1|GS), REML=FALSE, data = .x)
        # lr test
        test = anova(ExN, EplusN, N, E, Null)

      } else {
        # linear model
        ExN = lmer(mean ~  Gradient * Elevation_m + (1|GS), REML=FALSE, data = .x)
        EplusN = lmer(mean ~  Gradient + Elevation_m + (1|GS), REML=FALSE, data = .x)
        N = lmer(mean ~  Gradient + (1|GS), REML=FALSE, data = .x)
        E = lmer(mean ~  Elevation_m + (1|GS), REML=FALSE, data = .x)
        Null = lmer(mean ~  1 + (1|GS), REML=FALSE, data = .x)
        # lr test
        test = anova(ExN, EplusN, N, E, Null)
      }
    })) |>
      unnest(LTR)

}

#singular fit issue: LDMC, Dry mass, LA, dC13,
# dat <- trait_mean |> #distinct(trait_trans)
#   filter(trait_trans == "N_percent")
# ExN = lmer(mean ~  Gradient * poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = dat)
# EplusN = lmer(mean ~  Gradient + poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = dat)
# N = lmer(mean ~  Gradient + (1|GS), REML=FALSE, data = dat)
# E = lmer(mean ~  poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = dat)
# Null = lmer(mean ~  1 + (1|GS), REML=FALSE, data = dat)
# test = anova(ExN, EplusN, N, E, Null)



# make model selection to find
# make_trait_model_selection <- function(trait_mean){
#
#   model.sel <- trait_mean %>%
#     group_by(trait_trans) %>%
#     nest(data = -c(trait_trans)) %>%
#     mutate(model.set = map(data, ~{
#       mod <- lmer(mean ~  Gradient * Elevation_m + (1|GS), REML = FALSE, na.action = "na.fail", data = .x)
#       model.set = dredge(mod, rank = "AICc", extra = "R^2")
#     })) %>%
#     unnest(model.set)
#
#   return(model.sel)
#
# }


# make prediction for lmer
lmer_prediction <- function(dat, fit){

  newdat <- dat %>%
    select(Gradient, .continous_predictor)

  newdat$.response <- predict(fit, newdat, re.form = NA)

  mm <- model.matrix(terms(fit), newdat)

  prediction <- newdat %>%
    mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
           tvar1 = pvar1 + VarCorr(fit)$GS[1],  ## must be adapted for more complex models
           cmult = 1.96) %>%
    mutate(plo = .response - cmult*sqrt(pvar1),
           phi = .response + cmult*sqrt(pvar1),
           tlo = .response - cmult*sqrt(tvar1),
           thi = .response + cmult*sqrt(tvar1))

  return(prediction)
}

# produce model output, r2 and prediction
model_output_prediction <- function(model){

  model |>
    # make model output and prediction
    mutate(model_output = map(mod, tidy),
           r = map(mod, r.squaredGLMM),
           r = map(r, as.numeric),
           prediction = map2(.x = data, .y = mod, .f = ~ safely(lmer_prediction)(.x, .y)$result))

}



# trait output: model output and r2
make_trait_table <- function(community_model_output){

  community_model_output |>
    select(-data, -mod, -singular, -aic, -prediction, -r) |>
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
    left_join(community_model_output |>
                select(-data, -mod, -singular, -aic, -prediction, -model_output) |>
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
    select(Class = class, Trait, Term = term, Estimate = estimate, "Std error" = std.error, "t-value" = "statistic", "Marginal R2" = Rm, "Conditional R2" = Rc) %>%
    arrange(Trait) %>%
    write_csv(., file = "output/Community_trait_regression_output.csv")

}

# Run dN15 analysis

# make model
run_dN15_model <- function(dat, group, response, continous_predictor){

  dat |>
    rename(.response = {{response}},
           .continous_predictor = {{continous_predictor}}) |>
    group_by(across(all_of({{group}})))|>
    nest() |>
    mutate(mod_linear = map(data, ~ safely(lmer)(.response ~  Gradient * .continous_predictor + (1|GS), data = .x)$result),
           singular_linear = map_lgl(mod_linear, isSingular))
}


# # likelihood ratio test
# likelihood_ratio_test_variance <- function(trait_mean){
#
#   trait_mean |>
#     group_by(trait_trans) |>
#     nest() |>
#     mutate(LTR = map(data, ~{
#         # quadratic model
#         ExN = lmer(var ~  Gradient * poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = .x)
#         EplusN = lmer(var ~  Gradient + poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = .x)
#         N = lmer(var ~  Gradient + (1|GS), REML=FALSE, data = .x)
#         E = lmer(var ~  poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = .x)
#         Null = lmer(var ~  1 + (1|GS), REML=FALSE, data = .x)
#         # lr test
#         test = anova(ExN, EplusN, N, E, Null)
#     })) |>
#     unnest(LTR)
#
# }



# # trait output: model output and r2
# make_trait_variance_table <- function(community_variance_output){
#
#   community_variance_output |>
#     select(-data, -mod, -singular, -aic, -prediction, -r) |>
#     unnest(model_output) |>
#     fancy_trait_name_dictionary()  |>
#     filter(effect == "fixed")  |>
#     ungroup() |>
#     select(Trait = trait_fancy, text, term, estimate:statistic, class) %>%
#     mutate(term = recode(term,
#                          "(Intercept)" = "Intercept",
#                          "GradientB" = "N",
#                          ".continous_predictor" = "E",
#                          "GradientB:.continous_predictor" = "NxE" ,
#                          "poly(.continous_predictor, 2)1" = "E",
#                          "poly(.continous_predictor, 2)2" = "E2",
#                          "GradientB:poly(.continous_predictor, 2)1" = "NxE",
#                          "GradientB:poly(.continous_predictor, 2)2" = "NxE2")) %>%
#     # join R squared
#     left_join(community_variance_output |>
#                 select(-data, -mod, -singular, -aic, -prediction, -model_output) |>
#                 unnest_wider(col = r, names_sep = "_") |>
#                 rename(Rm = r_1, Rc = r_2) |>
#                 ungroup() %>%
#                 fancy_trait_name_dictionary(.) %>%
#                 select(trait_fancy, Rm, Rc),
#               by = c("Trait" = "trait_fancy")) %>%
#     mutate(estimate = round(estimate, digits = 2),
#            std.error = round(std.error, digits = 2),
#            statistic = round(statistic, digits = 2),
#            Rm = round(Rm, digits = 2),
#            Rc = round(Rc, digits = 2)) %>%
#     select(Class = class, Trait, Term = term, Estimate = estimate, "Std error" = std.error, "t-value" = "statistic", "Marginal R2" = Rm, "Conditional R2" = Rc) %>%
#     arrange(Trait) %>%
#     write_csv(., file = "output/Community_trait_variance_regression_output.csv")
#
# }

