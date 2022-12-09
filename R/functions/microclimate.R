### CLIMATE DATA

# make prediction for lmer
lmer_climate_prediction <- function(dat, fit){

  newdat <- dat %>%
    select(Gradient, Elevation_m) %>%
    mutate(Value = 0)
  newdat$Value <- predict(fit, newdat, re.form = NA)

  mm <- model.matrix(terms(fit), newdat)

  prediction <- newdat %>%
    mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
           tvar1 = pvar1 + VarCorr(fit)$GS[1],  ## must be adapted for more complex models
           cmult = 1.96) %>%
    mutate(plo = Value - cmult*sqrt(pvar1),
           phi = Value + cmult*sqrt(pvar1),
           tlo = Value - cmult*sqrt(tvar1),
           thi = Value + cmult*sqrt(tvar1)) |>
    rename(fitted = Value, Gradient_fit = Gradient, Elevation_fit = Elevation_m)

  return(prediction)
}


make_climate_analysis <- function(climate_data, coordinate){

  # quadratic model is better
  model_choice <- climate_data |>
    group_by(Variable) |>
    nest() |>
    mutate(lin_mod = map(data, ~ safely(lmer)(Value ~  Gradient * Elevation_m + (1|GS), data = .x)$result),
           poly_mod = map(data, ~ safely(lmer)(Value ~  Gradient * poly(Elevation_m, 2) + (1|GS), data = .x)$result),
           mod_aic = map2(.x = lin_mod, .y = poly_mod, .f = AIC)) |>
    unnest(col = mod_aic) |>
    mutate(Model = c("linear", "quadratic"))

  lrt <- climate_data|>
    group_by(Variable) |>
    nest() |>
    mutate(LTR = map(data, ~{
        ExN = lmer(Value ~  Gradient * poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = .x)
        EplusN = lmer(Value ~  Gradient + poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = .x)
        N = lmer(Value ~  Gradient + (1|GS), REML=FALSE, data = .x)
        E = lmer(Value ~  poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = .x)
        Null = lmer(Value ~  1 + (1|GS), REML=FALSE, data = .x)
        # lr test
        test = anova(ExN, EplusN, N, E, Null)
    })) |>
    unnest(LTR)

  out <- climate_data |>
    # run model by trait
    group_by(Variable) |>
    nest() |>
    mutate(model = map(data, ~ safely(lmer)(Value ~  Gradient * poly(Elevation_m, 2) + (1|GS), data = .x)$result),
           model_output = map(model, tidy),
           r = map(model, r.squaredGLMM),
           r = map(r, as.numeric),
           prediction = map2(.x = data, .y = model, ~ safely(lmer_climate_prediction)(.x, .y)$result))

  return(list(out, lrt, model_choice))

}


# climate figure
make_climate_figure <- function(climate_analysis){

  out <- climate_analysis |>
    # merge data and prediction
    mutate(output = map2(.x = data, .y = prediction, ~ bind_cols(.x, .y))) |>
    select(-data, -prediction, -model, -model_output, -r) |>
    unnest(output) |>
    rename(Value = Value...7, fitted = Value...14, Gradient = Gradient...4, Elevation_m = Elevation_m...9) |>
    select(-Elevation_m...13, -Gradient...12) |>
    mutate(Variable = recode(Variable, "SoilMoisture" = "Soil moisture in %", "SoilTemperature" = "Soil temperature in °C"))

  ggplot(out, aes(x = Elevation_m, y = Value, colour = Gradient)) +
    geom_point(alpha = 0.5) +
    geom_line(aes(y = fitted, colour = Gradient)) +
    geom_ribbon(aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_colour_manual(name = "", values = c("green4", "grey"), labels = c("Nutrient input", "Reference")) +
    scale_fill_manual(name = "", values = c("green4", "grey"), labels = c("Nutrient input", "Reference")) +
    labs(x = "Elevation m a.s.l.", y = "") +
    # add label
    geom_text(aes(x = Inf, y = Inf, label = "NxE"),
              size = 3, colour = "black", hjust = 1, vjust = 1) +
    facet_wrap(~ Variable, scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "top")

}

# dd <- climate_data |>
#   filter(Variable == "SoilTemperature") |>
#   mutate(GS = paste0(Gradient, Site))
# fit <- lmer(Value ~ Gradient * poly(Elevation_m, 2) + (1|GS), dd)
# check_model(fit)
# dd <- trait_mean |>
#   filter(trait_trans == "Thickness_mm_log")
# fit <- lmer(mean ~ Gradient * poly(Elevation_m, 2) + SoilTemperature + (1|GS), dd)
### SLA singular fit

# test linear vs polynomial model
lin_vs_poly_model <- function(trait_mean){

  trait_mean |>
    group_by(trait_trans) |>
    nest() |>
    mutate(lin_mod = map(data, ~ safely(lmer)(mean ~  Gradient * SoilTemperature + (1|GS), data = .x)$result),
           poly_mod = map(data, ~ safely(lmer)(mean ~  Gradient * poly(SoilTemperature, 2) + (1|GS), data = .x)$result),
           mod_aic = map2(.x = lin_mod, .y = poly_mod, .f = AIC)) |>
    unnest(col = mod_aic) |>
    mutate(Model = c("linear", "quadratic"))
}



likelihood_ratio_test_ST <- function(trait_mean){

  trait_mean |>
    #filter(trait_trans != "SLA_cm2_g") |>
    group_by(trait_trans) |>
    nest() |>
    mutate(LTR = map(data, ~{
      # models
      ExN = lmer(mean ~  Gradient * poly(SoilTemperature, 2) + (1|GS), REML=FALSE, data = .x)
      EplusN = lmer(mean ~  Gradient + poly(SoilTemperature, 2) + (1|GS), REML=FALSE, data = .x)
      N = lmer(mean ~  Gradient + (1|GS), REML=FALSE, data = .x)
      E = lmer(mean ~  poly(SoilTemperature, 2) + (1|GS), REML=FALSE, data = .x)
      Null = lmer(mean ~  1 + (1|GS), REML=FALSE, data = .x)
      # ltr
      test = anova(ExN, EplusN, N, E, Null)
    })) |>
    unnest(LTR)

}

likelihood_ratio_test_SM <- function(trait_mean){

  trait_mean |>
    #filter(trait_trans != "SLA_cm2_g") |>
    group_by(trait_trans) |>
    nest() |>
    mutate(LTR = map(data, ~{
      # models
      ExN = lmer(mean ~  Gradient * poly(SoilMoisture, 2) + (1|GS), REML=FALSE, data = .x)
      EplusN = lmer(mean ~  Gradient + poly(SoilMoisture, 2) + (1|GS), REML=FALSE, data = .x)
      N = lmer(mean ~  Gradient + (1|GS), REML=FALSE, data = .x)
      E = lmer(mean ~  poly(SoilMoisture, 2) + (1|GS), REML=FALSE, data = .x)
      Null = lmer(mean ~  1 + (1|GS), REML=FALSE, data = .x)
      # ltr
      test = anova(ExN, EplusN, N, E, Null)
    })) |>
    unnest(LTR)

}

# likelihood_ratio_test_ST(trait_mean) |> filter(`Pr(>Chisq)`<0.05)
# likelihood_ratio_test_SM(trait_mean) |> filter(`Pr(>Chisq)`<0.05)

# make prediction for lmer
lmer_prediction_ST <- function(dat, fit){

  newdat <- dat %>%
    select(Gradient, SoilTemperature) %>%
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

lmer_prediction_SM <- function(dat, fit){

  newdat <- dat %>%
    select(Gradient, SoilMoisture) %>%
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



#
#   out_ST <- fancy_trait_name_dictionary(trait_mean) |>
#     # text for full LRT
#     mutate(text = case_when(trait_trans %in% c("Plant_Height_cm_log", "dC13_permil") ~ "NxE",
#                             trait_trans %in% c("N_percent", "CN_ratio", "dN15_permil") ~ "N",
#                             TRUE ~ "Null")) |>
#     # run model
#     group_by(trait_trans, text) |>
#     nest() |>
#     mutate(
#       # run full model
#       model = map(data, ~ safely(lmer)(mean ~  Gradient * poly(SoilTemperature, 2) + (1|GS), data = .x)$result),
#       # model output
#       model_output = map(model, tidy),
#       r = map(model, r.squaredGLMM),
#       r = map(r, as.numeric),
#       # make model prediction
#       prediction = map2(.x = data, .y = model, ~ safely(lmer_prediction_ST)(.x, .y)$result))
#
#
#   out_SM <- fancy_trait_name_dictionary(trait_mean) |>
#     # text for full LRT
#     mutate(text = case_when(trait_trans %in% c("N_percent", "CN_ratio", "dN15_permil") ~ "N",
#                             TRUE ~ "Null")) |>
#     # run model
#     group_by(trait_trans, text) |>
#     nest() |>
#     mutate(
#       # run full model
#       model = map(data, ~ safely(lmer)(mean ~  Gradient * poly(SoilMoisture, 2) + (1|GS), data = .x)$result),
#       # model output
#       model_output = map(model, tidy),
#       r = map(model, r.squaredGLMM),
#       r = map(r, as.numeric),
#       # make model prediction
#       prediction = map2(.x = data, .y = model, ~ safely(lmer_prediction_SM)(.x, .y)$result))
#
#
# out2_ST <- out_ST |>
#   # merge data and prediction
#   mutate(output = map2(.x = data, .y = prediction, ~ bind_cols(.x, .y))) |>
#   select(-data, -prediction, -model, -model_output, -r) |>
#   unnest(output) |>
#   rename(mean = mean...4, fitted = mean...20, Gradient = Gradient...1, SoilTemperature = SoilTemperature...15) |>
#   select(-SoilTemperature...19, -Gradient...18)
#
# ggplot(out2_ST, aes(x = SoilTemperature, y = mean, colour = Gradient)) +
#   geom_point(alpha = 0.5) +
#   geom_line(aes(y = fitted, colour = Gradient)) +
#   geom_ribbon(aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
#   scale_colour_manual(name = "", values = c("green4", "grey"), labels = c("Nutrient input", "Reference")) +
#   scale_fill_manual(name = "", values = c("green4", "grey"), labels = c("Nutrient input", "Reference")) +
#   labs(x = "Soil temperature in °C", y = "Bootstrapped trait mean") +
#   # add label
#   geom_text(data = out2_ST |>
#               ungroup() |>
#               distinct(trait_trans, trait_fancy, text),
#             aes(x = Inf, y = Inf, label = text),
#             size = 3, colour = "black", hjust = 1, vjust = 1) +
#   facet_wrap(~ trait_fancy, scales = "free_y") +
#   theme_minimal() +
#   theme(legend.position = "top")
#
#
# out2_SM <- out_SM |>
#   # merge data and prediction
#   mutate(output = map2(.x = data, .y = prediction, ~ bind_cols(.x, .y))) |>
#   select(-data, -prediction, -model, -model_output, -r) |>
#   unnest(output) |>
#   rename(mean = mean...4, fitted = mean...20, Gradient = Gradient...1, SoilMoisture = SoilMoisture...14) |>
#   select(-SoilMoisture...19, -Gradient...18)
#
# ggplot(out2_SM, aes(x = SoilMoisture, y = mean, colour = Gradient)) +
#   geom_point(alpha = 0.5) +
#   geom_line(aes(y = fitted, colour = Gradient)) +
#   geom_ribbon(aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
#   scale_colour_manual(name = "", values = c("green4", "grey"), labels = c("Nutrient input", "Reference")) +
#   scale_fill_manual(name = "", values = c("green4", "grey"), labels = c("Nutrient input", "Reference")) +
#   labs(x = "Soil moisture", y = "Bootstrapped trait mean") +
#   # add label
#   geom_text(data = out2_SM |>
#               ungroup() |>
#               distinct(trait_trans, trait_fancy, text),
#             aes(x = Inf, y = Inf, label = text),
#             size = 3, colour = "black", hjust = 1, vjust = 1) +
#   facet_wrap(~ trait_fancy, scales = "free_y") +
#   theme_minimal() +
#   theme(legend.position = "top")
