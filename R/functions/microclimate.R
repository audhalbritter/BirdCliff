### CLIMATE DATA

climate_lrt <- function(climate_data){

  climate_data |>
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

}


# climate figure
make_climate_figure <- function(climate_model_output){

  out <- climate_model_output |>
    # merge data and prediction
    mutate(output = map2(.x = data, .y = prediction, ~ bind_cols(.x |> rename(Value = .response, Elevation_m = .continous_predictor),
                                                                 .y |> select(fitted = .response, plo, phi)))) |>
    select(-data, -prediction, -mod, -model_output, -r) |>
    unnest(output) |>
    mutate(Variable = recode(Variable, "SoilMoisture" = "Soil moisture in %", "SoilTemperature" = "Soil temperature in °C")) |>
    mutate(text = "NxE")

  ggplot(out, aes(x = Elevation_m, y = Value, colour = Gradient)) +
    geom_point(alpha = 0.5) +
    geom_line(aes(y = fitted, colour = Gradient)) +
    geom_ribbon(aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_colour_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
    scale_fill_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
    labs(x = "Elevation m a.s.l.", y = "") +
    # add label
    geom_text(aes(x = Inf, y = Inf, label = "NxE"),
              size = 3, colour = "black", hjust = 1, vjust = 1) +
    facet_wrap(~ Variable, scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "top")

}


# climate output: model output and r2
make_climate_table <- function(climate_model_output){

  climate_model_output |>
    select(-data, -mod, -singular, -aic, -prediction, -r) |>
    unnest(model_output) %>%
    filter(effect == "fixed") %>%
    ungroup() |>
    select(Variable, term, estimate:statistic) %>%
    mutate(term = recode(term,
                         "(Intercept)" = "Intercept",
                         "GradientB" = "N",
                         "poly(.continous_predictor, 2)1" = "E",
                         "poly(.continous_predictor, 2)2" = "E2",
                         "GradientB:poly(.continous_predictor, 2)1" = "NxE",
                         "GradientB:poly(.continous_predictor, 2)2" = "NxE2")) %>%
    # join R squared
    left_join(climate_model_output |>
                select(-data, -mod, -singular, -aic, -prediction, -model_output) |>
                unnest_wider(col = r) |>
                rename(Rm = ...1, Rc = ...2) |>
                ungroup() %>%
                select(Variable, Rm, Rc),
              by = "Variable") %>%
    mutate(estimate = round(estimate, digits = 2),
           std.error = round(std.error, digits = 2),
           statistic = round(statistic, digits = 2),
           Rm = round(Rm, digits = 2),
           Rc = round(Rc, digits = 2)) %>%
    select(Variable, Term = term, Estimate = estimate, "Std error" = std.error, "t-value" = "statistic", "Marginal R2" = Rm, "Conditional R2" = Rc) %>%
    arrange(Variable) %>%
    write_csv(., file = "output/Climate_output.csv")

}


### CLIMATE AND TRAITS

# APPROACH 1 comparing model with and without microclimate
# # run linear and quadratic model
# run_trait_climate_model <- function(dat, group, response, continous_predictor){
#
#   dat |>
#     rename(.response = {{response}},
#            .continous_predictor = {{continous_predictor}}) |>
#     group_by(across(all_of({{group}})))|>
#     nest() |>
#     mutate(
#       # linear model
#       mod_linear = map(data, ~ safely(lmer)(.response ~  Gradient * .continous_predictor + SoilTemperature + (1|GS), data = .x)$result),
#       singular_linear = map_lgl(mod_linear, isSingular),
#       aic_linear = map(.x = mod_linear, .f = ~ safely(AIC)(.x)$result),
#       aic_linear = as.numeric(aic_linear),
#
#       # quadratic
#       mod_quadratic = map(data, ~ safely(lmer)(.response ~  Gradient * poly(.continous_predictor, 2) + (1|GS), data = .x)$result),
#       singular_quadratic = map_lgl(mod_quadratic, isSingular),
#       aic_quadratic = map(.x = mod_quadratic, .f = ~ safely(AIC)(.x)$result),
#       aic_quadratic = as.numeric(aic_quadratic),
#
#       # quadratic plus SoilTemp
#       mod_quadraticPlus = map(data, ~ safely(lmer)(.response ~  Gradient * poly(.continous_predictor, 2) + SoilTemperature + (1|GS), data = .x)$result),
#       aic_quadraticPlus = map(.x = mod_quadraticPlus, .f = ~ safely(AIC)(.x)$result),
#       aic_quadraticPlus = as.numeric(aic_quadraticPlus))
# }
#
# # run linear and quadratic model
# tar_target(
#   name = trait_community_model,
#   command = run_trait_climate_model(dat = trait_mean,
#                             group = "trait_trans",
#                             response = mean,
#                             continous_predictor = Elevation_m) |>
#     pivot_longer(cols = -c(trait_trans, data),
#                  names_sep = "_",
#                  names_to = c(".value", "names")) |>
#     filter(singular == FALSE | is.na(singular)) |>
#     filter(aic == min(aic))
# ),


### APPROACH 2

likelihood_ratio_test_ST <- function(trait_mean){

  trait_mean |>
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


# trait climate figure
make_trait_soil_moisture_figure <- function(soil_moisture_model_output){

out <- soil_moisture_model_output |>
  # merge data and prediction
  mutate(output = map2(.x = data, .y = prediction, ~ bind_cols(.x |> rename(mean = .response, SoilMoisture = .continous_predictor),
                                                               .y |> select(fitted = .response, plo, phi)))) |>
  select(-data, -prediction, -mod, -model_output, -r) |>
  unnest(output) |>
  fancy_trait_name_dictionary() |>
  mutate(class = recode(class, "Leaf economics" = "LES", "Isotopes" = "I"),
         trait_fancy = paste(class, trait_fancy, sep = " - "),
         trait_fancy = factor(trait_fancy, levels = c("Size - Height cm", "Size - Dry mass g", "Size - Area cm2", "Size - Thickness mm", "LES - SLA cm2/g", "LES - LDMC", "LES - C %", "LES - N %", "LES - CN", "LES - P %", "LES - NP", "I - δC13 ‰", "I - δN15 ‰")))

ggplot(out, aes(x = SoilMoisture, y = mean, colour = Gradient)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = fitted, colour = Gradient)) +
  geom_ribbon(aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
  scale_fill_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
  scale_colour_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
  labs(x = "Soil moisture in %", y = "Bootstrapped trait mean") +
  # add label
  geom_text(data = out |>
              ungroup() |>
              distinct(trait_trans, trait_fancy, text),
            aes(x = Inf, y = Inf, label = text),
            size = 3, colour = "black", hjust = 1, vjust = 1) +
  facet_wrap(~ trait_fancy, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "top")
}


make_trait_soil_temp_figure <- function(soil_temp_model_output){

  out <- soil_temp_model_output |>
    # merge data and prediction
    mutate(output = map2(.x = data, .y = prediction, ~ bind_cols(.x |> rename(mean = .response, SoilTemperature = .continous_predictor),
                                                                 .y |> select(fitted = .response, plo, phi)))) |>
    select(-data, -prediction, -mod, -model_output, -r) |>
    unnest(output) |>
    fancy_trait_name_dictionary() |>
    mutate(class = recode(class, "Leaf economics" = "LES", "Isotopes" = "I"),
           trait_fancy = paste(class, trait_fancy, sep = " - "),
           trait_fancy = factor(trait_fancy, levels = c("Size - Height cm", "Size - Dry mass g", "Size - Area cm2", "Size - Thickness mm", "LES - SLA cm2/g", "LES - LDMC", "LES - C %", "LES - N %", "LES - CN", "LES - P %", "LES - NP", "I - δC13 ‰", "I - δN15 ‰")))

  ggplot(out, aes(x = SoilTemperature, y = mean, colour = Gradient)) +
    geom_point(alpha = 0.5) +
    geom_line(aes(y = fitted, colour = Gradient)) +
    geom_ribbon(aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_fill_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
    scale_colour_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
    labs(x = "Soil temperature in °C", y = "Bootstrapped trait mean") +
    # add label
    geom_text(data = out |>
                ungroup() |>
                distinct(trait_trans, trait_fancy, text),
              aes(x = Inf, y = Inf, label = text),
              size = 3, colour = "black", hjust = 1, vjust = 1) +
    facet_wrap(~ trait_fancy, scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "top")
}


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
