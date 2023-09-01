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
make_climate_figure <- function(environment_model_output){

  out <- environment_model_output |>
    # merge data and prediction
    mutate(output = map2(.x = data, .y = prediction, ~ bind_cols(.x |> rename(Value = .response, Elevation_m = .continous_predictor),
                                                                 .y |> select(fitted = .response, plo, phi)))) |>
    select(-data, -prediction, -mod, -model_output, -r) |>
    unnest(output) |>
    mutate(text = case_when(Variable == "C" ~ "N",
                            Variable == "N" ~ "N+E",
                            TRUE ~ "NxE")) |>
    mutate(Variable = recode(Variable,
                             "SoilMoisture" = "Soil moisture in %",
                             "SoilTemperature" = "Soil temperature in °C",
                             "C" = "Carbon content in %",
                             "N" = "Nitrogen content in %"),
           Variable = factor(Variable, levels = c("Soil temperature in °C", "Soil moisture in %", "Carbon content in %", "Nitrogen content in %")))

  p <- ggplot(out, aes(x = Elevation_m, y = Value, colour = Gradient)) +
    geom_point(alpha = 0.5) +
    geom_line(aes(y = fitted, colour = Gradient)) +
    geom_ribbon(aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_colour_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient")) +
    scale_fill_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient")) +
    labs(x = "Elevation m a.s.l.", y = "") +
    # add label
    geom_text(aes(x = Inf, y = Inf, label = text),
              size = 3, colour = "black", hjust = 1, vjust = 1) +
    facet_wrap(~ Variable, scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "top")

  p <- tag_facet(p)

  p + theme(strip.text = element_text())

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
                unnest_wider(col = r, names_sep = "_") |>
                rename(Rm = r_1, Rc = r_2) |>
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


### CLIMATE VS TRAIT ANALYSIS

# LRT soil temperature
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

# LRT soil moisture
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


# trait vs. climate figure
make_trait_soil_moisture_figure <- function(soil_moisture_model_output){

out <- soil_moisture_model_output |>
  # merge data and prediction
  mutate(output = map2(.x = data, .y = prediction, ~ bind_cols(.x |> rename(mean = .response, SoilMoisture = .continous_predictor),
                                                               .y |> select(fitted = .response, plo, phi)))) |>
  select(-data, -prediction, -mod, -model_output, -r) |>
  unnest(output) |>
  fancy_trait_name_dictionary() |>
  mutate(figure_names = factor(figure_names, levels = c("Size~-~Height~cm", "Size~-~Dry~mass~g", "Size~-~Area~cm^2", "Size~-~Thickness~mm", "LES~-~SLA~cm^2*g^{-1}", "LES~-~LDMC", "LES~-~C~'%'", "LES~-~N~'%'", "LES~-~CN", "LES~-~P~'%'", "LES~-~NP", "I~-~δ^{13}~C~'‰'", "I~-~δ^{15}~N~'‰'")))

ggplot(out, aes(x = SoilMoisture, y = mean, colour = Gradient)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = fitted, colour = Gradient)) +
  geom_ribbon(aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
  scale_fill_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient")) +
  scale_colour_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient")) +
  labs(x = "Soil moisture in %", y = "Bootstrapped trait mean") +
  # add label
  geom_text(data = out |>
              ungroup() |>
              distinct(trait_trans, figure_names, text),
            aes(x = Inf, y = Inf, label = text),
            size = 3, colour = "black", hjust = 1, vjust = 1) +
  facet_wrap(~ figure_names, scales = "free_y", labeller = label_parsed) +
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
    mutate(figure_names = factor(figure_names, levels = c("Size~-~Height~cm", "Size~-~Dry~mass~g", "Size~-~Area~cm^2", "Size~-~Thickness~mm", "LES~-~SLA~cm^2*g^{-1}", "LES~-~LDMC", "LES~-~C~'%'", "LES~-~N~'%'", "LES~-~CN", "LES~-~P~'%'", "LES~-~NP", "I~-~δ^{13}~C~'‰'", "I~-~δ^{15}~N~'‰'")))

  ggplot(out, aes(x = SoilTemperature, y = mean, colour = Gradient)) +
    geom_point(alpha = 0.5) +
    geom_line(aes(y = fitted, colour = Gradient)) +
    geom_ribbon(aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_fill_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient")) +
    scale_colour_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient")) +
    labs(x = "Soil temperature in °C", y = "Bootstrapped trait mean") +
    # add label
    geom_text(data = out |>
                ungroup() |>
                distinct(trait_trans, figure_names, text),
              aes(x = Inf, y = Inf, label = text),
              size = 3, colour = "black", hjust = 1, vjust = 1) +
    facet_wrap(~ figure_names, scales = "free_y", labeller = label_parsed) +
    theme_minimal() +
    theme(legend.position = "top")
}
