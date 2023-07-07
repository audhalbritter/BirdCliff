
run_simple_model <- function(dat, group, response, continous_predictor){

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

# make prediction for lmer
# lmer_prediction <- function(dat, fit){
#
#   newdat <- dat %>%
#     select(.continous_predictor)
#
#   newdat$.response <- predict(fit, newdat, re.form = NA)
#
#   mm <- model.matrix(terms(fit), newdat)
#
#   prediction <- newdat %>%
#     mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
#            tvar1 = pvar1 + VarCorr(fit)$GS[1],  ## must be adapted for more complex models
#            cmult = 1.96) %>%
#     mutate(plo = .response - cmult*sqrt(pvar1),
#            phi = .response + cmult*sqrt(pvar1),
#            tlo = .response - cmult*sqrt(tvar1),
#            thi = .response + cmult*sqrt(tvar1))
#
#   return(prediction)
# }

dat <- trait_mean |>
  select(Gradient:mean, GS) |>
  pivot_wider(names_from = trait_trans, values_from = mean) |>
  pivot_longer(cols = c(Plant_Height_cm_log:dC13_permil), names_to = "trait_trans", values_to = "mean")

model <- run_simple_model(dat = dat,
                group = "trait_trans",
                response = mean,
                continous_predictor = dN15_permil) |>
  pivot_longer(cols = -c(trait_trans, data),
               names_sep = "_",
               names_to = c(".value", "names")) |>
  # remove models that have singular fit
  filter(singular == FALSE) |>
  filter(aic == min(aic))


model_output <- model_output_prediction(model) |>
  # add LRT text
  mutate(text = case_when(trait_trans %in% c("LDMC", "P_percent") ~ "Null",
                          trait_trans %in% c("CN_ratio") ~ "N",
                          trait_trans %in% c("Thickness_mm_log", "N_percent", "NP_ratio") ~ "dN15",
                          trait_trans %in% c("Plant_Height_cm_log", "Leaf_Area_cm2_log", "dC13_permil") ~ "N+dN15",
                          TRUE ~ "NxdN15"))


dat |>
  #filter(trait_trans == "Plant_Height_cm_log") |>
  group_by(trait_trans) |>
  nest() |>
  mutate(LTR = map(data, ~{
      ExN = lmer(mean ~  Gradient * poly(dN15_permil, 2) + (1|GS), REML=FALSE, data = .x)
      EplusN = lmer(mean ~  Gradient + poly(dN15_permil, 2) + (1|GS), REML=FALSE, data = .x)
      N = lmer(mean ~  Gradient + (1|GS), REML=FALSE, data = .x)
      E = lmer(mean ~  poly(dN15_permil, 2) + (1|GS), REML=FALSE, data = .x)
      Null = lmer(mean ~  1 + (1|GS), REML=FALSE, data = .x)
      # lr test
      test = anova(ExN, EplusN, N, E, Null)
  })) |>
  unnest(LTR)



dat2 <- model_output |>
  # merge data and prediction
  mutate(output = map2(.x = data, .y = prediction, ~ bind_cols(.x, .y))) |>
  select(-data, -singular, -aic, -prediction, -mod, -model_output, -r) |>
  unnest(output) %>%
  #rename(mean = .response...6, dN15 = .continous_predictor...5, fitted = .response...8) |>
  rename(Gradient = Gradient...1, mean = .response...6, dN15 = .continous_predictor...5, fitted = .response...9) |>
  #select(-.continous_predictor...7) %>%
  select(-Gradient...7, -.continous_predictor...8) %>%
  fancy_trait_name_dictionary(.) |>
  mutate(class = recode(class, "Leaf economics" = "LES", "Isotopes" = "I"),
         trait_fancy = paste(class, trait_fancy, sep = " - "),
         trait_fancy = factor(trait_fancy, levels = c("Size - Height cm", "Size - Dry mass g", "Size - Area cm2", "Size - Thickness mm", "LES - SLA cm2/g", "LES - LDMC", "LES - C %", "LES - N %", "LES - CN", "LES - P %", "LES - NP", "I - δC13 ‰")))


dn15_figure <- ggplot(dat2, aes(x = dN15, y = mean, colour = Gradient)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = fitted, colour = Gradient)) +
  geom_ribbon(aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
  scale_colour_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
  scale_fill_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
  labs(x = "dN15 in permil", y = "Bootstrapped trait mean") +
  # add label
  geom_text(data = dat2 %>%
              ungroup() |>
              distinct(trait_trans, trait_fancy, text),
            aes(x = Inf, y = Inf, label = text),
            size = 3, colour = "black", hjust = 1, vjust = 1) +
  facet_wrap(~ trait_fancy, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "top")

ggsave("output/dn15_plot.jpg", dn15_figure, dpi = 300, height = 6, width = 8, bg = "white")
