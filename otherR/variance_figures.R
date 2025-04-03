#### Trait variance figure


t <- trait_mean |>
  filter(trait_trans %in% c("Dry_Mass_g_log", "SLA_cm2_g", "LDMC", "N_percent", "P_percent", "dN15_permil")) |>
  fancy_trait_name_dictionary()

# pca
pca_plot <- trait_pca[[1]] |>
  ggplot(aes(x = factor(Site), y = PC1, fill = Gradient, colour = Gradient)) +
  geom_violin(alpha = 0.4) +
  geom_point(position=position_jitterdodge()) +
  scale_fill_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
  scale_colour_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
  scale_x_discrete(labels = c("1" = "10.8", "2" = "41.8", "3" = "84.0", "4" = "123.0", "5" = "174.0", "6" = "224.0", "7" = "238.0")) +
  labs(x = "Elevation m a.s.l.",
       y = "PCA1 score") +
  geom_text(aes(x = Inf, y = Inf, label = "NxE"),
            size = 3, colour = "black", hjust = 1, vjust = 1) +
  theme_minimal() +
  theme(legend.position = "top")

# ridges
library()
ridge <- t |>
  filter(trait_trans %in% c("SLA_cm2_g", "dN15_permil")) |>
  ggplot(aes(x = mean, y = factor(Site), fill = Gradient)) +
  geom_density_ridges(alpha = 0.4, colour = NA) +
  scale_fill_manual(name = "", values = c("darkgrey", "green4"), labels = c("Reference", "Nutrient input")) +
  scale_y_discrete(labels = c("1" = "10.8", "2" = "41.8", "3" = "84.0", "4" = "123.0", "5" = "174.0", "6" = "224.0", "7" = "238.0")) +
  labs(x = "Bootstrapped trait mean",
       y = "Elevation in m a.s.l.") +
  facet_wrap( ~ trait_fancy, scales = "free") +
  theme_minimal() +
  theme(legend.position = "none")

variance_plot2 <- ridge / pca_plot + plot_layout(height = c(2, 1.5)) +
  plot_annotation(tag_levels = 'a', tag_suffix = ")")
ggsave("output/variance_plot2.jpg", variance_plot2, dpi = 300, height = 6, width = 6, bg = "white")





### FUNCTIONAL DIVERSITY
# MULTIVARIATE BOOTSTRAPPING
tar_target(
  name = multi_trait,
  command = make_multi_trait_impute(comm_raw, traits_raw |>
                                      # without ratios
                                      filter(!Trait %in% c("CN_ratio", "NP_ratio")))
)

# no chem traits
tar_target(
  name = multi_trait_leaf,
  command = make_multi_trait_impute(comm_raw,
                                    # remove all chemical traits
                                    traits_raw |>
                                      filter(!Trait %in% c("C_percent", "N_percent", "CN_ratio", "dN15_permil", "dC13_permil", "P_percent", "NP_ratio")))
)

# make bootstrap
tar_target(
  name = func_diversity,
  command = make_multi_bootstrap(multi_trait)
)

# make bootstrap
tar_target(
  name = func_diversity_leaf,
  command = make_multi_bootstrap(multi_trait_leaf)
)


# data wrangling
multi_out <- func_diversity %>%
  mutate(result = map(result, as.data.frame)) %>%
  unnest(result) |>
  pivot_longer(cols = c(FEve:RaoQ), names_to = "variable", values_to = "value") |>
  mutate(Gradient = factor(Gradient, levels = c("C", "B")),
         GS = paste0(Gradient, Site))


multi_model <- run_trait_model(dat = multi_out,
                group = "variable",
                response = value,
                continous_predictor = Elevation_m) |>
  pivot_longer(cols = -c(variable, data),
               names_sep = "_",
               names_to = c(".value", "names"))

multi_model2 <- multi_model |>
  #filter(names == "linear")
  # remove models that have singular fit
  filter(singular == FALSE) |>
  filter(aic == min(aic))

likelihood_ratio_test(multi_out |> rename(trait_trans = variable))

dat <- model_output_prediction(multi_model2) |>
  # merge data and prediction
  mutate(output = map2(.x = data, .y = prediction, ~ bind_cols(.x, .y))) |>
  select(-data, -singular, -aic, -prediction, -mod, -model_output, -r) |>
  unnest(output) %>%
  rename(Gradient = Gradient...2, mean = .response...8, Elevation_m = .continous_predictor...4, fitted = .response...12) |>
  select(-Gradient...10, -.continous_predictor...11) %>%
  mutate(variable = recode(variable,
                           "Evenness" = "FEve",
                           "Dispersion" = "FDis"))

ggplot(dat, aes(x = Elevation_m, y = mean, colour = Gradient)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = fitted, colour = Gradient)) +
  geom_ribbon(aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
  scale_colour_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
  scale_fill_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
  labs(x = "Elevation m a.s.l.", y = "Bootstrapped trait mean") +
  # add label
  # geom_text(data = dat %>%
  #             ungroup() |>
  #             distinct(trait_trans, trait_fancy, text),
  #           aes(x = Inf, y = Inf, label = text),
  #           size = 3, colour = "black", hjust = 1, vjust = 1) +
  facet_wrap(~ variable, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "top")

multi_out |>
  filter(variable == "FDis") |>
  ggplot(aes(x = Elevation_m, y = value, colour = Gradient)) +
  geom_point(alpha = 0.5) +
  labs(x = "Elevation in m a.s.l.", y = "Dispersion") +
  scale_colour_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
  #facet_wrap(~ variable, scales = "free_y") +
  theme_minimal()


dat <- multi_out |>
  filter(variable == "FDis")

fit <- lm(value ~ Gradient * Elevation_m, dat)
summary(fit)
dat$fitted = predict(fit, dat)

dat |>
  ggplot(aes(x = Elevation_m, y = value, colour = Gradient)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = fitted)) +
  labs(x = "Elevation in m a.s.l.", y = "Dispersion") +
  scale_colour_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
  #facet_wrap(~ variable, scales = "free_y") +
  theme_minimal()
