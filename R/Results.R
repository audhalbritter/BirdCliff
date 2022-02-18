### RESULT text


trait_mean %>%
  filter(Gradient == "B", trait_trans == "P_percent") %>%
  mutate(hl = if_else(Site == "5", "H", "L")) %>%
  group_by(hl) %>% summarise(mean(mean))

0.328/0.158 # times higher



# Regression output
tar_load(model_output)
fancy_trait_name_dictionary(model_output[[1]]) %>%
  filter(effect == "fixed") %>%
  select(Trait = trait_fancy, Model, term, estimate:statistic) %>%
  mutate(term = recode(term, "(Intercept)" = "Intercept", "Elevation_m" = "E", "GradientC" = "G", "GradientC:Elevation_m" = "GxE")) %>%
  left_join(model_output[[2]] %>%
              select(trait_trans, Rm, Rc), by = "trait_trans") %>%
  mutate(estimate = round(estimate, digits = 2),
         std.error = round(std.error, digits = 2),
         statistic = round(statistic, digits = 2),
         Rm = round(Rm, digits = 2),
         Rc = round(Rc, digits = 2)) %>%
  ungroup() %>%
  select(-trait_trans, Term = term, Estimate = estimate, "Std error" = std.error, "t value" = "statistic") %>%
  write_csv(., file = "output/Regression_output.csv")



# trait ordination output
tar_load(trait_pca_B)
tar_load(trait_pca_C)
bind_rows(Birdcliff = trait_pca_B[[2]],
          Reference = trait_pca_C[[2]],
          .id = "Gradient") %>%
  select(Gradient, Trait = trait_fancy, PC1:PC4) %>%
  mutate(PC1 = round(PC1, digits = 2),
         PC2 = round(PC2, digits = 2),
         PC3 = round(PC3, digits = 2),
         PC4 = round(PC4, digits = 2)) %>%
  mutate(Gradient = recode(Gradient, Birdcliff = "Bird cliff")) %>%
  write_csv(., file = "output/Loadings_trait_PCA.csv")


