### RESULT text


trait_mean %>%
  filter(Gradient == "B", trait_trans == "P_percent") %>%
  mutate(hl = if_else(Site == "5", "H", "L")) %>%
  group_by(hl) %>% summarise(mean(mean))

0.328/0.158 # times higher

fancy_trait_name_dictionary(trait_mean) %>%
  group_by(Gradient, trait_fancy) %>%
  summarise(se = round(sd(mean)/sqrt(n()), 2),
            mean = round(mean(mean), 2),
            se_var = round(sd(var)/sqrt(n()), 2),
            var = round(mean(var), 2)) %>%
  select(Trait = trait_fancy, Gradient, Mean = mean, se, Variance = var, se_var) %>%
  write_csv(file = "output/Mean_var.csv")



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
  select(-trait_trans, Term = term, Estimate = estimate, "Std error" = std.error, "t-value" = "statistic", "Marginal R2" = Rm, "Conditional R2" = Rc) %>%
  arrange(Trait) %>%
  write_csv(., file = "output/Regression_output.csv")



# trait ordination output
tar_load(trait_pca_B)
tar_load(trait_pca_C)
tar_load(trait_pca)
bind_rows(Birdcliff = trait_pca_B[[2]],
          Reference = trait_pca_C[[2]] %>%
            mutate(PC1 = PC1 * -1,
                   PC2 = PC2 * -1,
                   PC3 = PC3 * -1,
                   PC4 = PC4 * -1),
          Both = trait_pca[[2]],
          .id = "Gradient") %>%
  select(Gradient, Trait = trait_fancy, PC1:PC4) %>%
  mutate(PC1 = round(PC1, digits = 2),
         PC2 = round(PC2, digits = 2),
         PC3 = round(PC3, digits = 2),
         PC4 = round(PC4, digits = 2)) %>%
  mutate(Gradient = recode(Gradient, Birdcliff = "Bird cliff"),
         Gradient = factor(Gradient, levels = c("Bird cliff", "Reference", "Both"))) %>%
  arrange(Gradient, -PC1) %>%
  write_csv(., file = "output/Loadings_trait_PCA.csv")



bind_rows(Birdcliff = tidy(trait_pca_B[[4]]),
          Reference = tidy(trait_pca_C[[4]]),
          Both = tidy(trait_pca[[4]]),
          .id = "Gradient") %>%
  write_csv(., file = "output/adonis_pca.csv")



### Ind Species
best_ind_model <- ind_traits %>%
  filter(Functional_group == "vascular") %>%
  distinct(Taxon, trait_trans) %>%
  mutate(best = c("NA", "Null", "G+E", "Null", "G", "G+E", "G+E", "GxE", "G", "Gx"))

tar_load(ind_traits_output)
fancy_trait_name_dictionary(ind_traits_output[[1]]) %>%
  filter(effect == "fixed") %>%
  left_join(best_ind_model, by = c("Taxon", "trait_trans")) %>%
  select(Trait = trait_fancy, Model = best, term, estimate:statistic) %>%
  mutate(term = recode(term, "(Intercept)" = "Intercept", "Elevation_m" = "E", "GradientC" = "G", "GradientC:Elevation_m" = "GxE")) %>%
  left_join(model_output[[2]] %>%
              select(trait_trans, Rm, Rc), by = "trait_trans") %>%
  mutate(estimate = round(estimate, digits = 2),
         std.error = round(std.error, digits = 2),
         statistic = round(statistic, digits = 2),
         Rm = round(Rm, digits = 2),
         Rc = round(Rc, digits = 2)) %>%
  ungroup() %>%
  select(-trait_trans, Term = term, Estimate = estimate, "Std error" = std.error, "t-value" = "statistic", "Marginal R2" = Rm, "Conditional R2" = Rc) %>%
  arrange(Taxon, Trait) %>%
  write_csv(., file = "output/Ind_sp_regression_output.csv")



