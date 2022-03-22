### RESULT text


trait_mean %>%
  filter(Gradient == "B", trait_trans == "P_percent") %>%
  mutate(hl = if_else(Site == "5", "H", "L")) %>%
  group_by(hl) %>% summarise(mean(mean))

0.328/0.158 # times higher





mutate(Trait = factor(Trait, levels = c("Height cm", "Dry mass g", "Area cm2", "Thickness mm", "LDMC", "SLA cm2/g", "C %", "N %", "CN", "P %", "NP", "δC13 ‰", "δN15 ‰"))) |>



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


tar_load(top_site)
fancy_trait_name_dictionary(top_site) |>
  ungroup() |>
  mutate(term = if_else(npar == 4, "Site", "Null"),
         AIC = round(AIC, digits = 1),
         `Pr(>Chisq)` = round(`Pr(>Chisq)`, digits = 3)) |>
  select(trait_fancy, term, AIC, df = Df, "p value" = `Pr(>Chisq)`) %>%
  write_csv(., file = "output/top_site.csv")

