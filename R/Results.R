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





tar_load(top_site)
fancy_trait_name_dictionary(top_site) |>
  ungroup() |>
  mutate(term = if_else(npar == 4, "Site", "Null"),
         AIC = round(AIC, digits = 1),
         `Pr(>Chisq)` = round(`Pr(>Chisq)`, digits = 3)) |>
  select(trait_fancy, term, AIC, df = Df, "p value" = `Pr(>Chisq)`) %>%
  write_csv(., file = "output/top_site.csv")

