#### INDIVIDUAL LEVEL TRAITS

# combine vascular and bryo traits
combine_traits <- function(traits_raw, bryo_traits_raw){

  # combine species vascular and bryophytes
  ind_traits <- bind_rows(traits_raw %>%
                            filter(Taxon %in% c("luzula confusa", "salix polaris"),
                                   trait_trans %in% c("Plant_Height_cm_log", "Leaf_Area_cm2_log", "LDMC", "N_percent", "dN15_permil")),
                          bryo_traits_raw %>%
                            filter(Taxon %in% c("aulacomnium turgidum", "hylocomium splendens", "sanionia sp"),
                                   trait_trans %in% c("Shoot_Length_cm_log", "Shoot_ratio", "SSL_cm_g", "WHC_g_g", "P_percent"))) %>%
    mutate(trait_trans = factor(trait_trans, levels = c("Plant_Height_cm_log", "Leaf_Area_cm2_log", "LDMC", "N_percent", "dN15_permil", "Shoot_Length_cm_log", "Shoot_ratio", "SSL_cm_g", "WHC_g_g", "P_percent")),
           Taxon = factor(Taxon, levels = c("luzula confusa", "salix polaris", "aulacomnium turgidum", "hylocomium splendens", "sanionia sp"))) %>%
    mutate(GS = paste0(Gradient, Site))

  return(ind_traits)

}


### BRYOPHYTES

# run linear and quadratic model
run_bryophyte_model <- function(dat, group, response, continous_predictor){

  dat |>
    rename(.response = {{response}},
           .continous_predictor = {{continous_predictor}}) |>
    group_by(across(all_of({{group}})))|>
    nest() |>
    mutate(mod_linear = map(data, ~ safely(lm)(.response ~  Gradient * .continous_predictor, data = .x)$result),
           mod_quadratic = map(data, ~ safely(lm)(.response ~  Gradient * poly(.continous_predictor, 2), data = .x)$result),
           aic_linear = map(.x = mod_linear, .f = ~ safely(AIC)(.x)$result),
           aic_linear = as.numeric(aic_linear),
           aic_quadratic = map(.x = mod_quadratic, .f = ~ safely(AIC)(.x)$result),
           aic_quadratic = as.numeric(aic_quadratic))
}



# Bryophyte table
make_bryo_table <- function(bryophyte_model_output){

  fancy_trait_name_dictionary(bryophyte_model_output) |>
    unnest(model_output) %>%
    ungroup() |>
    select(-data, -names, -mod, -aic, -prediction) |>
    mutate(term = recode(term,
                         "(Intercept)" = "Intercept",
                         ".continous_predictor" = "E",
                         "GradientB" = "G",
                         "GradientB:.continous_predictor" = "GxE",
                         "poly(.continous_predictor, 2)1" = "E",
                         "poly(.continous_predictor, 2)2" = "E^2^",
                         "GradientB:poly(.continous_predictor, 2)1" = "GxE",
                         "GradientB:poly(.continous_predictor, 2)2" = "GxE^2^")

           ) %>%
    mutate(estimate = round(estimate, digits = 2),
           std.error = round(std.error, digits = 2),
           statistic = round(statistic, digits = 2),
           p.value = round(p.value, digits = 3)) %>%
    select(Class = class, Taxon, Trait = trait_fancy, Term = term, Estimate = estimate, "Std. error" = std.error, "t value" = statistic, "p value" = p.value) %>%
    arrange(Taxon, Trait) %>%
    write_csv(., file = "output/bryo_trait_output.csv")

}

# make bryophyte figure
make_bryo_figure <- function(bryophyte_model_output){

  out <- fancy_trait_name_dictionary(bryophyte_model_output) %>%
    select(-mod, -aic, -model_output) |>
    unnest(prediction) |>
    mutate(text = case_when(trait_trans == "Shoot_ratio" ~ "",
                            Taxon == "aulacomnium turgidum" & trait_trans == "Shoot_Length_cm_log" ~ "",
                            Taxon == "sanionia sp" & trait_trans == "P_percent" ~ "E",
                            Taxon == "hylocomium splendens" & trait_trans == "Shoot_Length_cm_log" ~ "E",
                            Taxon == "hylocomium splendens" & trait_trans == "WHC_g_g" ~ "N",
                            TRUE ~ "NxE")) |>
    mutate(Taxon = recode(Taxon, "sanionia sp" = "Sanionia sp", "hylocomium splendens" = "Hylocomium splendens", "aulacomnium turgidum" = "Aulacomnium turgidum"))


 ggplot(out, aes(x = .continous_predictor, y = .response, colour = Gradient)) +
    geom_point(alpha = 0.5) +
    geom_line(aes(y = .fitted, colour = Gradient)) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_colour_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
    scale_fill_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
    labs(x = "Elevation in m a.s.l.", y = "Trait value") +
    geom_text(data = out |>
                distinct(Taxon, trait_fancy, text),
              aes(x = Inf, y = Inf, label = text),
              size = 3, colour = "black", hjust = 1, vjust = 1) +
    facet_grid(trait_fancy ~ Taxon, scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "top",
          plot.title = element_text(size = 10),
          strip.text.x = element_text(face = "italic"),
          aspect.ratio = 0.7)

}
