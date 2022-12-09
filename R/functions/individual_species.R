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


#### VASCULAR PLANTS

likelihood_ratio_test_vasc <- function(ind_traits){

  ind_traits |>
    filter(Functional_group == "vascular") |>
    group_by(Taxon, trait_trans) |>
    nest() |>
    mutate(LTR = map(data, ~{

      if (Taxon == "salix polaris" & trait_trans == "dN15_permil") {

        # linear model
        ExN = lmer(Value ~  Gradient * Elevation_m + (1|GS), REML=FALSE, data = .x)
        EplusN = lmer(Value ~  Gradient + Elevation_m + (1|GS), REML=FALSE, data = .x)
        N = lmer(Value ~  Gradient + (1|GS), REML=FALSE, data = .x)
        E = lmer(Value ~  Elevation_m + (1|GS), REML=FALSE, data = .x)
        Null = lmer(Value ~  1 + (1|GS), REML=FALSE, data = .x)
        # lr test
        test = anova(ExN, EplusN, N, E, Null)

      } else {

        # quadratic model
        ExN = lmer(Value ~  Gradient * poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = .x)
        EplusN = lmer(Value ~  Gradient + poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = .x)
        N = lmer(Value ~  Gradient + (1|GS), REML=FALSE, data = .x)
        E = lmer(Value ~  poly(Elevation_m, 2) + (1|GS), REML=FALSE, data = .x)
        Null = lmer(Value ~  1 + (1|GS), REML=FALSE, data = .x)
        # lr test
        test = anova(ExN, EplusN, N, E, Null)

      }
    })) |>
    unnest(LTR)

}


# trait figure
make_vascular_figure <- function(vascular_model_output){

  out <- vascular_model_output |>
    # merge data and prediction
    mutate(output = map2(.x = data, .y = prediction, ~ bind_cols(.x, .y))) |>
    select(-data, -singular, -aic, -prediction, -mod, -model_output, -r) |>
    unnest(output) |>
    rename(Value = .response...11, fitted = .response...19, Gradient = Gradient...4, Elevation_m = .continous_predictor...12) |>
    select(-.continous_predictor...18, -Gradient...17) %>%
    fancy_trait_name_dictionary(.) |>
    mutate(Taxon = recode(Taxon,
                          "luzula confusa" = "Luzula confusa",
                          "salix polaris" = "Salix polaris"))

  out |>
    ggplot(aes(x = Elevation_m, y = Value, colour = Gradient)) +
    geom_point(alpha = 0.5) +
    geom_line(aes(y = fitted, colour = Gradient)) +
    geom_ribbon(aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_colour_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
    scale_fill_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
    labs(x = "Elevation m a.s.l.", y = "Trait mean") +
    # add label
    geom_text(data = out |>
                ungroup() |>
                distinct(Taxon, trait_trans, trait_fancy, text),
              aes(x = Inf, y = Inf, label = text),
              size = 3, colour = "black", hjust = 1, vjust = 1) +
    facet_grid(trait_fancy ~ Taxon, scales = "free") +
    theme_minimal() +
    theme(legend.position = "top",
          strip.text.x = element_text(face = "italic"))

}


### BRYOPHYTES

lin_poly_bryophytes <- function(ind_traits){

  ind_traits |>
    filter(Functional_group == "bryophyte") |>
    group_by(Taxon, trait_trans) |>
    nest() |>
    mutate(lin_mod = map(data, ~ safely(lm)(value_trans ~  Gradient * Elevation_m, data = .x)$result),
           poly_mod = map(data, ~ safely(lm)(value_trans ~  Gradient * poly(Elevation_m, 2), data = .x)$result),
           mod_aic = map2(.x = lin_mod, .y = poly_mod, .f = AIC)) |>
    unnest(col = mod_aic) |>
    mutate(Model = c("linear", "quadratic")) |>
    filter(AIC == min(AIC))

}


# bryophyte models
make_bryo_trait_model <- function(ind_traits){

  bryo_trait_output <- fancy_trait_name_dictionary(ind_traits) %>%
    filter(Functional_group == "bryophyte") %>%
    group_by(Taxon, trait_fancy) %>%
    nest(data = -c(Taxon, trait_fancy)) %>%
    mutate(estimate = map(data, ~{
      mod <- lm(value_trans ~ Gradient*Elevation_m, data =  .x)
      estimates = tidy(mod)
    })) %>%
    unnest(estimate) %>%
    mutate(term = recode(term, "(Intercept)" = "Intercept", "Elevation_m" = "E", "GradientC" = "G", "GradientC:Elevation_m" = "GxE")) %>%
    mutate(estimate = round(estimate, digits = 2),
           std.error = round(std.error, digits = 2),
           statistic = round(statistic, digits = 2),
           p.value = round(p.value, digits = 3)) %>%
    select(Taxon, Trait = trait_fancy, Term = term, Estimate = estimate, "Std. error" = std.error, "t value" = statistic, "p value" = p.value) %>%
    arrange(Taxon, Trait)

  write_csv(bryo_trait_output, file = "output/bryo_trait_output.csv")

  return(bryo_trait_output)

}



make_bryo_figure <- function(ind_traits, bryo_trait_output){

  ### BRYO FIGURE

  term = tibble(term = c("E", " ", "NxE", "NxE", "NxE", " ", " ", "N", " ", "NxE", " ", " ", " ", " ", " "))

  result <- bryo_trait_output |>
    distinct(Taxon, Trait) |>
   bind_cols(term = term) |>
    mutate(pvalue = if_else(term == " ", "non-sign", "sign"))

  b_dat <- fancy_trait_name_dictionary(ind_traits) %>%
    filter(Functional_group == "bryophyte") |>
    left_join(result, by = c("Taxon", "trait_fancy" = "Trait")) |>
    mutate(pvalue = factor(pvalue, levels = c("sign", "non-sign")),
           Taxon = recode(Taxon, "sanionia sp" = "Sanionia sp", "hylocomium splendens" = "Hylocomium splendens", "aulacomnium turgidum" = "Aulacomnium turgidum"))

  bryo_plot <- b_dat |>
    ggplot(aes(x = Elevation_m, y = value_trans, colour = Gradient, linetype = pvalue, shape = class)) +
    geom_point(alpha = 0.5) +
    geom_smooth(mapping = aes(fill = Gradient), b_dat |> filter(pvalue == "sign"), method = "lm", se = TRUE) +
    scale_colour_manual(name = "", values = c("green4", "grey"), labels = c("Nutrient input", "Reference")) +
    scale_fill_manual(name = "", values = c("green4", "grey"), labels = c("Nutrient input", "Reference")) +
    scale_shape_manual(values = c(16, 17)) +
    labs(x = "Elevation in m a.s.l.", y = "Trait value") +
    geom_text(ata = b_dat |> distinct(Taxon, trait_fancy, term) , aes(x = Inf, y = Inf, label = term), size = 3, colour = "black", hjust = 1, vjust = 1) +
    facet_grid(trait_fancy ~ Taxon, scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "top",
          plot.title = element_text(size = 10),
          strip.text.x = element_text(face = "italic"),
          aspect.ratio = 0.7) +
    guides(linetype = "none")

  # bryo_plot_leaf_eco <- b_dat |>
  #   filter(class == "Leaf economic") |>
  #   ggplot(aes(x = Elevation_m, y = value_trans, colour = Gradient, linetype = pvalue)) +
  #   geom_point(alpha = 0.5) +
  #   geom_smooth(mapping = aes(fill = Gradient), b_dat |> filter(class == "Leaf economic", pvalue == "sign"), method = "lm", se = TRUE) +
  #   scale_colour_manual(name = "", values = c("green4", "grey"), labels = c("Nutrient input", "Reference")) +
  #   scale_fill_manual(name = "", values = c("green4", "grey"), labels = c("Nutrient input", "Reference")) +
  #   labs(x = "Elevation in m a.s.l.", y = "Trait value") +
  #   geom_text(ata = b_dat |> distinct(Taxon, trait_fancy, term) , aes(x = Inf, y = Inf, label = term), size = 3, colour = "black", hjust = 1, vjust = 1) +
  #   facet_grid(trait_fancy ~ Taxon, scales = "free_y") +
  #   theme_minimal() +
  #   theme(legend.position = "top",
  #         plot.title = element_text(size = 10),
  #         strip.text.x = element_text(face = "italic"),
  #         aspect.ratio = 0.7) +
  #   guides(linetype = "none")
  #
  # bryo_plot_size / bryo_plot_leaf_eco +
  #   plot_layout(guides = 'collect',
  #               widths = c(1, 5),
  #               heights = c(1, 4)) &
  #   theme(legend.position = 'top')

  return(bryo_plot)
}


# ind_traits %>% distinct(Taxon)
# ind_traits %>% distinct(trait_trans)
#
# # manual model selection
# dd <- ind_traits %>%
#   filter(trait_trans == "dN15_permil",
#          Taxon == "poa arctica")
# fit <- lmer(value_trans ~ Gradient * Elevation_m + (1|GS), na.action = "na.fail", REML = FALSE, data = dd)
# model.set = dredge(fit, rank = "AICc", extra = "R^2")
# model.set

# Results from model selection
# sal.pol A: G+E, LDMC: G+E, N: G, N15: G
# #cer.arc LDMC: Null (removed!)
# luz.conf H: G+E, A: G+E, LDMC: GxE, N: G, N15: G
#
# Bryos: too little data and most are Null model!!!
# san.unc WHC: G
# aul.tur P: Null
# hyl spl SL: Null, SSL: Null




