# traits_raw %>% distinct(Taxon)
#
# # choosing species and traits
# traits_raw %>%
#   distinct(Gradient, Site, Taxon) %>%
#   group_by(Gradient, Taxon) %>%
#   count() %>%
#   filter(n > 3)
#
# bryo_traits_raw %>%
#   distinct(Gradient, Site, Taxon) %>%
#   group_by(Gradient, Taxon) %>%
#   count()


combine_traits <- function(traits_raw, bryo_traits_raw){

  # cobmine species vascular and bryophytes
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

make_ind_sp_plot <- function(ind_traits){

  # best model from model selection
  best_model <- ind_traits %>%
    filter(Functional_group == "vascular") %>%
    distinct(Taxon, trait_trans) %>%
    mutate(best = c("NA", "Null", "G+E", "Null", "G", "G+E", "G+E", "GxE", "G", "Gx"),
           x.pos = c(rep(10, 10)),
           y.pos = c(2.8, 0.8, 0.5, 3.8, 12, 2.8, 0.8, 0.5, 3.8, 12))


  v <- fancy_trait_name_dictionary(ind_traits) %>%
    filter(Functional_group == "vascular") %>%
    left_join(best_model, by = c("Taxon", "trait_trans")) %>%
    ggplot(aes(x = Elevation_m, y = value_trans, colour = Gradient)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_text(aes(label = best, x = x.pos, y = y.pos), parse = TRUE, hjust = 0, size = 3, colour = "black") +
    scale_colour_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
    facet_grid(trait_fancy ~ Taxon, scales = "free_y") +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(aspect.ratio = 0.8,
          plot.margin = margin(0.5, 0.5, 0.5, 0.5))


  b <- fancy_trait_name_dictionary(ind_traits) %>%
    filter(Functional_group == "bryophyte") %>%
    ggplot(aes(x = Elevation_m, y = value_trans, colour = Gradient)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE) +
    scale_colour_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
    facet_grid(trait_fancy ~ Taxon, scales = "free_y") +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(aspect.ratio = 0.8,
          plot.margin = margin(0.5, 0.5, 0.5, 0.5))

  ind_sp_traits <- v + b + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

  return(ind_sp_traits)

}



make_bryo_trait_model <- function(){
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
    select(Taxon, Trait = trait_fancy, Term = term, "Std. error" = std.error, "t value" = statistic, "p value" = p.value) %>%
    arrange(Taxon, Trait)

  return(bryo_trait_output)

}


# Model selection (DOES NOT RUN!!!)
run_ind_model_sel <- function(ind_traits){

  ind_traits %>%
    filter(Taxon == "salix polaris") %>%
    group_by(trait_trans) %>%
    nest(data = -c(trait_trans)) %>%
    mutate(model.set = map(data, ~{
      mod <- lmer(mean ~  Gradient * Elevation_m + (1|GS), REML = FALSE, na.action = "na.fail", data = .x)
      model.set = dredge(mod, rank = "AICc", extra = "R^2")
      model.set
    })) %>%
    unnest(model.set)

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




# Run models for single species
run_ind_models <- function(ind_traits){

  ### NULL MODEL
  dat <- ind_traits %>%
    filter(Taxon == "salix polaris" &  trait_trans %in% c("Leaf_Area_cm2_log", "N_percent"))

  estimate_n <- dat %>%
    group_by(Taxon, trait_trans) %>%
    nest(data = -c(Taxon, trait_trans)) %>%
    mutate(estimate = map(data, ~{
      mod <- lmer(value_trans ~ 1 + (1|GS), data =  .x)
      estimates = broom.mixed::tidy(mod)
    })) %>%
    unnest(estimate)

  r_square_n <- dat %>%
    group_by(Taxon, trait_trans) %>%
    nest(data = -c(Taxon, trait_trans)) %>%
    mutate(r = map(data, ~{
      mod <- lmer(value_trans ~ 1 + (1|GS), data =  .x)
      r = as.numeric(r.squaredGLMM(mod))
    })) %>%
    unnest_wider(col = r) %>%
    select(trait_trans, "Rm" = "...1", "Rc" = "...2")


  ### GRADIENT MODEL
  dat <- ind_traits %>%
    filter(Taxon == "salix polaris" &  trait_trans == "dN15_permil" |
           Taxon == "luzula confusa" & trait_trans == "N_percent")

  estimate_g <- dat %>%
    group_by(Taxon, trait_trans) %>%
    nest(data = -c(Taxon, trait_trans)) %>%
    mutate(estimate = map(data, ~{
      mod <- lmer(value_trans ~ Gradient + (1|GS), data =  .x)
      estimates = broom.mixed::tidy(mod)
    })) %>%
    unnest(estimate)

  r_square_g <- dat %>%
    group_by(Taxon, trait_trans) %>%
    nest(data = -c(Taxon, trait_trans)) %>%
    mutate(r = map(data, ~{
      mod <- lmer(value_trans ~ Gradient + (1|GS), data =  .x)
      r = as.numeric(r.squaredGLMM(mod))
    })) %>%
    unnest_wider(col = r) %>%
    select(trait_trans, "Rm" = "...1", "Rc" = "...2")


  ### GRADIENT + ELEVATION MODEL
  dat <- ind_traits %>%
    filter(Taxon == "salix polaris" & trait_trans == "LDMC" |
             Taxon == "luzula confusa"  & trait_trans %in% c("Plant_Height_cm_log", "Leaf_Area_cm2_log"))

  estimate_ge <- dat %>%
    group_by(Taxon, trait_trans) %>%
    nest(data = -c(Taxon, trait_trans)) %>%
    mutate(estimate = map(data, ~{
      mod <- lmer(value_trans ~ Gradient + Elevation_m + (1|GS), data =  .x)
      estimates = broom.mixed::tidy(mod)
    })) %>%
    unnest(estimate)

  r_square_ge <- dat %>%
    group_by(Taxon, trait_trans) %>%
    nest(data = -c(Taxon, trait_trans)) %>%
    mutate(r = map(data, ~{
      mod <- lmer(value_trans ~ Gradient + Elevation_m + (1|GS), data =  .x)
      r = as.numeric(r.squaredGLMM(mod))
    })) %>%
    unnest_wider(col = r) %>%
    select(trait_trans, "Rm" = "...1", "Rc" = "...2")


  ### GRADIENT * ELEVATION MODEL
  dat <- ind_traits %>%
    filter(Taxon %in% c("luzula confusa")  & trait_trans %in% c("LDMC", "dN15_permil"))

  estimate_gxe <- dat %>%
    group_by(Taxon, trait_trans) %>%
    nest(data = -c(Taxon, trait_trans)) %>%
    mutate(estimate = map(data, ~{
      mod <- lmer(value_trans ~ Gradient * Elevation_m + (1|GS), data =  .x)
      estimates = broom.mixed::tidy(mod)
    })) %>%
    unnest(estimate)

  r_square_gxe <- dat %>%
    group_by(Taxon, trait_trans) %>%
    nest(data = -c(Taxon, trait_trans)) %>%
    mutate(r = map(data, ~{
      mod <- lmer(value_trans ~ Gradient * Elevation_m + (1|GS), data =  .x)
      r = as.numeric(r.squaredGLMM(mod))
    })) %>%
    unnest_wider(col = r) %>%
    select(trait_trans, "Rm" = "...1", "Rc" = "...2")


  estimate <- bind_rows(
    Null = estimate_n,
    G = estimate_g,
    GE = estimate_ge,
    GxE = estimate_gxe,
    .id = "Best_model"
  )

  r <- bind_rows(
    Null = r_square_n,
    G = r_square_g,
    GE = r_square_ge,
    GxE = r_square_gxe,
    .id = "Best_model"
  )

  return(list(estimate, r))
}

