

make_trait_model_selection <- function(trait_mean){

  model.sel <- trait_mean %>%
    # remove singular fit
    #filter(!trait_trans %in% c("SLA_cm2_g")) %>%
    #filter(!trait_trans %in% c("SLA_cm2_g", "dC13_permil")) %>%
    group_by(trait_trans) %>%
    nest(data = -c(trait_trans)) %>%
    mutate(model.set = map(data, ~{
      mod <- lmer(mean ~  Gradient * Elevation_m + (1|GS), REML = FALSE, na.action = "na.fail", data = .x)
      model.set = dredge(mod, rank = "AICc", extra = "R^2")
    })) %>%
    unnest(model.set)

  return(model.sel)

}

# Model selection
#library(MuMIn)
# dd <- trait_mean %>%
#   filter(trait_trans == "SLA_cm2_g")
# fit1 <- lmer(mean ~ Gradient * Elevation_m + (1|GS), REML = TRUE, na.action = "na.fail", data = dd)
# model.set <- dredge(fit1, rank = "AICc", extra = "R^2")
# model.set
#
# model.sel %>% filter(trait_trans == "Thickness_mm_log")


test_top_site <- function(trait_mean){

  top_site <- trait_mean %>%
    filter(Gradient == "B") %>%
    mutate(Site2 = if_else(Site == "5", "high", "low"),
           Site2 = factor(Site2, levels = c("low", "high"))) %>%
    group_by(trait_trans) %>%
    nest(data = -c(trait_trans)) %>%
    mutate(best = map(data, ~{
      mod1 = lmer(mean ~  Site2 + (1|GS), data = .)
      mod2 = lmer(mean ~  1 + (1|GS), data = .)
      best = anova(mod1, mod2)
        })) %>%
    unnest(best)

  fancy_trait_name_dictionary(top_site) |>
    ungroup() |>
    mutate(term = if_else(npar == 4, "Site", "Null"),
           AIC = round(AIC, digits = 1),
           `Pr(>Chisq)` = round(`Pr(>Chisq)`, digits = 3)) |>
    select(trait_fancy, term, AIC, df = Df, "p value" = `Pr(>Chisq)`) %>%
    write_csv(., file = "output/top_site.csv")

  return(top_site)
}


# dd <- trait_mean %>%
#   filter(Gradient == "B") %>%
#   mutate(Site2 = if_else(Site == "5", "high", "low"),
#          Site2 = factor(Site2, levels = c("low", "high"))) %>%
#   filter(trait_trans == "CN_ratio")
# fit <- lmer(mean ~  Site2 + (1|GS), data = dd)
# fit2 <- lmer(mean ~  1 + (1|GS), data = dd)
# best <- anova(fit, fit2)
# tidy(best)
# performance::check_model(fit)



# to do
# non parametric test for proportions, is mean trait value iTV different N traits vs growth traits




make_trait_output <- function(trait_mean){

  # NULL MODEL
  # estimate
  null_est <- trait_mean %>%
    filter(trait_trans %in% c("C_percent", "Dry_Mass_g_log", "P_percent", "Thickness_mm_log", "NP_ratio")) %>%
    group_by(trait_trans) %>%
    nest(data = -c(trait_trans)) %>%
    mutate(estimate = map(data, ~{
      mod <- lmer(mean ~ 1 + (1|Site), data = .x)
      estimates = broom.mixed::tidy(mod)
    })) %>%
    unnest(estimate)


  # r squared
  null_r <- trait_mean %>%
    filter(trait_trans %in% c("C_percent", "Dry_Mass_g_log", "P_percent", "Thickness_mm_log", "NP_ratio")) %>%
    group_by(trait_trans) %>%
    nest(data = -c(trait_trans)) %>%
    mutate(r = map(data, ~{
      mod <- lmer(mean ~ 1 + (1|Site), data = .x)
      r = as.numeric(r.squaredGLMM(mod))
    })) %>%
    unnest_wider(col = r) %>%
    select(trait_trans, "Rm" = "...1", "Rc" = "...2")


  # ELEVATION MODEL
  # estimate
  e_est <- trait_mean %>%
    filter(trait_trans %in% c("Plant_Height_cm_log")) %>%
    group_by(trait_trans) %>%
    nest(data = -c(trait_trans)) %>%
    mutate(estimate = map(data, ~{
      mod <- lmer(mean ~ Elevation_m + (1|Site), data = .x)
      estimates = broom.mixed::tidy(mod)
    })) %>%
    unnest(estimate)


  # r squared
  e_r <- trait_mean %>%
    filter(trait_trans %in% c("Plant_Height_cm_log")) %>%
    group_by(trait_trans) %>%
    nest(data = -c(trait_trans)) %>%
    mutate(r = map(data, ~{
      mod <- lmer(mean ~ Elevation_m + (1|Site), data = .x)
      r = as.numeric(r.squaredGLMM(mod))
    })) %>%
    unnest_wider(col = r) %>%
    select(trait_trans, "Rm" = "...1", "Rc" = "...2")


  # GRADIENT MODEL
  # estimate
  g_est <- trait_mean %>%
    filter(trait_trans %in% c("CN_ratio", "N_percent")) %>%
    group_by(trait_trans) %>%
    nest(data = -c(trait_trans)) %>%
    mutate(estimate = map(data, ~{
      mod <- lmer(mean ~ Gradient + (1|Site), data = .x)
      estimates = broom.mixed::tidy(mod)
    })) %>%
    unnest(estimate)


  # r squared
  g_r <- trait_mean %>%
    filter(trait_trans %in% c("CN_ratio", "N_percent")) %>%
    group_by(trait_trans) %>%
    nest(data = -c(trait_trans)) %>%
    mutate(r = map(data, ~{
      mod <- lmer(mean ~ Gradient + (1|Site), data = .x)
      r = as.numeric(r.squaredGLMM(mod))
    })) %>%
    unnest_wider(col = r) %>%
    select(trait_trans, "Rm" = "...1", "Rc" = "...2")


  # GRADIENT + ELEVATION MODEL
  # estimate
  ge_est <- trait_mean %>%
    filter(trait_trans %in% c("Leaf_Area_cm2_log", "LDMC")) %>%
    group_by(trait_trans) %>%
    nest(data = -c(trait_trans)) %>%
    mutate(estimate = map(data, ~{
      mod <- lmer(mean ~ Gradient + Elevation_m + (1|Site), data = .x)
      estimates = broom.mixed::tidy(mod)
    })) %>%
    unnest(estimate)


  # r squared
  ge_r <- trait_mean %>%
    filter(trait_trans %in% c("Leaf_Area_cm2_log", "LDMC")) %>%
    group_by(trait_trans) %>%
    nest(data = -c(trait_trans)) %>%
    mutate(r = map(data, ~{
      mod <- lmer(mean ~ Gradient + Elevation_m + (1|Site), data = .x)
      r = as.numeric(r.squaredGLMM(mod))
    })) %>%
    unnest_wider(col = r) %>%
    select(trait_trans, "Rm" = "...1", "Rc" = "...2")


  # GRADIENT * ELEVATION MODEL
  # estimate
  gxe_est <- trait_mean %>%
    filter(trait_trans %in% c("dN15_permil")) %>%
    group_by(trait_trans) %>%
    nest(data = -c(trait_trans)) %>%
    mutate(estimate = map(data, ~{
      mod <- lmer(mean ~ Gradient * Elevation_m + (1|Site), data = .x)
      estimates = broom.mixed::tidy(mod)
    })) %>%
    unnest(estimate)


  # r squared
  gxe_r <- trait_mean %>%
    filter(trait_trans %in% c("dN15_permil")) %>%
    group_by(trait_trans) %>%
    nest(data = -c(trait_trans)) %>%
    mutate(r = map(data, ~{
      mod <- lmer(mean ~ Gradient * Elevation_m + (1|Site), data = .x)
      r = as.numeric(r.squaredGLMM(mod))
    })) %>%
    unnest_wider(col = r) %>%
    select(trait_trans, "Rm" = "...1", "Rc" = "...2")

  estimate <- bind_rows(
    Null = null_est,
    E = e_est,
    G = g_est,
    "G+E" = ge_est,
    GxE = gxe_est,
    .id = "Model"
    )

  r_squared <- bind_rows(
    Null = null_r,
    E = e_r,
    G = g_r,
    "G+E" = ge_r,
    GxE = gxe_r,
    .id = "Model"
  )

  # save regression model outputs
  fancy_trait_name_dictionary(estimate) %>%
    filter(effect == "fixed") %>%
    select(Trait = trait_fancy, Model, term, estimate:statistic) %>%
    mutate(term = recode(term, "(Intercept)" = "Intercept", "Elevation_m" = "E", "GradientC" = "G", "GradientC:Elevation_m" = "GxE")) %>%
    left_join(r_squared %>%
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


  return(list(estimate, r_squared))

}

