

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


# model.sel <- trait_mean %>%
#   filter(Gradient == "B") %>%
#   mutate(Site2 = if_else(Site == "5", "high", "low")) |>
#   group_by(trait_trans) %>%
#   nest(data = -c(trait_trans)) %>%
#   mutate(model = map(data, ~{
#     mod <- lmer(mean ~  Site2 + (1|GS), data = .x)
#   })) %>%
#   unnest(model)


# to do
#sum total divide each blue and yellow by 100
# x axis title: relative contribution
# table with all values
# split ind specie figure in 2
# non parametric test for proportions, is mean trait value iTV different N traits vs growth traits
# Climate figure: soil temp lines with G model
# check shoot length log transformation





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

  return(list(estimate, r_squared))

}





