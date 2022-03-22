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


#### VASCULAR PLANTS

# Model selection (DOES NOT RUN!!!)
# run_ind_model_sel <- function(ind_traits){
#
#   ind_traits %>%
#     filter(Functional_group == "vascular", Taxon == "salix polaris") %>%
#     group_by(Taxon, trait_trans) %>%
#     mutate(model.set = map(data, ~{
#       mod <- lmer(Value ~  Gradient * Elevation_m + (1|GS), REML = FALSE, na.action = "na.fail", data = .x)
#       model.set = dredge(mod, rank = "AICc", extra = "R^2")
#     })) %>%
#     unnest(model.set)
#
# }


# Run models for single species
run_vascular_plant_models <- function(ind_traits){

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

  # save output
  # best model
  best_ind_model <- ind_traits %>%
    filter(Functional_group == "vascular") %>%
    distinct(Taxon, trait_trans) %>%
    mutate(best = c("NA", "Null", "G+E", "Null", "G", "G+E", "G+E", "GxE", "G", "Gx"))

  fancy_trait_name_dictionary(estimate) %>%
    filter(effect == "fixed") %>%
    left_join(best_ind_model, by = c("Taxon", "trait_trans")) %>%
    select(Trait = trait_fancy, Model = best, term, estimate:statistic) %>%
    mutate(term = recode(term, "(Intercept)" = "Intercept", "Elevation_m" = "E", "GradientC" = "G", "GradientC:Elevation_m" = "GxE")) %>%
    left_join(r %>%
                select(trait_trans, Rm, Rc), by = c("trait_trans", "Taxon")) %>%
    mutate(estimate = round(estimate, digits = 2),
           std.error = round(std.error, digits = 2),
           statistic = round(statistic, digits = 2),
           Rm = round(Rm, digits = 2),
           Rc = round(Rc, digits = 2)) %>%
    ungroup() %>%
    select(-trait_trans, Term = term, Estimate = estimate, "Std error" = std.error, "t-value" = "statistic", "Marginal R2" = Rm, "Conditional R2" = Rc) %>%
    arrange(Taxon, Trait) %>%
    write_csv(., file = "output/Ind_sp_regression_output.csv")

  return(list(estimate, r))
}




make_ind_vascular_plant_plot <- function(ind_traits){

  #### VASCULAR PLANTS PLOT
  v_dat <- fancy_trait_name_dictionary(ind_traits) %>%
    filter(Functional_group == "vascular")

  ### Salix - LDMC
  dd <- v_dat %>%
    filter(Taxon == "salix polaris",
           trait_trans == "LDMC")
  fit <- lmer(value_trans ~ Gradient + Elevation_m + (1|Site), data = dd)

  newdat <- dd %>%
    distinct(Elevation_m, Gradient) %>%
    mutate(value_trans = 0)
  newdat$value_trans <-  predict(fit, newdat, re.form = NA)

  mm <- model.matrix(terms(fit), newdat)

  s_ldmc <- newdat %>%
    mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
           tvar1 = pvar1 + VarCorr(fit)$Site[1],  ## must be adapted for more complex models
           cmult = 1.96) %>%
    mutate(plo = value_trans - cmult*sqrt(pvar1),
           phi = value_trans + cmult*sqrt(pvar1),
           tlo = value_trans - cmult*sqrt(tvar1),
           thi = value_trans + cmult*sqrt(tvar1))

  g0 <- v_dat |>
    filter(Taxon == "salix polaris",
           trait_trans == "LDMC") |>
    ggplot(aes(x = Elevation_m, y = value_trans, colour = Gradient, fill = Gradient)) +
    geom_point(alpha = 0.5) +
    scale_colour_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
    scale_fill_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 10),
          text = element_text(size = 8),
          aspect.ratio = 0.7)

  gs_ldmc <- g0 +
    geom_line(data = s_ldmc) +
    geom_ribbon(data = s_ldmc, aes(ymin = plo, ymax = phi), alpha = 0.3, linetype = 0) +
    annotate("text", x = Inf, y = Inf, label = "G+E", size = 3, hjust = 1, vjust = 1) +
    theme(axis.text.x = element_blank())


  ### Salix - height
  fake_data <- s_ldmc |>
    slice(0)

  gs_height <- g0 %+% subset(v_dat, trait_trans == "Plant_Height_cm_log") +
    geom_line(data = fake_data) +
    geom_ribbon(data = fake_data, aes(ymin = plo, ymax = phi), alpha = 0.3, linetype = 0) +
    labs(title = expression(italic("Salix polaris"))) +
    #annotate("text", x = -Inf, y = Inf, label = "(b)", size = 3, hjust = 0, vjust = 1) +
    theme(axis.text.x = element_blank())


  ### Salix - area
  gs_area <- g0 %+% subset(v_dat, trait_trans == "Leaf_Area_cm2_log") +
    geom_line(data = fake_data) +
    geom_ribbon(data = fake_data, aes(ymin = plo, ymax = phi), alpha = 0.3, linetype = 0) +
    #annotate("text", x = -Inf, y = Inf, label = "(d)", size = 3, hjust = 0, vjust = 1) +
    annotate("text", x = Inf, y = Inf, label = "Null", size = 3, hjust = 1, vjust = 1) +
    theme(axis.text.x = element_blank())



  ### Salix - N_percent
  gs_N <- g0 %+% subset(v_dat, trait_trans == "N_percent") +
    geom_line(data = fake_data) +
    geom_ribbon(data = fake_data, aes(ymin = plo, ymax = phi), alpha = 0.3, linetype = 0) +
    annotate("text", x = Inf, y = Inf, label = "Null", size = 3, hjust = 1, vjust = 1)
    #annotate("text", x = -Inf, y = Inf, label = "(h)", size = 3, hjust = 0, vjust = 1)



  ### Salix - dN15
  dd <- v_dat %>%
    filter(Taxon == "salix polaris",
           trait_trans == "dN15_permil")
  fit <- lmer(value_trans ~ Gradient + (1|Site), data = dd)

  newdat <- dd %>%
    distinct(Gradient, Elevation_m) %>%
    mutate(value_trans = 0)
  newdat$value_trans <-  predict(fit, newdat, re.form = NA)

  mm <- model.matrix(terms(fit), newdat)

  s_dN15 <- newdat %>%
    mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
           tvar1 = pvar1 + VarCorr(fit)$Site[1],  ## must be adapted for more complex models
           cmult = 1.96) %>%
    mutate(plo = value_trans - cmult*sqrt(pvar1),
           phi = value_trans + cmult*sqrt(pvar1),
           tlo = value_trans - cmult*sqrt(tvar1),
           thi = value_trans + cmult*sqrt(tvar1))


  gs_dN15 <- g0 %+% subset(v_dat, trait_trans == "dN15_permil") +
    geom_line(data = s_dN15) +
    geom_ribbon(data = s_dN15, aes(ymin = plo, ymax = phi), alpha = 0.3, linetype = 0) +
    annotate("text", x = Inf, y = Inf, label = "G", size = 3, hjust = 1, vjust = 1)


  ### LUZULA - Plant_Height_cm_log
  dd <- v_dat %>%
    filter(Taxon == "luzula confusa",
           trait_trans == "Plant_Height_cm_log")
  fit <- lmer(value_trans ~ Gradient + Elevation_m + (1|Site), data = dd)

  newdat <- dd %>%
    distinct(Elevation_m, Gradient) %>%
    mutate(value_trans = 0)
  newdat$value_trans <-  predict(fit, newdat, re.form = NA)

  mm <- model.matrix(terms(fit), newdat)

  l_height <- newdat %>%
    mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
           tvar1 = pvar1 + VarCorr(fit)$Site[1],  ## must be adapted for more complex models
           cmult = 1.96) %>%
    mutate(plo = value_trans - cmult*sqrt(pvar1),
           phi = value_trans + cmult*sqrt(pvar1),
           tlo = value_trans - cmult*sqrt(tvar1),
           thi = value_trans + cmult*sqrt(tvar1))


  gl_height <- g0 %+% subset(v_dat, trait_trans == "Plant_Height_cm_log") +
    geom_line(data = l_height) +
    geom_ribbon(data = l_height, aes(ymin = plo, ymax = phi), alpha = 0.3, linetype = 0) +
    labs(y = "Height cm", title = expression(italic("Luzula confusa"))) +
    annotate("text", x = Inf, y = Inf, label = "G+E", size = 3, hjust = 1, vjust = 1) +
    theme(axis.text.x = element_blank())


  ### luzula - Leaf_Area_cm2_log
  dd <- v_dat %>%
    filter(Taxon == "luzula confusa",
           trait_trans == "Leaf_Area_cm2_log")
  fit <- lmer(value_trans ~ Gradient + Elevation_m + (1|Site), data = dd)

  newdat <- dd %>%
    distinct(Elevation_m, Gradient) %>%
    mutate(value_trans = 0)
  newdat$value_trans <-  predict(fit, newdat, re.form = NA)

  mm <- model.matrix(terms(fit), newdat)

  l_area <- newdat %>%
    mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
           tvar1 = pvar1 + VarCorr(fit)$Site[1],  ## must be adapted for more complex models
           cmult = 1.96) %>%
    mutate(plo = value_trans - cmult*sqrt(pvar1),
           phi = value_trans + cmult*sqrt(pvar1),
           tlo = value_trans - cmult*sqrt(tvar1),
           thi = value_trans + cmult*sqrt(tvar1))


  gl_area <- g0 %+% subset(v_dat, trait_trans == "Leaf_Area_cm2_log") +
    geom_line(data = l_area) +
    geom_ribbon(data = l_area, aes(ymin = plo, ymax = phi), alpha = 0.3, linetype = 0) +
    labs(y = "Area cm2") +
    annotate("text", x = Inf, y = Inf, label = "G+E", size = 3, hjust = 1, vjust = 1) +
    theme(axis.text.x = element_blank())



  ### luzula - LDMC
  dd <- v_dat %>%
    filter(Taxon == "luzula confusa",
           trait_trans == "LDMC")
  fit <- lmer(value_trans ~ Gradient * Elevation_m + (1|Site), data = dd)

  newdat <- dd %>%
    distinct(Elevation_m, Gradient) %>%
    mutate(value_trans = 0)
  newdat$value_trans <-  predict(fit, newdat, re.form = NA)

  mm <- model.matrix(terms(fit), newdat)

  l_ldmc <- newdat %>%
    mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
           tvar1 = pvar1 + VarCorr(fit)$Site[1],  ## must be adapted for more complex models
           cmult = 1.96) %>%
    mutate(plo = value_trans - cmult*sqrt(pvar1),
           phi = value_trans + cmult*sqrt(pvar1),
           tlo = value_trans - cmult*sqrt(tvar1),
           thi = value_trans + cmult*sqrt(tvar1))


  gl_ldmc <- g0 %+% subset(v_dat, trait_trans == "LDMC") +
    geom_line(data = l_ldmc) +
    geom_ribbon(data = l_ldmc, aes(ymin = plo, ymax = phi), alpha = 0.3, linetype = 0) +
    labs(y = "LDMC") +
    annotate("text", x = Inf, y = Inf, label = "GxE", size = 3, hjust = 1, vjust = 1) +
    theme(axis.text.x = element_blank())


  ### luzula - N
  dd <- v_dat %>%
    filter(Taxon == "luzula confusa",
           trait_trans == "N_percent")
  fit <- lmer(value_trans ~ Gradient + (1|Site), data = dd)

  newdat <- dd %>%
    distinct(Elevation_m, Gradient) %>%
    mutate(value_trans = 0)
  newdat$value_trans <-  predict(fit, newdat, re.form = NA)

  mm <- model.matrix(terms(fit), newdat)

  l_N <- newdat %>%
    mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
           tvar1 = pvar1 + VarCorr(fit)$Site[1],  ## must be adapted for more complex models
           cmult = 1.96) %>%
    mutate(plo = value_trans - cmult*sqrt(pvar1),
           phi = value_trans + cmult*sqrt(pvar1),
           tlo = value_trans - cmult*sqrt(tvar1),
           thi = value_trans + cmult*sqrt(tvar1))


  gl_N <- g0 %+% subset(v_dat, trait_trans == "N_percent") +
    geom_line(data = l_N) +
    geom_ribbon(data = l_N, aes(ymin = plo, ymax = phi), alpha = 0.3, linetype = 0) +
    labs(y = "N %") +
    annotate("text", x = Inf, y = Inf, label = "G", size = 3, hjust = 1, vjust = 1) +
    theme(axis.text.x = element_blank())


  ### luzula - dN15_permil
  dd <- v_dat %>%
    filter(Taxon == "luzula confusa",
           trait_trans == "dN15_permil")
  fit <- lmer(value_trans ~ Gradient * Elevation_m + (1|Site), data = dd)

  newdat <- dd %>%
    distinct(Elevation_m, Gradient) %>%
    mutate(value_trans = 0)
  newdat$value_trans <-  predict(fit, newdat, re.form = NA)

  mm <- model.matrix(terms(fit), newdat)

  l_dN15 <- newdat %>%
    mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
           tvar1 = pvar1 + VarCorr(fit)$Site[1],  ## must be adapted for more complex models
           cmult = 1.96) %>%
    mutate(plo = value_trans - cmult*sqrt(pvar1),
           phi = value_trans + cmult*sqrt(pvar1),
           tlo = value_trans - cmult*sqrt(tvar1),
           thi = value_trans + cmult*sqrt(tvar1))


  gl_dN15 <- g0 %+% subset(v_dat, trait_trans == "dN15_permil") +
    geom_line(data = l_dN15) +
    geom_ribbon(data = l_dN15, aes(ymin = plo, ymax = phi), alpha = 0.3, linetype = 0) +
    labs(y = "δN15 ‰") +
  annotate("text", x = Inf, y = Inf, label = "GxE", size = 3, hjust = 1, vjust = 1)


  ind_sp_traits <- wrap_plots(gl_height, gs_height,
             gl_area, gs_area,
             gl_ldmc, gs_ldmc,
             gl_N, gs_N,
             gl_dN15, gs_dN15, ncol = 2) +
    plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

  return(ind_sp_traits)

}


### BROPHYTES


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
    select(Taxon, Trait = trait_fancy, Term = term, "Std. error" = std.error, "t value" = statistic, "p value" = p.value) %>%
    arrange(Taxon, Trait)

  write_csv(bryo_trait_output, file = "output/bryo_trait_output.csv")

  return(bryo_trait_output)

}



make_bryo_figure <- function(ind_traits, bryo_trait_output){

  ### BRYO FIGURE

  term = tibble(term = c("E", "Null", "GxE", "GxE", "GxE", "Null", "Null", "G", "Null", "GxE", "Null", "Null", "Null", "Null", "Null"))

  result <- bryo_trait_output |>
    distinct(Taxon, Trait) |>
   bind_cols(term = term) |>
    mutate(pvalue = if_else(term == "Null", "non-sign", "sign"))

  b_dat <- fancy_trait_name_dictionary(ind_traits) %>%
    filter(Functional_group == "bryophyte") |>
    left_join(result, by = c("Taxon", "trait_fancy" = "Trait")) |>
    mutate(pvalue = factor(pvalue, levels = c("sign", "non-sign")),
           Taxon = recode(Taxon, "sanionia sp" = "Sanionia sp", "hylocomium splendens" = "Hylocomium splendens", "aulacomnium turgidum" = "Aulacomnium turgidum"))

  bryo_plot <- b_dat |>
    ggplot(aes(x = Elevation_m, y = value_trans, colour = Gradient, linetype = pvalue)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = TRUE, aes(fill = Gradient)) +
    scale_colour_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
    scale_fill_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
    labs(x = "Elevation in m a.s.l.", y = "Trait value") +
    geom_text(ata = b_dat |> distinct(Taxon, trait_fancy, term) , aes(x = Inf, y = Inf, label = term), size = 3, colour = "black", hjust = 1, vjust = 1) +
    facet_grid(trait_fancy ~ Taxon, scales = "free_y") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10),
          strip.text.x = element_text(face = "italic"),
          aspect.ratio = 0.7) +
    guides(linetype = "none")

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




