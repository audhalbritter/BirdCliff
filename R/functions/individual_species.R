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

make_ind_sp_plot <- function(ind_traits, bryo_trait_output){

  # best model from model selection
  # best_model <- ind_traits %>%
  #   filter(Functional_group == "vascular") %>%
  #   distinct(Taxon, trait_trans) %>%
  #   mutate(best = c("", "Null", "G+E", "Null", "G", "G+E", "G+E", "GxE", "G", "Gx"),
  #          x.pos = c(rep(10, 10)),
  #          y.pos = c(2.8, 0.8, 0.5, 3.8, 12, 2.8, 0.8, 0.5, 3.8, 12))
  #
  #
  # v <- fancy_trait_name_dictionary(ind_traits) %>% distinct(trait_fancy)
  #   filter(Functional_group == "vascular") %>%
  #   left_join(best_model, by = c("Taxon", "trait_trans")) %>%
  #   ggplot(aes(x = Elevation_m, y = value_trans, colour = Gradient)) +
  #   geom_point(alpha = 0.5) +
  #   geom_smooth(method = "lm", se = FALSE) +
  #   #geom_text(aes(label = best, x = x.pos, y = y.pos), parse = TRUE, hjust = 0, size = 3, colour = "black") +
  #   scale_colour_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
  #   facet_grid(trait_fancy ~ Taxon, scales = "free_y") +
  #   labs(x = "", y = "") +
  #   theme_minimal() +
  #   theme(aspect.ratio = 0.8,
  #         plot.margin = margin(0.5, 0.5, 0.5, 0.5))


  #### VASCULAR PLANTS PLOT
  v_dat <- fancy_trait_name_dictionary(ind_traits) %>%
    filter(Functional_group == "vascular")

  ### Salix - height
  g0 <- v_dat |>
    filter(Taxon == "salix polaris",
           trait_trans == "Plant_Height_cm_log") |>
    ggplot(aes(x = Elevation_m, y = value_trans, colour = Gradient)) +
    geom_point(alpha = 0.5) +
    scale_colour_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10),
          text = element_text(size = 8),
          aspect.ratio = 0.7)

  gs_height <- g0 +
    labs(title = expression(italic("Salix polaris"))) +
    #annotate("text", x = -Inf, y = Inf, label = "(b)", size = 3, hjust = 0, vjust = 1) +
    theme(axis.text.x = element_blank())


  ### Salix - area
  gs_area <- g0 %+% subset(v_dat, trait_trans == "Leaf_Area_cm2_log") +
    #annotate("text", x = -Inf, y = Inf, label = "(d)", size = 3, hjust = 0, vjust = 1) +
    theme(axis.text.x = element_blank())



  ### Salix - N_percent
  gs_N <- g0 %+% subset(v_dat, trait_trans == "N_percent")
    #annotate("text", x = -Inf, y = Inf, label = "(h)", size = 3, hjust = 0, vjust = 1)

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


  gs_ldmc <- g0 %+% subset(v_dat, trait_trans == "LDMC") +
    geom_line(data = s_ldmc) +
    geom_ribbon(data = s_ldmc, aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_fill_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
    annotate("text", x = Inf, y = Inf, label = "G+E", size = 3, hjust = 1, vjust = 1) +
    theme(axis.text.x = element_blank())

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
    geom_ribbon(data = s_dN15, aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_fill_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
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
    geom_ribbon(data = l_height, aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_fill_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
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
    geom_ribbon(data = l_area, aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_fill_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
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
    geom_ribbon(data = l_ldmc, aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_fill_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
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
    geom_ribbon(data = l_N, aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_fill_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
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
    geom_ribbon(data = l_dN15, aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_fill_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
    labs(y = "dN15 ‰") +
  annotate("text", x = Inf, y = Inf, label = "GxE", size = 3, hjust = 1, vjust = 1)


  ### BRYO FIGURE

  b_dat <- fancy_trait_name_dictionary(ind_traits) %>%
    filter(Functional_group == "bryophyte")

  ### Aulacomnium turgidum - length
  b0 <- b_dat |>
    filter(Taxon == "aulacomnium turgidum",
           trait_trans == "Shoot_Length_cm_log") |>
    ggplot(aes(x = Elevation_m, y = value_trans, colour = Gradient)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = TRUE, aes(fill = Gradient)) +
    scale_colour_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
    scale_fill_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10),
          text = element_text(size = 8),
          aspect.ratio = 0.7)


  gat_length <- b0 +
    labs(y = "Shoot length cm", title = expression(italic("Aulacomnium turgidum"))) +
    annotate("text", x = Inf, y = Inf, label = "E", size = 3, hjust = 1, vjust = 1) +
    theme(axis.text.x = element_blank())


  ### Aulacomnium turgidum - ratio
  gat_ratio <- b0 %+% subset(b_dat, trait_trans == "Shoot_ratio") +
    labs(y = "Shoot ratio") +
    theme(axis.text.x = element_blank())

  ### Aulacomnium turgidum - WHC_g_g
  gat_WHC <- b0 %+% subset(b_dat, trait_trans == "WHC_g_g") +
    labs(y = "WHC g/g") +
    annotate("text", x = Inf, y = Inf, label = "GxE", size = 3, hjust = 1, vjust = 1) +
    theme(axis.text.x = element_blank())

  ### Aulacomnium turgidum - SSL_cm_g
  gat_SSL <- b0 %+% subset(b_dat, trait_trans == "SSL_cm_g") +
    labs(y = "SSL cm/g") +
    annotate("text", x = Inf, y = Inf, label = "GxE", size = 3, hjust = 1, vjust = 1) +
    theme(axis.text.x = element_blank())

  ### Aulacomnium turgidum - P_percent
  gat_P <- b0 %+% subset(b_dat, trait_trans == "P_percent") +
    labs(y = "P %") +
    annotate("text", x = Inf, y = Inf, label = "GxE", size = 3, hjust = 1, vjust = 1)



  ### hylocomium splendens - length
  ghs_length <- b0 %+% subset(b_dat, Taxon == "hylocomium splendens") +
    labs(title = expression(italic("Hylocomium splendens"))) +
    theme(axis.text.x = element_blank())

  ### hylocomium splendens - ratio
  ghs_ratio <- b0 %+% subset(b_dat, Taxon == "hylocomium splendens" &
                                trait_trans == "Shoot_ratio") +
    theme(axis.text.x = element_blank())

  ### hylocomium splendens - WHC_g_g
  ghs_WHC <- b0 %+% subset(b_dat, Taxon == "hylocomium splendens" &
                               trait_trans == "WHC_g_g") +
    theme(axis.text.x = element_blank())

  ### hylocomium splendens - SSL_cm_g
  ghs_SSL <- b0 %+% subset(b_dat, Taxon == "hylocomium splendens" &
                               trait_trans == "SSL_cm_g") +
    annotate("text", x = Inf, y = Inf, label = "G", size = 3, hjust = 1, vjust = 1) +
    theme(axis.text.x = element_blank())

  ### hylocomium splendens - P_percent
  ghs_P <- b0 %+% subset(b_dat, Taxon == "hylocomium splendens" &
                               trait_trans == "P_percent") +
    annotate("text", x = Inf, y = Inf, label = "GxE", size = 3, hjust = 1, vjust = 1)



  ### SANONIA SP - length
  gss_length <- b0 %+% subset(b_dat, Taxon == "sanionia sp") +
    labs(title = expression(italic("Sanionia sp"))) +
    theme(axis.text.x = element_blank())

  ### hylocomium splendens - ratio
  gss_ratio <- b0 %+% subset(b_dat, Taxon == "sanionia sp" &
                               trait_trans == "Shoot_ratio") +
    theme(axis.text.x = element_blank())

  ### hylocomium splendens - WHC_g_g
  gss_WHC <- b0 %+% subset(b_dat, Taxon == "sanionia sp" &
                             trait_trans == "WHC_g_g") +
    theme(axis.text.x = element_blank())

  ### hylocomium splendens - SSL_cm_g
  gss_SSL <- b0 %+% subset(b_dat, Taxon == "sanionia sp" &
                             trait_trans == "SSL_cm_g") +
    theme(axis.text.x = element_blank())

  ### hylocomium splendens - P_percent
  gss_P <- b0 %+% subset(b_dat, Taxon == "sanionia sp" &
                           trait_trans == "P_percent")

  library(cowplot)
  legend <- cowplot::get_legend(gss_SSL + theme(legend.position = "bottom"))


layout <- "
  ABCDE
  FGHIJ
  KLMNO
  PQRST
  UVWXY
  ZZZZZ
"

ind_sp_traits <- wrap_plots(gl_height, gs_height, gat_length, ghs_length, gss_length,
           gl_area, gs_area, gat_ratio, ghs_ratio, gss_ratio,
           gl_ldmc, gs_ldmc, gat_WHC, ghs_WHC, gss_WHC,
           gl_N, gs_N, gat_SSL, ghs_SSL, gss_SSL,
           gl_dN15, gs_dN15, gat_P, ghs_P, gss_P,
           legend) +
  plot_layout(design = layout)

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

