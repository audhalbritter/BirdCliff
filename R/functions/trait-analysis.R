make_trait_analysis <- function(trait_mean){

  dat <- fancy_trait_name_dictionary(trait_mean)

  #C_percent
  dd <- dat %>%
    filter(trait_trans == "C_percent")
  fit <- lmer(mean ~ Gradient * Elevation_m + (1|Site), data = dd)

  newdat <- dd %>%
    distinct(Elevation_m, Gradient) %>%
    mutate(mean = predict(fit, newdat, re.form = NA))

  mm <- model.matrix(terms(fit), newdat)

  c <- newdat %>%
    mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
           tvar1 = pvar1 + VarCorr(fit)$Site[1],  ## must be adapted for more complex models
           cmult = 1.96) %>%
    mutate(plo = mean - cmult*sqrt(pvar1),
           phi = mean + cmult*sqrt(pvar1),
           tlo = mean - cmult*sqrt(tvar1),
           thi = mean + cmult*sqrt(tvar1))

  #plot confidence
  g0 <-  dat %>%
    filter(trait_trans == "C_percent") %>%
    ggplot(aes(x = Elevation_m, y = mean, colour = Gradient)) +
    geom_point(alpha = 0.5) +
    scale_colour_manual(values = c("green4", "grey"), labels = c("Birdcliff", "Reference")) +
    facet_wrap(~ trait_fancy, scales = "free_y") +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none",
          aspect.ratio = 0.7,
          plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5))

  gc <- g0 +
    geom_line(data = c) +
    geom_ribbon(data = c, aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_fill_manual(values = c("green4", "grey"), labels = c("Birdcliff", "Reference")) +
    theme(axis.text.x = element_blank())


  #CN_ratio
  dd <- dat %>%
    filter(trait_trans == "CN_ratio")
  fit <- lmer(mean ~ Gradient * Elevation_m + (1|Site), data = dd)

  newdat <- dd %>%
    distinct(Elevation_m, Gradient) %>%
    mutate(mean = predict(fit, newdat, re.form = NA))

  mm <- model.matrix(terms(fit), newdat)

  cn <- newdat %>%
    mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
           tvar1 = pvar1 + VarCorr(fit)$Site[1],  ## must be adapted for more complex models
           cmult = 1.96) %>%
    mutate(plo = mean - cmult*sqrt(pvar1),
           phi = mean + cmult*sqrt(pvar1),
           tlo = mean - cmult*sqrt(tvar1),
           thi = mean + cmult*sqrt(tvar1))

  gcn <- g0 %+% subset(dat, trait_trans == "CN_ratio") +
    geom_line(data = cn) +
    geom_ribbon(data = cn, aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_fill_manual(values = c("green4", "grey"), labels = c("Birdcliff", "Reference")) +
    annotate("text", x = 50, y = 24, label = "E", size = 3)




  #dN15_permil
  dd <- dat %>%
    filter(trait_trans == "dN15_permil")
  fit <- lmer(mean ~ Gradient * Elevation_m + (1|Site), data = dd)

  newdat <- dd %>%
    distinct(Elevation_m, Gradient) %>%
    mutate(mean = predict(fit, newdat, re.form = NA))

  mm <- model.matrix(terms(fit), newdat)

  dn <- newdat %>%
    mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
           tvar1 = pvar1 + VarCorr(fit)$Site[1],  ## must be adapted for more complex models
           cmult = 1.96) %>%
    mutate(plo = mean - cmult*sqrt(pvar1),
           phi = mean + cmult*sqrt(pvar1),
           tlo = mean - cmult*sqrt(tvar1),
           thi = mean + cmult*sqrt(tvar1))

  gdn <- g0 %+% subset(dat, trait_trans == "dN15_permil") +
    geom_line(data = dn) +
    geom_ribbon(data = dn, aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_fill_manual(values = c("green4", "grey"), labels = c("Birdcliff", "Reference")) +
    annotate("text", x = 50, y = 15, label = "GxE", size = 3)



  #Dry_Mass_g_log
  dd <- dat %>%
    filter(trait_trans == "Dry_Mass_g_log")
  fit <- lmer(mean ~ Gradient * Elevation_m + (1|Site), data = dd)

  newdat <- dd %>%
    distinct(Elevation_m, Gradient) %>%
    mutate(mean = predict(fit, newdat, re.form = NA))

  mm <- model.matrix(terms(fit), newdat)

  dry <- newdat %>%
    mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
           tvar1 = pvar1 + VarCorr(fit)$Site[1],  ## must be adapted for more complex models
           cmult = 1.96) %>%
    mutate(plo = mean - cmult*sqrt(pvar1),
           phi = mean + cmult*sqrt(pvar1),
           tlo = mean - cmult*sqrt(tvar1),
           thi = mean + cmult*sqrt(tvar1))

  gdry <- g0 %+% subset(dat, trait_trans == "Dry_Mass_g_log") +
    geom_line(data = dry) +
    geom_ribbon(data = dry, aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_fill_manual(values = c("green4", "grey"), labels = c("Birdcliff", "Reference")) +
    theme(axis.text.x = element_blank())


  #LDMC
  dd <- dat %>%
    filter(trait_trans == "LDMC")
  fit <- lmer(mean ~ Gradient * Elevation_m + (1|Site), data = dd)

  newdat <- dd %>%
    distinct(Elevation_m, Gradient) %>%
    mutate(mean = predict(fit, newdat, re.form = NA))

  mm <- model.matrix(terms(fit), newdat)

  ldmc <- newdat %>%
    mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
           tvar1 = pvar1 + VarCorr(fit)$Site[1],  ## must be adapted for more complex models
           cmult = 1.96) %>%
    mutate(plo = mean - cmult*sqrt(pvar1),
           phi = mean + cmult*sqrt(pvar1),
           tlo = mean - cmult*sqrt(tvar1),
           thi = mean + cmult*sqrt(tvar1))

  gldmc <- g0 %+% subset(dat, trait_trans == "LDMC") +
    geom_line(data = ldmc) +
    geom_ribbon(data = ldmc, aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_fill_manual(values = c("green4", "grey"), labels = c("Birdcliff", "Reference")) +
    annotate("text", x = 50, y = 0.12, label = "G+E", size = 3) +
    theme(axis.text.x = element_blank())


  #Leaf_Area_cm2_log
  dd <- dat %>%
    filter(trait_trans == "Leaf_Area_cm2_log")
  fit <- lmer(mean ~ Gradient * Elevation_m + (1|Site), data = dd)

  newdat <- dd %>%
    distinct(Elevation_m, Gradient) %>%
    mutate(mean = predict(fit, newdat, re.form = NA))

  mm <- model.matrix(terms(fit), newdat)

  area <- newdat %>%
    mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
           tvar1 = pvar1 + VarCorr(fit)$Site[1],  ## must be adapted for more complex models
           cmult = 1.96) %>%
    mutate(plo = mean - cmult*sqrt(pvar1),
           phi = mean + cmult*sqrt(pvar1),
           tlo = mean - cmult*sqrt(tvar1),
           thi = mean + cmult*sqrt(tvar1))

  garea <- g0 %+% subset(dat, trait_trans == "Leaf_Area_cm2_log") +
    geom_line(data = area) +
    geom_ribbon(data = area, aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_fill_manual(values = c("green4", "grey"), labels = c("Birdcliff", "Reference")) +
    annotate("text", x = 50, y = 0.5, label = "G+E", size = 3) +
    theme(axis.text.x = element_blank())


  #N_percent
  dd <- dat %>%
    filter(trait_trans == "N_percent")
  fit <- lmer(mean ~ Gradient * Elevation_m + (1|Site), data = dd)

  newdat <- dd %>%
    distinct(Elevation_m, Gradient) %>%
    mutate(mean = predict(fit, newdat, re.form = NA))

  mm <- model.matrix(terms(fit), newdat)

  n <- newdat %>%
    mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
           tvar1 = pvar1 + VarCorr(fit)$Site[1],  ## must be adapted for more complex models
           cmult = 1.96) %>%
    mutate(plo = mean - cmult*sqrt(pvar1),
           phi = mean + cmult*sqrt(pvar1),
           tlo = mean - cmult*sqrt(tvar1),
           thi = mean + cmult*sqrt(tvar1))

  gn <- g0 %+% subset(dat, trait_trans == "N_percent") +
    geom_line(data = n) +
    geom_ribbon(data = n, aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_fill_manual(values = c("green4", "grey"), labels = c("Birdcliff", "Reference")) +
    annotate("text", x = 50, y = 3.8, label = "G+E", size = 3) +
    theme(axis.text.x = element_blank())


  #NP_ratio
  dd <- dat %>%
    filter(trait_trans == "NP_ratio")
  fit <- lmer(mean ~ Gradient * Elevation_m + (1|Site), data = dd)

  newdat <- dd %>%
    distinct(Elevation_m, Gradient) %>%
    mutate(mean = predict(fit, newdat, re.form = NA))

  mm <- model.matrix(terms(fit), newdat)

  np <- newdat %>%
    mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
           tvar1 = pvar1 + VarCorr(fit)$Site[1],  ## must be adapted for more complex models
           cmult = 1.96) %>%
    mutate(plo = mean - cmult*sqrt(pvar1),
           phi = mean + cmult*sqrt(pvar1),
           tlo = mean - cmult*sqrt(tvar1),
           thi = mean + cmult*sqrt(tvar1))

  gnp <- g0 %+% subset(dat, trait_trans == "NP_ratio") +
    geom_line(data = np) +
    geom_ribbon(data = np, aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_fill_manual(values = c("green4", "grey"), labels = c("Birdcliff", "Reference"))


  #P_percent
  dd <- dat %>%
    filter(trait_trans == "P_percent")
  fit <- lmer(mean ~ Gradient * Elevation_m + (1|Site), data = dd)

  newdat <- dd %>%
    distinct(Elevation_m, Gradient) %>%
    mutate(mean = predict(fit, newdat, re.form = NA))

  mm <- model.matrix(terms(fit), newdat)

  p <- newdat %>%
    mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
           tvar1 = pvar1 + VarCorr(fit)$Site[1],  ## must be adapted for more complex models
           cmult = 1.96) %>%
    mutate(plo = mean - cmult*sqrt(pvar1),
           phi = mean + cmult*sqrt(pvar1),
           tlo = mean - cmult*sqrt(tvar1),
           thi = mean + cmult*sqrt(tvar1))

  gp <- g0 %+% subset(dat, trait_trans == "P_percent") +
    geom_line(data = p) +
    geom_ribbon(data = p, aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_fill_manual(values = c("green4", "grey"), labels = c("Birdcliff", "Reference")) +
    annotate("text", x = 50, y = 0.38, label = "GxE", size = 3)


  #Plant_Height_cm_log
  dd <- dat %>%
    filter(trait_trans == "Plant_Height_cm_log")
  fit <- lmer(mean ~ Gradient * Elevation_m + (1|Site), data = dd)

  newdat <- dd %>%
    distinct(Elevation_m, Gradient) %>%
    mutate(mean = predict(fit, newdat, re.form = NA))

  mm <- model.matrix(terms(fit), newdat)

  height <- newdat %>%
    mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
           tvar1 = pvar1 + VarCorr(fit)$Site[1],  ## must be adapted for more complex models
           cmult = 1.96) %>%
    mutate(plo = mean - cmult*sqrt(pvar1),
           phi = mean + cmult*sqrt(pvar1),
           tlo = mean - cmult*sqrt(tvar1),
           thi = mean + cmult*sqrt(tvar1))

  gheight <- g0 %+% subset(dat, trait_trans == "Plant_Height_cm_log") +
    geom_line(data = height) +
    geom_ribbon(data = height, aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_fill_manual(values = c("green4", "grey"), labels = c("Birdcliff", "Reference")) +
    annotate("text", x = 50, y = 1.8, label = "E", size = 3) +
    theme(axis.text.x = element_blank())


  #SLA_cm2_g
  dd <- dat %>%
    filter(trait_trans == "SLA_cm2_g")
  fit <- lmer(mean ~ Gradient * Elevation_m + (1|Site), data = dd)

  newdat <- dd %>%
    distinct(Elevation_m, Gradient) %>%
    mutate(mean = predict(fit, newdat, re.form = NA))

  mm <- model.matrix(terms(fit), newdat)

  sla <- newdat %>%
    mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
           tvar1 = pvar1 + VarCorr(fit)$Site[1],  ## must be adapted for more complex models
           cmult = 1.96) %>%
    mutate(plo = mean - cmult*sqrt(pvar1),
           phi = mean + cmult*sqrt(pvar1),
           tlo = mean - cmult*sqrt(tvar1),
           thi = mean + cmult*sqrt(tvar1))

  gsla <- g0 %+% subset(dat, trait_trans == "SLA_cm2_g") +
    geom_line(data = sla) +
    geom_ribbon(data = sla, aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_fill_manual(values = c("green4", "grey"), labels = c("Birdcliff", "Reference")) +
    labs(y = "Bootstrapped trait mean") +
    annotate("text", x = 50, y = 249, label = "GxE", size = 3) +
    theme(axis.text.x = element_blank(),
          axis.title.y = element_text(margin = margin(r = 20)))


  #Thickness_mm_log
  dd <- dat %>%
    filter(trait_trans == "Thickness_mm_log")
  fit <- lmer(mean ~ Gradient * Elevation_m + (1|Site), data = dd)

  newdat <- dd %>%
    distinct(Elevation_m, Gradient) %>%
    mutate(mean = predict(fit, newdat, re.form = NA))

  mm <- model.matrix(terms(fit), newdat)

  thick <- newdat %>%
    mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
           tvar1 = pvar1 + VarCorr(fit)$Site[1],  ## must be adapted for more complex models
           cmult = 1.96) %>%
    mutate(plo = mean - cmult*sqrt(pvar1),
           phi = mean + cmult*sqrt(pvar1),
           tlo = mean - cmult*sqrt(tvar1),
           thi = mean + cmult*sqrt(tvar1))

  gthick <- g0 %+% subset(dat, trait_trans == "Thickness_mm_log") +
    geom_line(data = thick) +
    geom_ribbon(data = thick, aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_fill_manual(values = c("green4", "grey"), labels = c("Birdcliff", "Reference")) +
    theme(axis.text.x = element_blank())



  # patchwork
  p2 <- ggplot() +
    annotate("text", x = 1, y = 2, label = "Elevation m a.s.l.") +
    theme_void()

  library(cowplot)
  legend <- gheight + theme(legend.position = "bottom")
  legend <- cowplot::get_legend(legend)


  layout <- "
  ABCD
  ABCD
  EFGH
  EFGH
  JIKL
  JIKL
  MMMM
  NNNN
"

  trait_plot <- wrap_plots(gheight, gdry, garea, gthick, gsla, gldmc, gc, gn, gcn, gp, gnp, gdn) + p2 + legend  + plot_layout(design = layout)

  return(trait_plot)

}


#
# library(MuMIn)
# dd <- dat %>%
#   filter(trait_trans == "dC13_permil")
# fit1 <- lmer(mean ~ Gradient * Elevation_m + (1|Site), data = dd, na.action = "na.fail", REML = FALSE)
# model.set <- dredge(fit1, rank = "AICc", extra = "R^2")
# model.set
#
#
# # Likelihoot ratio test result
# result <- dat %>%
#   distinct(trait_trans) %>%
#   mutate(model = c("", "G", "???", "GxE", "", "G+E", "G+E", "G", "", "G*E", "E", "G*E", ""))
#
#
# dat %>%
#   group_by(trait_trans) %>%
#   nest(data = -c(trait_trans)) %>%
#   mutate(mod = map(data, ~lmer(mean ~ Gradient * Elevation_m + (1|Site), data = .x)))

