si_figures_plan <- list(

  # climate data
  tar_target(
    name = climate_analysis,
    command = {

      dat <- climate_data %>%
      mutate(Variable = recode(Variable, "SoilMoisture" = "soil moisture in %", "SoilTemperature" = "soil temperature in °C"),
             GS = paste0(Gradient, Site)) %>%
      left_join(coordinates, by = c("Gradient", "Site", "PlotID"))

      # model selection
      dat %>%
        group_by(Variable) %>%
        nest() %>%
        mutate(model.set = map(data, ~{
          mod <- lmer(Value ~  Gradient * Elevation_m + (1|GS), REML = FALSE, na.action = "na.fail", data = .x)
          model.set = dredge(mod, rank = "AICc", extra = "R^2")
        })) %>%
        unnest(model.set)

      # moisture: GxE, temperature G
      fit_m <- lmer(Value ~  Gradient * Elevation_m + (1|GS), data = dat %>% filter(Variable == "soil moisture in %"))

      fit_t <- lmer(Value ~  Gradient + (1|GS), data = dat %>% filter(Variable == "soil temperature in °C"))

      bind_rows(
        moisture = tidy(fit_m),
        temperature = tidy(fit_t),
        .id = "Variable") %>%
        left_join(tibble(Rm = r.squaredGLMM(fit_m) %>% as.vector(),
                               Rc = r.squaredGLMM(fit_t) %>% as.vector(),
                               Variable = c("moisture", "temperature")),
                  by = "Variable")
    }
  ),

  tar_target(
    name = climate_plot,
    command = {

      # climate data
      dat <- climate_data %>%
        mutate(Variable = recode(Variable, "SoilMoisture" = "soil moisture in %", "SoilTemperature" = "soil temperature in °C"),
               GS = paste0(Gradient, Site)) %>%
        left_join(coordinates, by = c("Gradient", "Site", "PlotID"))

      # soil moisture
      dd <- dat %>%
        filter(Variable == "soil moisture in %")
      fit <- lmer(Value ~  Gradient * Elevation_m + (1|GS), data = dd)

      newdat <- dd %>%
        distinct(Elevation_m, Gradient) %>%
        mutate(Value = 0)
      newdat$Value <-  predict(fit, newdat, re.form = NA)

      mm <- model.matrix(terms(fit), newdat)

      soil_m <- newdat %>%
        mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
               tvar1 = pvar1 + VarCorr(fit)$GS[1],  ## must be adapted for more complex models
               cmult = 1.96) %>%
        mutate(plo = Value - cmult*sqrt(pvar1),
               phi = Value + cmult*sqrt(pvar1),
               tlo = Value - cmult*sqrt(tvar1),
               thi = Value + cmult*sqrt(tvar1))

      g0 <- dat |>
        filter(Variable == "soil moisture in %") |>
        ggplot(aes(x = Elevation_m, y = Value, colour = Gradient, fill = Gradient)) +
        geom_point(alpha = 0.5) +
        scale_colour_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
        labs(x = "Elevation in m a.s.l.", y = "") +
        theme_minimal() +
        theme(legend.position = "bottom")

        soil_m_plot <- g0 +
        geom_line(data = soil_m) +
        geom_ribbon(data = soil_m, aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
        scale_fill_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
        annotate("text", x = Inf, y = Inf, label = "GxE", size = 3, hjust = 1, vjust = 1) +
        labs(y = "Soil moisture in %")



      # soil temperature
      dd <- dat %>%
        filter(Variable == "soil temperature in °C")
      fit <- lmer(Value ~  Gradient + (1|GS), data = dd)

      newdat <- dd %>%
        distinct(Elevation_m, Gradient) %>%
        mutate(Value = 0)
      newdat$Value <-  predict(fit, newdat, re.form = NA)

      mm <- model.matrix(terms(fit), newdat)

      soil_t <- newdat %>%
        mutate(pvar1 = diag(mm %*% tcrossprod(vcov(fit), mm)),
               tvar1 = pvar1 + VarCorr(fit)$GS[1],  ## must be adapted for more complex models
               cmult = 1.96) %>%
        mutate(plo = Value - cmult*sqrt(pvar1),
               phi = Value + cmult*sqrt(pvar1),
               tlo = Value - cmult*sqrt(tvar1),
               thi = Value + cmult*sqrt(tvar1))

      soil_t_plot <- g0 %+% subset(dat, Variable == "soil temperature in °C") +
        geom_line(data = soil_t) +
        geom_ribbon(data = soil_t, aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
        scale_fill_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
        annotate("text", x = Inf, y = Inf, label = "G", size = 3, hjust = 1, vjust = 1) +
        labs(y = "Soil temperature in °C")

      climate_plot <- soil_m_plot + soil_t_plot + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

    }),


  # species ordination
  tar_target(
    name = ordination_plot,
    command = make_ordination_plot(comm_raw,
                                   NMDS = sp_ordination[[1]],
                                   fNMDS = sp_ordination[[2]])),


  # check the nr of dimensions for NMDS
  tar_target(
    name = stress_plot,
    command = check_dimensions_NMDS(comm_raw)
  ),


  # check the nr of dimensions for NMDS
  tar_target(
    name = imputation_plot,
    command = {

      gradient <- c(
        B = "Bird cliff",
        C = "Reference")

      trait_names <- c(
        "Plant_Height_cm_log" = "Height cm",
        "Dry_Mass_g_log" = "Dry mass g",
        "Thickness_mm_log" = "Thickness mm",
        "Leaf_Area_cm2_log" = "Area cm2",
        "SLA_cm2_g" = "SLA cm2/g",
        "LDMC" = "LDMC",
        "C_percent" = "C %",
        "N_percent" = "N % ",
        "CN_ratio" = "CN",
        "dN15_permil" = "δN15 ‰",
        "dC13_permil" = "δC13 ‰",
        "P_percent" = "P %",
        "NP_ratio" = "NP")

      #check trait coverage
      imputation_plot <- trait_impute %>%
        autoplot(.) +
        scale_y_continuous(breaks = c(0, 0.5, 1)) +
        scale_x_discrete(breaks = c("B_1_D", "B_2_D", "B_3_D", "B_4_D", "B_5_D", "C_1_D", "C_2_D", "C_3_D", "C_4_D", "C_5_D", "C_6_D", "C_7_D"),
                         labels = c("1", "2", "3", "4", "5", "1", "2", "3", "4", "5", "6", "7")) +
        facet_grid(trait_trans ~ Gradient, scales = "free_x",
                   labeller = labeller(Gradient = gradient, trait_trans = trait_names)) +
        labs(x = "Site") +
        theme_minimal() +
        theme(strip.text.y = element_text(angle = 0))

      return(imputation_plot)

    }
  )



  # correlation plot
  # tar_target(
  #   name = correlation_plot,
  #   command = {
  #     # trait correlation both gradients
  #     trait_corr_all = get_trait_correlations(trait_mean) %>%
  #       mutate(r = if_else(r == 1, NA_real_, r)) %>%
  #       ggplot(aes(x = trait1, y = trait2, fill = r, label = round(r_sig, 2))) +
  #       geom_tile() +
  #       labs(x = NULL, y = NULL,
  #            title = "Both gradients",
  #            fill = "Pearson's\nCorrelation") +
  #       scale_fill_gradient2(mid = "#f7f7f7", low = "#f1a340",
  #                            high = "#998ec3",
  #                            na.value = "grey80",
  #                            limits = c(-1, 1)) +
  #       geom_text() +
  #       theme_minimal() +
  #       theme(legend.position = "top") +
  #       scale_x_discrete(expand=c(0,0)) +
  #       scale_y_discrete(expand=c(0,0))
  #
  #     # Bird cliff
  #     trait_corr_B = get_trait_correlations(trait_mean %>%
  #                                             filter(Gradient == "B")) %>%
  #       mutate(r = if_else(r == 1, NA_real_, r)) %>%
  #       ggplot(aes(x = trait1, y = trait2, fill = r, label = round(r_sig, 2))) +
  #       geom_tile() +
  #       labs(x = NULL, y = NULL,
  #            title = "Bird cliff",
  #            fill = "Pearson's\nCorrelation") +
  #       scale_fill_gradient2(mid = "#f7f7f7", low = "#f1a340",
  #                            high = "#998ec3",
  #                            na.value = "grey80",
  #                            limits = c(-1, 1)) +
  #       geom_text() +
  #       theme_minimal() +
  #       theme(legend.position = "none") +
  #       scale_x_discrete(expand=c(0,0)) +
  #       scale_y_discrete(expand=c(0,0))
  #
  #     # Control
  #     trait_corr_C = get_trait_correlations(trait_mean %>%
  #                                             filter(Gradient == "C")) %>%
  #       mutate(r = if_else(r == 1, NA_real_, r)) %>%
  #       ggplot(aes(x = trait1, y = trait2, fill = r, label = round(r_sig, 2))) +
  #       geom_tile() +
  #       labs(x = NULL, y = NULL,
  #            title = "Reference",
  #            fill = "Pearson's\nCorrelation") +
  #       scale_fill_gradient2(mid = "#f7f7f7", low = "#f1a340",
  #                            high = "#998ec3",
  #                            na.value = "grey80",
  #                            limits = c(-1, 1)) +
  #       geom_text() +
  #       theme_minimal() +
  #       theme(legend.position = "none") +
  #       scale_x_discrete(expand=c(0,0)) +
  #       scale_y_discrete(expand=c(0,0))
  #
  #     trait_corr = trait_corr_all / trait_corr_B / trait_corr_C
  #   })


)
