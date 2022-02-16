figures_plan <- list(

  # FIGURE 2a: trait change along gradients
  tar_target(
    name = trait_plot,
    command = {
      fancy_trait_name_dictionary(trait_mean) %>%
        ggplot(aes(x = Elevation_m, y = mean, colour = Gradient)) +
        geom_point(alpha = 0.3) +
        geom_smooth(method = "lm") +
        scale_colour_manual(values = c("green4", "grey"), labels = c("Birdcliff", "Reference")) +
        labs(x = "Elevation in m a.s.l.", y = "Bootstrapped trait mean") +
        facet_wrap(~ trait_fancy, scales = "free_y") +
        theme_minimal() +
        theme(legend.position = c(0.9, 0.05))

    }),

  # FIGURE 2B: trait ordination
  tar_target(
    name = trait_ordination_plot,
    command = {

      plot <- trait_pca[[1]] %>%
        ggplot(aes(x = PC1, y = PC2, colour = Site, shape = Gradient)) +
        geom_point(size = 3) +
        coord_equal() +
        scale_colour_manual(values = c("grey50", "pink", "lightblue", "red", "blue", "orange")) +
        scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1) +
        scale_shape_manual(values = c(16, 1), , labels = c("Birdcliff", "Control")) +
        labs(x = "PC 1", y = "PC 2") +
        theme_minimal()

  arrow <- trait_pca[[1]] %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_segment(data = trait_pca[[2]],
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "cm")),
                 colour = "grey50",
                 inherit.aes = FALSE) +
    geom_text(data = trait_pca[[2]],
              aes(x = PC1 * 1.1,y = PC2 * 1.1, label = Label),
              size = 3,
              inherit.aes = FALSE, colour = "black") +
    labs(x = "", y = "") +
    scale_x_continuous(expand = c(.2, 0)) +
    theme_minimal()

  plot + arrow

    }),

  # bryo plot
  tar_target(
    name = bryophyte_plot,
    command = {

      bryo_traits_raw %>%
        ggplot(aes(x = Elevation_m, y = Value, colour = Gradient)) +
        geom_point() +
        geom_smooth(method = "lm") +
        scale_colour_manual(values = c("green4", "grey"), labels = c("Birdcliff", "Control")) +
        labs(x = "Elevation in m a.s.l.", y = "Trait value") +
        facet_wrap(~ Trait, scales = "free_y")
    })

)



#### plot 4: ITV plot ####

make_itv_figure <- function(trait_mean){

  trait_mean <- trait_mean %>%
    mutate(itv_diff = mean-mean_noitv)

  t_test_itv <- trait_mean %>%
    group_by(Gradient, Site, trait_trans) %>%
    summarise(P = t.test(itv_diff, mu = 0)$p.value,
              Sig = ifelse(P < 0.05, "*", ifelse(P<0.1 & P > 0.05, "+", "")),
              MaxWidth = max(itv_diff)) %>%
    ungroup()

trait_mean %>%
    ggplot() +
    geom_boxplot(aes(x = Site, y = itv_diff, fill = Gradient)) +
    geom_hline(aes(yintercept = 0)) +
    geom_blank(aes(x = Site, y = itv_diff + itv_diff*0.6)) +
    scale_fill_manual(values = c("darkgray", "red")) +
    theme_classic() +
    theme(text = element_text(size = 15),
          strip.background = element_blank(),
          legend.position = "top") +
    facet_wrap(~trait_trans, scales = "free_y") +
    geom_text(aes(label = Sig, y = Inf, x = Site, group = Gradient, vjust = 1), position = position_dodge(0.75), data = t_test_itv, size = 4.5) +
    xlab("Habitat Type") +
    ylab("Mean Trait Value - Mean Trait Value (no ITV)")

  return(itv_plot)

}



make_intra_vs_inter_figure <- function(var_split_exp, var_split){

  varpart_graph <- var_split_exp %>%
    mutate(level = trimws(level)) %>%
    filter(RelSumSq.Turnover < 999) %>%
    rename(Turnover = RelSumSq.Turnover, Intraspecific = RelSumSq.Intraspec., Covariation = RelSumSq.Covariation, Total = RelSumSq.Total) %>%
    mutate(level = plyr::mapvalues(level, from = c("Site", "Site:Treatment"), to = c("Habitat", "Habitat:Treatment"))) %>%
    gather(key = variable, value = value, -trait, -level) %>%
    filter(variable != "Covariation", level != "Total", variable != "Total") %>%
    mutate(level = factor(level, levels = c("Habitat", "Treatment", "Habitat:Treatment", "Residuals"))) %>%
    mutate(level = plyr::mapvalues(level, from = c("Habitat", "Treatment", "Habitat:Treatment", "Residuals"), to = c("H", "T", "HxT", "Resid"))) %>%
    mutate(trait = plyr::mapvalues(trait, from = c("SLA_cm2_g", "LDMC", "Leaf_Area_cm2", "Leaf_Thickness_mm", "N_percent", "C_percent", "P_Ave", "CN_ratio", "dC13_percent", "dN15_percent", "Dry_Mass_g", "Plant_Height_cm"), to = c("`SLA`*` `*(cm^2/g)", "`LDMC`*` `*(g/g)", "'Leaf'*' '*'Area'*' '*(cm^2)", "'Leaf'*' '*'Thickness'*' '*(mm)", "'N'*' '*'(%)'", "'C'*' '*'(%)'", "'P'*' '*'(%)'", "'C'*':'*'N'", "paste(delta^13, 'C'*' '*'(\u2030)')", "paste(delta^15, 'N'*' '*'(\u2030)')", "'Dry'*' '*'Mass'*' '*'(g)'", "'Plant'*' '*'Height'*' '*'(cm)'"))) %>%
    mutate(trait = factor(trait, levels = c("'Plant'*' '*'Height'*' '*'(cm)'", "`SLA`*` `*(cm^2/g)", "'Dry'*' '*'Mass'*' '*'(g)'","'Leaf'*' '*'Area'*' '*(cm^2)",  "'Leaf'*' '*'Thickness'*' '*(mm)", "`LDMC`*` `*(g/g)",  "'N'*' '*'(%)'", "'C'*' '*'(%)'", "'P'*' '*'(%)'", "'C'*':'*'N'", "paste(delta^13, 'C'*' '*'(\u2030)')", "paste(delta^15, 'N'*' '*'(\u2030)')" ))) %>%
    ggplot() +
    geom_bar(aes(x = level, y = value, fill = variable), stat = "identity") +
    geom_point(aes(x = level, y  = value), data = var_split, size = 1) +
    facet_wrap(~trait, nrow = 3, labeller = label_parsed) +
    theme_classic() +
    theme(text = element_text(size = 15), legend.position = "top",
          strip.background = element_blank()) +
    xlab(NULL) +
    ylab("Proportion Variation Explained") +
    scale_fill_manual(values = c("blue", "darkorange"), name = "Source of Variation") +
    scale_x_discrete(drop = FALSE)

}
