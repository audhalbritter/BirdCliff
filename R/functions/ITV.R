#### Intraspecific vs. interspecific variation ####


Intra_vs_Inter <- function(traits_raw, trait_mean){

  #prepare trait data
  trait <- traits_raw %>%
    select(Gradient, Site, PlotID, Individual_nr, Taxon, trait_trans, value_trans) %>%
    mutate(GS = paste0(Gradient, Site))

  var_res <- data.frame()

  for(i in unique(trait$trait_trans)){
    v <- varcomp(lme(value_trans~1, random=~1|Taxon, data=trait %>% filter(trait_trans == i), na.action = na.omit), 1)[c(1,2)]

    v$trait <- i

    v <- unlist(v)

    var_res <- bind_rows(var_res, v)

  }

  var_split <- trait_mean %>%
    group_by(trait_trans) %>%
    do(test = trait.flex.anova(~Gradient*Elevation_m, mean, mean_noitv, data = .))

  var_split_exp <- data.frame(RelSumSq.Turnover = 1000, RelSumSq.Intraspec. = 1000, RelSumSq.Covariation = 1000, RelSumSq.Total = 1000, trait = "E", level = "F")

  for(i in 1:nrow(var_split)){
    out <- as.data.frame(var_split$test[[i]][2])
    out$trait <- as.factor(rep(var_split[i,1], 5))
    out$level <- rownames(out)
    var_split_exp <- rbind(var_split_exp, out)
  }

  return(var_split_exp)

}


Intra_vs_Inter_var_split <- function(var_split_exp){

  var_split <- var_split_exp %>%
    mutate(level = trimws(level)) %>%
    filter(RelSumSq.Turnover < 999) %>%
    rename(Turnover = RelSumSq.Turnover, Intraspecific = RelSumSq.Intraspec., Covariation = RelSumSq.Covariation, Total = RelSumSq.Total) %>%
    gather(key = variable, value = value, -trait, -level) %>%
    filter(variable != "Covariation", level != "Total", variable != "Total") %>%
    filter(variable != "Total") %>%
    filter(level != "Total") %>%
    mutate(level = factor(level, levels = c("Gradient", "Elevation_m", "Gradient:Elevation_m", "Residuals"))) %>%
    mutate(level = plyr::mapvalues(level, from = c("Gradient", "Elevation_m", "Gradient:Elevation_m", "Residuals"), to = c("G", "E", "GxE", "Resid"))) #%>%
    # mutate(trait = plyr::mapvalues(trait, from = c("SLA_cm2_g", "LDMC", "Leaf_Area_cm2", "Leaf_Thickness_mm", "N_percent", "C_percent", "P_Ave", "CN_ratio", "dC13_percent", "dN15_percent", "Dry_Mass_g", "Plant_Height_cm"), to = c("`SLA`*` `*(cm^2/g)", "`LDMC`*` `*(g/g)", "'Leaf'*' '*'Area'*' '*(cm^2)", "'Leaf'*' '*'Thickness'*' '*(mm)", "'N'*' '*'(%)'", "'C'*' '*'(%)'", "'P'*' '*'(%)'", "'C'*':'*'N'", "paste(delta^13, 'C'*' '*'(\u2030)')", "paste(delta^15, 'N'*' '*'(\u2030)')", "'Dry'*' '*'Mass'*' '*'(g)'", "'Plant'*' '*'Height'*' '*'(cm)'"))) %>%
    # mutate(trait = factor(trait, levels = c("'Plant'*' '*'Height'*' '*'(cm)'", "`SLA`*` `*(cm^2/g)", "'Dry'*' '*'Mass'*' '*'(g)'","'Leaf'*' '*'Area'*' '*(cm^2)",  "'Leaf'*' '*'Thickness'*' '*(mm)", "`LDMC`*` `*(g/g)",  "'N'*' '*'(%)'", "'C'*' '*'(%)'", "'P'*' '*'(%)'", "'C'*':'*'N'", "paste(delta^13, 'C'*' '*'(\u2030)')", "paste(delta^15, 'N'*' '*'(\u2030)')" )))

  return(var_split)
}




#### ITV plot - importance intra specific trait variation ####
# fixed (species turnover): inter specific trait variation
# specific: inter + intra specific trait variation
# Difference (intra specific)  = fixed - specific

make_itv_figure <- function(trait_mean){

  traitMean <- fancy_trait_name_dictionary(trait_mean) %>%
    mutate(itv_diff = mean-mean_noitv) %>%
    mutate(trait_fancy = factor(trait_fancy, levels = c("Height cm", "Dry mass g", "Area cm2", "Thickness mm", "SLA cm2/g", "LDMC", "C %", "N %", "CN", "P %", "NP", "dC13 ‰", "dN15 ‰")))

  # t test if intra is important (different from zero)
  t_test_itv <- traitMean %>%
    group_by(trait_fancy, Gradient) %>%
    summarise(P = t.test(itv_diff, mu = 0)$p.value,
              Sig = ifelse(P < 0.05, "*", ifelse(P < 0.1 & P > 0.05, "+", "")),
              MaxWidth = max(itv_diff)) %>%
    ungroup()

  itv_plot <- traitMean %>%
    ggplot() +
    geom_boxplot(aes(x = Gradient, y = itv_diff)) +
    geom_hline(aes(yintercept = 0)) +
    geom_blank(aes(x = Gradient, y = itv_diff + itv_diff*0.6)) +
    scale_fill_manual(values = c("darkgray", "red")) +
    theme_classic() +
    theme(text = element_text(size = 15),
          strip.background = element_blank(),
          legend.position = "top") +
    facet_wrap(~trait_trans, scales = "free_y", labeller = label_parsed) +
    #geom_text(aes(label = Sig, y = Inf, x = Gradient, vjust = 1), position = position_dodge(0.75), data = t_test_itv, size = 4.5) +
    labs(y = "Mean Trait Value - Mean Trait Value (no ITV)")

  return(itv_plot)

}




make_intra_vs_inter_figure <- function(trait_mean, var_split_exp){

  names <- fancy_trait_name_dictionary(trait_mean) %>%
    distinct(trait_fancy) %>%
    mutate(trait = factor(1:13))


  varpart_graph <- var_split_exp %>%
    mutate(level = trimws(level)) %>%
    filter(RelSumSq.Turnover < 999) %>%
    rename(Turnover = RelSumSq.Turnover, Intraspecific = RelSumSq.Intraspec., Covariation = RelSumSq.Covariation, Total = RelSumSq.Total) %>%
    gather(key = variable, value = value, -trait, -level) %>%
    filter(variable != "Covariation", level != "Total", variable != "Total") %>%
    mutate(level = factor(level, levels = c("Gradient", "Elevation_m", "Gradient:Elevation_m", "Residuals"))) %>%
    mutate(level = plyr::mapvalues(level, from = c("Gradient", "Elevation_m", "Gradient:Elevation_m", "Residuals"), to = c("G", "E", "GxE", "Resid"))) %>%
    left_join(names, by = "trait") %>%
    ggplot() +
    geom_bar(aes(x = level, y = value, fill = variable), stat = "identity") +
    geom_point(aes(x = level, y  = value), data = variation_split, size = 1) +
    facet_wrap(~ trait_fancy) +
    #facet_wrap(~trait, nrow = 3, labeller = label_parsed) +
    theme_minimal() +
    theme(text = element_text(size = 15), legend.position = "top",
          strip.background = element_blank()) +
    xlab(NULL) +
    ylab("Proportion Variation Explained") +
    scale_fill_manual(values = c("blue", "darkorange"), name = "Source of Variation") +
    scale_x_discrete(drop = FALSE)

}


# mutate(trait = plyr::mapvalues(trait, from = c("SLA_cm2_g", "LDMC", "Leaf_Area_cm2", "Leaf_Thickness_mm", "N_percent", "C_percent", "P_Ave", "CN_ratio", "dC13_percent", "dN15_percent", "Dry_Mass_g", "Plant_Height_cm"), to = c("`SLA`*` `*(cm^2/g)", "`LDMC`*` `*(g/g)", "'Leaf'*' '*'Area'*' '*(cm^2)", "'Leaf'*' '*'Thickness'*' '*(mm)", "'N'*' '*'(%)'", "'C'*' '*'(%)'", "'P'*' '*'(%)'", "'C'*':'*'N'", "paste(delta^13, 'C'*' '*'(\u2030)')", "paste(delta^15, 'N'*' '*'(\u2030)')", "'Dry'*' '*'Mass'*' '*'(g)'", "'Plant'*' '*'Height'*' '*'(cm)'"))) %>%
# mutate(trait = factor(trait, levels = c("'Plant'*' '*'Height'*' '*'(cm)'", "'Dry'*' '*'Mass'*' '*'(g)'","'Leaf'*' '*'Area'*' '*(cm^2)", "`SLA`*` `*(cm^2/g)",  "'Leaf'*' '*'Thickness'*' '*(mm)", "`LDMC`*` `*(g/g)",  "'N'*' '*'(%)'", "'C'*' '*'(%)'", "'P'*' '*'(%)'", "'C'*':'*'N'", "paste(delta^13, 'C'*' '*'(\u2030)')", "paste(delta^15, 'N'*' '*'(\u2030)')" ))) %>%
