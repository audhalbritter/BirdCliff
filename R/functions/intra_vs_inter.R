#### Intraspecific vs. interspecific trait variation ####

# some fancy old school functions to calculate intra vs. interspecific trait variation

#Intra_vs_Inter <- function(traits_raw, trait_mean){
Intra_vs_Inter <- function(trait_mean){

  # #prepare trait data
  # trait <- traits_raw %>%
  #   select(Gradient, Site, PlotID, Taxon, trait_trans, value_trans)
  #
  # var_res <- data.frame()
  #
  # for(i in unique(trait$trait_trans)){
  #   v <- varcomp(lme(value_trans ~ 1, random = ~ 1|Taxon, data = trait %>% filter(trait_trans == i), na.action = na.omit), 1)[c(1,2)]
  #
  #   v$trait <- i
  #
  #   v <- unlist(v)
  #
  #   var_res <- bind_rows(var_res, v)
  #
  # }

  var_split <- trait_mean %>%
    group_by(trait_trans) %>%
    do(test = trait.flex.anova(~Site * Gradient, mean, mean_noitv, data = .))

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
    as_tibble() %>%
    mutate(level = trimws(level)) %>%
    filter(RelSumSq.Turnover < 999) %>%
    rename(Turnover = RelSumSq.Turnover, Intraspecific = RelSumSq.Intraspec., Covariation = RelSumSq.Covariation, Total = RelSumSq.Total) %>%
    pivot_longer(cols = c(-trait, -level), names_to = "variable", values_to = "value") %>%
    filter(variable == "Total") %>%
    filter(level != "Total") %>%
    #mutate(level = factor(level, levels = c("Habitat", "Treatment", "Habitat:Treatment", "Residuals"))) %>%
    #mutate(level = plyr::mapvalues(level, from = c("Habitat", "Treatment", "Habitat:Treatment", "Residuals"), to = c("H", "T", "HxT", "Resid"))) %>%
    #mutate(trait = plyr::mapvalues(trait, from = c("SLA_cm2_g", "LDMC", "Leaf_Area_cm2", "Leaf_Thickness_mm", "N_percent", "C_percent", "P_Ave", "CN_ratio", "dC13_percent", "dN15_percent", "Dry_Mass_g", "Plant_Height_cm"), to = c("`SLA`*` `*(cm^2/g)", "`LDMC`*` `*(g/g)", "'Leaf'*' '*'Area'*' '*(cm^2)", "'Leaf'*' '*'Thickness'*' '*(mm)", "'N'*' '*'(%)'", "'C'*' '*'(%)'", "'P'*' '*'(%)'", "'C'*':'*'N'", "paste(delta^13, 'C'*' '*'(\u2030)')", "paste(delta^15, 'N'*' '*'(\u2030)')", "'Dry'*' '*'Mass'*' '*'(g)'", "'Plant'*' '*'Height'*' '*'(cm)'"))) %>%
    #mutate(trait = factor(trait, levels = c("'Plant'*' '*'Height'*' '*'(cm)'", "`SLA`*` `*(cm^2/g)", "'Dry'*' '*'Mass'*' '*'(g)'","'Leaf'*' '*'Area'*' '*(cm^2)",  "'Leaf'*' '*'Thickness'*' '*(mm)", "`LDMC`*` `*(g/g)",  "'N'*' '*'(%)'", "'C'*' '*'(%)'", "'P'*' '*'(%)'", "'C'*':'*'N'", "paste(delta^13, 'C'*' '*'(\u2030)')", "paste(delta^15, 'N'*' '*'(\u2030)')" )))

  return(var_split)
}
