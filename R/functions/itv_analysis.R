### ITV ANALYSIS ###

# total <- specif.avg <- mean (local mean)
# turnover <- const.avg <- mean_noitv (global mean)
# intraspecific <- specific.avg - const.avg

#tar_load(trait_mean)

make_ITV_analysis <- function(trait_mean){

  trait_long <- trait_mean %>%
    # calculate ITV = specific - constant mean
    mutate(diff = mean - mean_noitv) %>%
    # make long table
    pivot_longer(cols = c(mean, mean_noitv, diff), names_to = "mean", values_to = "value")

  itv_output <- trait_long %>%
    group_by(Gradient, trait_trans, mean) %>%
    nest() %>%
    # anova for each mean (3x)
    mutate(estimate = map(data, ~{
      mod <- aov(value ~ Site, data =  .x)
      # output tidy results
      estimates = tidy(mod)
    })) %>%
    unnest(estimate)

  return(itv_output)

}



make_ITV_plot <- function(itv_output){

  variance_part <- itv_output %>%
    # select important columns: sumsq = SS
    select(Gradient, trait_trans, mean, term, sumsq) %>%
    # make wide table
    pivot_wider(names_from = mean, values_from = sumsq) %>%
    # rename columns
    rename("total" = mean, "turnover" = mean_noitv, "intraspecific" = diff) %>%
    # calculate total SS for total variation
    group_by(Gradient, trait_trans) %>%
    mutate(total_var = sum(total)) %>%
    ungroup() %>%
    # calculate covariation
    mutate(covariation = total - turnover - intraspecific,
           # calculate proportion explained variation
           total_p = total/total_var,
           turnover_p = turnover/total_var,
           intra_p = intraspecific/total_var,
           covariation_p = covariation/total_var) %>%
    # make long table
    pivot_longer(cols = c(total_p:covariation_p), names_to = "process", values_to = "value")


  ITV_plot <- fancy_trait_name_dictionary(variance_part) %>%
    # filter only turnover and ITV
    filter(process %in% c("turnover_p", "intra_p"),
           # remove residuals
           term == "Site") %>%
    mutate(process = recode(process, turnover_p = "turnover", intra_p = "ITV")) %>%
    ggplot(aes(x = Gradient, y = value, fill = process)) +
    geom_col() +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "inferno") +
    labs(x = "", y = "Variation explained") +
    facet_wrap(~trait_fancy, scales = "free_y") +
    theme_minimal()

  return(ITV_plot)

}








# test code with one trait
# m <- trait_mean %>%
#   filter(trait_trans %in% c("LDMC"),
#          Gradient == "B") %>%
#   # specific - constant
#   mutate(diff = mean - mean_noitv)
#
# # With effect of Site
# # specific
# fit_specific <- aov(mean ~ Site, data = m)
#
# # constant
# fit_constant <- aov(mean_noitv ~ Site, data = m)
#
# # diff
# fit_diff <- aov(diff ~ Site, data = m)
#
# out <- bind_rows(
#   total = tidy(fit_specific),
#   turnover = tidy(fit_constant),
#   intraspecific = tidy(fit_diff),
#   .id = "model"
# )
#
# out %>%
#   select(model, term, sumsq) %>%
#   pivot_wider(names_from = model, values_from = sumsq) %>%
#   mutate(covariation = total - turnover - intraspecific) %>%
#   pivot_longer(cols = c(total:covariation), names_to = "process", values_to = "value") %>%
#   filter(process %in% c("turnover", "intraspecific")) %>%
#   ggplot(aes(x = term, y = value, fill = process)) +
#   geom_col() +
#   scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "inferno") +
#   labs(x = "", y = "Variation explained")






#centre elevation to fix singular fix

# trait_mean %>%
#   mutate(Elevation_cent = scale(Elevation_m, center = TRUE, scale = FALSE)[1] %>% as.vector())
#
#
# scale(trait_mean$Elevation_m, center = TRUE, scale = FALSE)[,1] %>% as_tibble()



# make screeplot self


# PCAs do rda test
# could also use PCA for species composition, square root/tranform species data
