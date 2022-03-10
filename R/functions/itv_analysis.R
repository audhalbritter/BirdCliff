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
    group_by(trait_trans, mean) %>%
    nest() %>%
    # anova for each mean (3x)
    mutate(estimate = map(data, ~{
      #mod <- aov(value ~ Site, data =  .x)
      mod <- aov(value ~ Gradient, data =  .x)
      # output tidy results
      estimates = tidy(mod)
    })) %>%
    unnest(estimate)

  return(itv_output)

}


make_ITV_plot <- function(itv_output){

  variance_part <- itv_output %>%
    # select important columns: sumsq = SS
    select(trait_trans, mean, term, sumsq) %>%
    # make wide table
    pivot_wider(names_from = mean, values_from = sumsq) %>%
    # rename columns
    rename("total" = mean, "turnover" = mean_noitv, "intraspecific" = diff) %>%
    # calculate total SS for total variation
    group_by(trait_trans) %>%
    mutate(total_var = sum(total),
           total_turnover = sum(turnover),
           total_ITV = sum(intraspecific)) %>%
    ungroup() %>%
    # calculate covariation
    mutate(covariation = total - turnover - intraspecific,
           # calculate proportion explained variation
           #total_p = total/total_var,
           turnover_p = turnover/total_var,
           intra_p = intraspecific/total_var,
           #covariation_p = covariation/total_var,
           total_turnover_p = total_turnover/total_var,
           total_ITV_p = total_ITV/total_var) %>%
    # make long table
    pivot_longer(cols = c(turnover_p:total_ITV_p), names_to = "process", values_to = "value") |>
    mutate(term = if_else(str_detect(process, "total_turnover|total_ITV"), "Total", term),
           process = recode(process, total_turnover_p = "turnover_p", total_ITV_p = "intra_p"))


  ITV_plot <- fancy_trait_name_dictionary(variance_part) %>%
    # filter only turnover and ITV
    #filter(process %in% c("turnover_p", "intra_p")) %>%
    mutate(process = recode(process, turnover_p = "turnover", intra_p = "ITV"),
           term = factor(term, levels = c("Total", "Gradient", "Residuals"))) %>%
    ggplot(aes(x = trait_fancy, y = value, fill = process)) +
    geom_col() +
    coord_flip() +
    scale_fill_viridis_d(begin = 0.25, end = 1, option = "viridis") +
    labs(x = "", y = "% explained") +
    facet_wrap(~ term, scales = "free_x") +
    theme_minimal()

  return(ITV_plot)

}





# itv_output <- trait_long %>%
#   group_by(trait_trans, mean) %>%
#   nest() %>%
#   # anova for each mean (3x)
#   mutate(estimate = map(data, ~{
#     mod <- aov(value ~ 1, data =  .x)
#     # output tidy results
#     estimates = tidy(mod)
#   })) %>%
#   unnest(estimate)
#
# variance_part <- itv_output %>%
#   # select important columns: sumsq = SS
#   select(trait_trans, mean, term, sumsq) %>%
#   # make wide table
#   pivot_wider(names_from = mean, values_from = sumsq) %>%
#   # rename columns
#   rename("total" = mean, "turnover" = mean_noitv, "intraspecific" = diff) %>%
#   ungroup() %>%
#   # calculate covariation
#   mutate(covariation = total - turnover - intraspecific,
#          # calculate proportion explained variation
#          turnover_p = turnover/total,
#          intra_p = intraspecific/total) %>%
#   # make long table
#   pivot_longer(cols = c(turnover_p:intra_p), names_to = "process", values_to = "value")
#
#
# ITV <- fancy_trait_name_dictionary(variance_part) %>%
#   mutate(process = recode(process, turnover_p = "turnover", intra_p = "ITV")) %>%
#   ggplot(aes(x = trait_fancy, y = value, fill = process)) +
#   geom_col() +
#   coord_flip() +
#   scale_x_discrete(limits = rev) +
#   scale_fill_viridis_d(begin = 0.25, end = 1, option = "viridis") +
#   labs(x = "", y = "% explained") +
#   theme_minimal()
# ggsave(ITV, filename = "ITV.jpg", dpi = 300, height = 6, width = 8)
