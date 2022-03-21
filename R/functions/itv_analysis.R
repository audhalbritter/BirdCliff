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
      mod <- aov(value ~ 1, data =  .x)
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
    group_by(Gradient, trait_trans) %>%
    mutate(total_var = sum(total),
           total_turnover = sum(turnover),
           total_ITV = sum(intraspecific)) %>%
    ungroup() %>%
    # calculate covariation
    mutate(covariation = total - turnover - intraspecific,
           # calculate proportion explained variation
           turnover_p = turnover/total,
           intra_p = intraspecific/total) %>%
    # make long table
    pivot_longer(cols = c(turnover_p:intra_p), names_to = "process", values_to = "value") |>
    group_by(Gradient, trait_trans) |>
    mutate(sum = sum(value),
           proportion = value / sum)


  ITV_plot <- fancy_trait_name_dictionary(variance_part) %>%
    # filter only turnover and ITV
    #filter(process %in% c("turnover_p", "intra_p")) %>%
    mutate(process = recode(process, turnover_p = "turnover", intra_p = "ITV"),
           Gradient = recode(Gradient, B = "Bird cliff", C = "Reference")
           ) %>%
    ggplot(aes(x = trait_fancy, y = proportion, fill = process)) +
    geom_col() +
    geom_hline(yintercept = 0.5, colour = "grey", linetype = "dashed") +
    scale_x_discrete(limits = rev) +
    coord_flip() +
    scale_fill_viridis_d(begin = 0.25, end = 1, option = "viridis") +
    labs(x = "", y = "relative contribution") +
    facet_wrap(~ Gradient, scales = "free_x") +
    theme_minimal()

  return(ITV_plot)

}

