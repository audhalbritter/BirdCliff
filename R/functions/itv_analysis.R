### ITV ANALYSIS ###

# total <- specif.avg <- mean (local mean)
# turnover <- const.avg <- mean_noitv (global mean)
# intraspecific <- specific.avg - const.avg

#tar_load(trait_mean)

make_ITV_analysis <- function(trait_mean){

  trait_long <- trait_mean %>%
    # remove nutrient ratio traits
    #filter(!trait_trans %in% c("CN_ratio", "NP_ratio")) |>
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
    rename("total_ss" = mean, "turnover_ss" = mean_noitv, "intraspecific_ss" = diff) %>%
    # calculate covariation SS
    mutate(covariation_ss = total_ss- turnover_ss - intraspecific_ss) |>

    # calculate proportion explained variation (divide ss by total_ss)
    mutate(total_p = total_ss/total_ss,
           turnover_p = turnover_ss/total_ss,
           intraspecific_p = intraspecific_ss/total_ss,
           covariation_p = covariation_ss/total_ss) %>%

    # make long table
    pivot_longer(cols = c(total_ss:covariation_p), names_to = c("process", "variable"), names_sep = "_",values_to = "value") |>
    mutate(variable = recode(variable, "ss" = "sumsq", "p" = "proportion")) |>
    pivot_wider(names_from = variable, values_from = value)

  # proportion ITV
  # variance_part |>
  #   filter(process  == "intraspecific") |>
  #   ungroup() |>
  #   group_by(Gradient) |>
  #   summarise(range(proportion))

  # write table with results
  fancy_trait_name_dictionary(variance_part) |>
    mutate(sumsq = round(sumsq, digits = 1),
           proportion = round(proportion, digits = 1)) |>
    ungroup() |>
    select(Gradient, trait_fancy, term:proportion) |>
    write_csv(, file = "output/ITV_output.csv")

  # test difference
  # dd <- variance_part |>
  #   filter(process %in% c("turnover", "intraspecific")) |>
  #   select(Gradient, trait_trans, process, proportion) |>
  #   #pivot_wider(names_from = process, values_from = proportion) |>
  #   mutate(type = if_else(trait_trans %in% c("Plant_Height_cm_log", "Dry_Mass_g_log", "Leaf_Area_cm2_log", "Thickness_mm_log", "LDMC", "SLA_cm2_g"), "size", "nutrient"))
  #
  #
  # dd |>
  #   ungroup() |>
  #   group_by(Gradient, type) %>%
  #   nest() %>%
  #   mutate(test = map(data, ~{
  #     mod = aov(proportion ~ process, data = .)
  #     result = tidy(mod)
  #   })) |>
  #   unnest(test)

  Group_plot <- fancy_trait_name_dictionary(variance_part) %>%
    # filter for processes (turnover and ITV) we are interested in and standardize to 1
    filter(process %in% c("turnover", "intraspecific")) |>
    group_by(Gradient, trait_trans) |>
    mutate(sum = sum(proportion),
           proportion_standardized = proportion / sum) |>
    ungroup() |>
    group_by(Gradient, class, process) |>
    summarise(proportion_standardized = mean(proportion_standardized)) |>
    mutate(process = recode(process, intraspecific = "ITV"),
           Gradient = recode(Gradient, B = "Nutrient", C = "Reference"),
           var = "Total") |>
    mutate(new = paste(class, Gradient, sep = "_"),
           new = factor(new, levels = c("Size_Reference", "Size_Nutrient", "Leaf economics_Reference", "Leaf economics_Nutrient", "Isotopes_Reference", "Isotopes_Nutrient"))) |>
    # make turnover negative
    mutate(proportion_standardized = if_else(process == "turnover", -1*proportion_standardized, proportion_standardized)) |>
    ggplot(aes(x = class, y = proportion_standardized, fill = process)) +
    geom_col() +
    geom_hline(yintercept = 0, colour = "grey", linetype = "dashed") +
    scale_x_discrete(limits = rev) +
    coord_flip() +
    scale_fill_manual(name = "Process", values = c("#005BBB", "#FFD500")) +
    scale_alpha_manual(values = c(1, 0.5)) +
    lims(y = c(-1, 1)) +
    labs(x = "", y = "Relative contribution",
         tag = "(a)") +
    facet_wrap( ~ Gradient, scales = "free_x") +
    theme_minimal() +
    theme(text = element_text(size = 18),
          panel.spacing = unit(1, "cm"))


  ITV_plot <- fancy_trait_name_dictionary(variance_part) %>%
    # remove class part
    mutate(figure_names = str_remove(figure_names, "Size~-~|LES~-~|I~-~")) |>
    mutate(figure_names = factor(figure_names, levels = c("Height~cm", "Dry~mass~g", "Area~cm^2", "Thickness~mm", "SLA~cm^2*g^{-1}", "LDMC", "C~'%'", "N~'%'", "CN", "P~'%'", "NP", "δ^{13}~C~'‰'", "δ^{15}~N~'‰'"))) |>
    # filter for processes (turnover and ITV) we are interested in and standardize to 1
    filter(process %in% c("turnover", "intraspecific")) |>
    group_by(Gradient, trait_trans) |>
    mutate(sum = sum(proportion),
           proportion_standardized = proportion / sum) |>
    ungroup() |>
    mutate(process = recode(process, intraspecific = "ITV"),
           Gradient = recode(Gradient, B = "Nutrient", C = "Reference")) |>
    # make turnover negative
    mutate(proportion_standardized = if_else(process == "turnover", -1*proportion_standardized, proportion_standardized)) |>
    ggplot(aes(x = figure_names, y = proportion_standardized, fill = process)) +
    geom_col() +
    geom_hline(yintercept = 0, colour = "grey", linetype = "dashed") +
    scale_x_discrete(limits = rev, labels = ggplot2:::parse_safe) +
    coord_flip() +
    scale_fill_manual(name = "Process", values = c("#005BBB", "#FFD500")) +
    scale_alpha_manual(values = c(1, 0.5)) +
    lims(y = c(-1, 1)) +
    labs(x = "", y = "Relative contribution",
         tag = "(b)") +
    facet_grid(class ~ Gradient, scales = "free", space = "free_y") +
    theme_minimal() +
    theme(text = element_text(size = 18),
          panel.spacing = unit(0.3, "cm"),
          strip.text = element_blank())

  ITV_plot2 <- Group_plot / ITV_plot + plot_layout(guides = 'collect', heights = c(1, 4)) & theme(legend.position = 'top')

  return(ITV_plot2)

}




