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
           Gradient = recode(Gradient, B = "Nutrient input", C = "Reference"),
           var = "Total") |>
    mutate(new = paste(class, Gradient, sep = "_"),
           new = factor(new, levels = c("Size_Reference", "Size_Nutrient input", "Leaf economics_Reference", "Leaf economics_Nutrient input", "Isotopes_Reference", "Isotopes_Nutrient input"))) |>
    # make turnover negative
    mutate(proportion_standardized = if_else(process == "turnover", -1*proportion_standardized, proportion_standardized)) |>
    ggplot(aes(x = class, y = proportion_standardized, fill = process)) +
    #ggplot(aes(x = new, y = proportion_standardized, fill = process, alpha = Gradient)) +
    geom_col() +
    geom_hline(yintercept = 0, colour = "grey", linetype = "dashed") +
    scale_x_discrete(limits = rev) +
    scale_x_discrete(labels = c("", "Isotopes", "", "Leaf economics", "", "Size")) +
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
    # filter for processes (turnover and ITV) we are interested in and standardize to 1
    filter(process %in% c("turnover", "intraspecific")) |>
    group_by(Gradient, trait_trans) |>
    mutate(sum = sum(proportion),
           proportion_standardized = proportion / sum) |>
    ungroup() |>
    mutate(process = recode(process, intraspecific = "ITV"),
           Gradient = recode(Gradient, B = "Nutrient input", C = "Reference")) |>
    # make turnover negative
    mutate(proportion_standardized = if_else(process == "turnover", -1*proportion_standardized, proportion_standardized)) |>
    # mutate(new = paste(trait_fancy, Gradient, sep = "_"),
    #        new = factor(new, levels = c("Height cm_Reference", "Height cm_Nutrient input",
    #                                     "Dry mass g_Reference", "Dry mass g_Nutrient input",
    #                                     "Area cm2_Reference", "Area cm2_Nutrient input",
    #                                     "Thickness mm_Reference", "Thickness mm_Nutrient input",
    #
    #                                     "LDMC_Reference", "LDMC_Nutrient input",
    #                                     "SLA cm2/g_Reference", "SLA cm2/g_Nutrient input",
    #                                     "C %_Reference", "C %_Nutrient input",
    #                                     "N %_Reference", "N %_Nutrient input",
    #                                     "CN_Reference", "CN_Nutrient input",
    #                                     "P %_Reference", "P %_Nutrient input",
    #                                     "NP_Reference", "NP_Nutrient input",
    #
    #                                     "δC13 ‰_Reference", "δC13 ‰_Nutrient input",
    #                                     "δN15 ‰_Reference", "δN15 ‰_Nutrient input"))) |>
    ggplot(aes(x = trait_fancy, y = proportion_standardized, fill = process)) +
    #ggplot(aes(x = new, y = proportion_standardized, fill = process, alpha = Gradient)) +
    geom_col() +
    geom_hline(yintercept = 0, colour = "grey", linetype = "dashed") +
    scale_x_discrete(limits = rev) +
    # scale_x_discrete(labels = c("Height cm", "",
    #                             "Dry mass g", "",
    #                             "Area cm2", "",
    #                             "Thickness mm", "",
    #
    #                             "LDMC", "",
    #                             "SLA cm2/g", "",
    #                             "C e", "",
    #                             "N %", "",
    #                             "CN", "",
    #                             "P %", "",
    #                             "NP", "",
    #
    #                             "δC13 ‰", "",
    #                             "δN15 ‰", "")) +
    coord_flip() +
    scale_fill_manual(name = "Process", values = c("#005BBB", "#FFD500")) +
    scale_alpha_manual(values = c(1, 0.5)) +
    lims(y = c(-1, 1)) +
    labs(x = "", y = "Relative contribution",
         tag = "(b)") +
    facet_grid(class ~ Gradient, scales = "free", space = "free_y") +
    #facet_grid(class ~ 1, scales = "free", space = "free_y") +
    theme_minimal() +
    theme(text = element_text(size = 18),
          panel.spacing = unit(0.3, "cm"),
          strip.text = element_blank())

  ITV_plot2 <- Group_plot / ITV_plot + plot_layout(guides = 'collect', heights = c(1, 4)) & theme(legend.position = 'top')

  return(ITV_plot2)

}




