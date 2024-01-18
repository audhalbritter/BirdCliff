#### TRAIT FIGURES ####

# trait figure
make_trait_figure <- function(community_model_output){

  dat <- community_model_output |>
    # merge data and prediction
    mutate(output = map2(.x = data, .y = prediction, ~ bind_cols(.x, .y))) |>
    select(-data, -singular, -aic, -prediction, -mod, -model_output, -r) |>
    unnest(output) %>%
    rename(Gradient = Gradient...1, mean = .response...4, Elevation_m = .continous_predictor...9, fitted = .response...19) |>
    select(-Gradient...17, -.continous_predictor...18) %>%
    fancy_trait_name_dictionary(.) %>%
    mutate(figure_names = factor(figure_names, levels = c("Size~-~Height~(cm)", "Size~-~Dry~mass~(g)", "Size~-~Area~(cm^2)", "Size~-~Thickness~(mm)", "LES~-~SLA~(cm^2*g^{-1})", "LES~-~LDMC", "LES~-~C~('%')", "LES~-~N~('%')", "LES~-~CN", "LES~-~P~('%')", "LES~-~NP", "I~-~δ^{13}~C~'(‰)'", "I~-~δ^{15}~N~'(‰)'")))



  ggplot(dat, aes(x = Elevation_m, y = mean, colour = Gradient)) +
    geom_point(alpha = 0.5) +
    geom_line(aes(y = fitted, colour = Gradient)) +
    geom_ribbon(aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_colour_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient")) +
    scale_fill_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient")) +
    labs(x = "Elevation m a.s.l.", y = "Bootstrapped trait mean") +
    # add label
    geom_text(data = dat %>%
                ungroup() |>
                distinct(trait_trans, figure_names, text),
              aes(x = Inf, y = Inf, label = text),
              size = 3, colour = "black", hjust = 1, vjust = 1) +
    facet_wrap(~ figure_names, scales = "free_y", labeller = label_parsed) +
    theme_minimal() +
    theme(legend.position = "top")

}

make_dN15_figure <- function(dN15_model_output){

  # prep data
  dat <- dN15_model_output |>
    # merge data and prediction
    mutate(output = map2(.x = data, .y = prediction, ~ bind_cols(.x, .y))) |>
    select(-data, -singular, -prediction, -mod, -model_output, -r) |>
    unnest(output) %>%
    rename(Gradient = Gradient...1, mean = .response...6, dN15 = .continous_predictor...5, fitted = .response...9) |>
    select(-Gradient...7, -.continous_predictor...8) %>%
    fancy_trait_name_dictionary(.) %>%
    mutate(figure_names = factor(figure_names, levels = c("Size~-~Height~(cm)", "Size~-~Dry~mass~(g)", "Size~-~Area~(cm^2)", "Size~-~Thickness~(mm)", "LES~-~SLA~(cm^2*g^{-1})", "LES~-~LDMC", "LES~-~C~('%')", "LES~-~N~('%')", "LES~-~CN", "LES~-~P~('%')", "LES~-~NP", "I~-~δ^{13}~C~'(‰)'", "I~-~δ^{15}~N~'(‰)'"))) |>
    # fix stats
    mutate(text = case_match(text,
                             "δN15" ~ "δ^{15}~N",
                             "N+δN15" ~ "N+δ^{15}~N",
                             "NxδN15" ~ "Nxδ^{15}~N",
                             "Null" ~ "Null"))

  # make figure
  dn15_figure <- ggplot(dat, aes(x = dN15, y = mean, colour = Gradient)) +
    geom_point(alpha = 0.5) +
    geom_line(aes(y = fitted, colour = Gradient)) +
    geom_ribbon(aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_colour_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient")) +
    scale_fill_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient")) +
    labs(x = bquote(δ^'15'~"N trait mean"), y = "Bootstrapped trait mean") +
    # add label
    geom_text(data = dat %>%
                ungroup() |>
                distinct(trait_trans, figure_names, text),
              aes(x = Inf, y = Inf, label = text),
              size = 3, colour = "black", hjust = 1, vjust = 1, parse = TRUE) +
    facet_wrap(~ figure_names, scales = "free_y", labeller = label_parsed) +
    theme_minimal() +
    theme(legend.position = "top")

}




# # trait figure
# make_trait_variance_figure <- function(community_variance_output){
#
#   dat <- community_variance_output |>
#     # merge data and prediction
#     mutate(output = map2(.x = data, .y = prediction, ~ bind_cols(.x, .y))) |>
#     select(-data, -singular, -aic, -prediction, -mod, -model_output, -r) |>
#     unnest(output) %>%
#     rename(Gradient = Gradient...1, var = .response...5, Elevation_m = .continous_predictor...9, fitted = .response...19) |>
#     select(-Gradient...17, -.continous_predictor...18) %>%
#     fancy_trait_name_dictionary(.) |>
#     mutate(class = recode(class, "Leaf economics" = "LES", "Isotopes" = "I"),
#            trait_fancy = paste(class, trait_fancy, sep = " - "),
#            trait_fancy = factor(trait_fancy, levels = c("Size - Height cm", "Size - Dry mass g", "Size - Area cm2", "Size - Thickness mm", "LES - SLA cm2/g", "LES - LDMC", "LES - C %", "LES - N %", "LES - CN", "LES - P %", "I - δC13 ‰", "I - δN15 ‰")))
#
#   ggplot(dat, aes(x = Elevation_m, y = var, colour = Gradient)) +
#     geom_point(alpha = 0.5) +
#     geom_line(aes(y = fitted, colour = Gradient)) +
#     geom_ribbon(aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
#     scale_colour_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
#     scale_fill_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
#     labs(x = "Elevation m a.s.l.", y = "Bootstrapped trait variance") +
#     # add label
#     geom_text(data = dat %>%
#                 ungroup() |>
#                 distinct(trait_trans, trait_fancy, text),
#               aes(x = Inf, y = Inf, label = text),
#               size = 3, colour = "black", hjust = 1, vjust = 1) +
#     facet_wrap(~ trait_fancy, scales = "free_y") +
#     theme_minimal() +
#     theme(legend.position = "top")
#
# }
