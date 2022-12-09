#### TRAIT FIGURES ####

# trait figure
make_trait_figure <- function(community_model_output){

  dat <- community_model_output |>
    # merge data and prediction
    mutate(output = map2(.x = data, .y = prediction, ~ bind_cols(.x, .y))) |>
    select(-data, -singular, -aic, -prediction, -mod, -model_output, -r) |>
    unnest(output) %>%
    rename(Gradient = Gradient...1, mean = .response...4, Elevation_m = .continous_predictor...7, fitted = .response...18) |>
    select(-Gradient...16, -.continous_predictor...17 ) %>%
    fancy_trait_name_dictionary(.) |>
    mutate(class = recode(class, "Leaf economics" = "LES", "Nutrient cycling" = "NS"),
           trait_fancy = paste(class, trait_fancy, sep = " - "),
           trait_fancy = factor(trait_fancy, levels = c("Size - Height cm", "Size - Dry mass g", "Size - Area cm2", "Size - Thickness mm", "LES - SLA cm2/g", "LES - LDMC", "LES - C %", "LES - N %", "LES - CN", "LES - P %", "LES - NP", "NS - δC13 ‰", "NS - δN15 ‰")))

  ggplot(dat, aes(x = Elevation_m, y = mean, colour = Gradient)) +
    geom_point(alpha = 0.5) +
    geom_line(aes(y = fitted, colour = Gradient)) +
    geom_ribbon(aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_colour_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
    scale_fill_manual(name = "", values = c("grey", "green4"), labels = c("Reference", "Nutrient input")) +
    labs(x = "Elevation m a.s.l.", y = "Bootstrapped trait mean") +
    # add label
    geom_text(data = dat %>%
                ungroup() |>
                distinct(trait_trans, trait_fancy, text),
              aes(x = Inf, y = Inf, label = text),
              size = 3, colour = "black", hjust = 1, vjust = 1) +
    facet_wrap(~ trait_fancy, scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "top")

}
