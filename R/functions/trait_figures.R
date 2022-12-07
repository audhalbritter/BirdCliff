#### TRAIT FIGURES ####

# trait figure prep
trait_fig_prep <- function(community_trait_model){

  community_trait_model |>
    # merge data and prediction
    mutate(output = map2(.x = data, .y = prediction, ~ bind_cols(.x, .y))) |>
    select(-data, -prediction, -model, -model_output, -r) |>
    unnest(output) |>
  rename(mean = mean...4, fitted = mean...20, Gradient = Gradient...1, Elevation_m = Elevation_m...7) |>
  select(-Elevation_m...19, -Gradient...18)

}


# trait figure
make_trait_figure <- function(community_trait_model){

  community_trait_model |>
    # remove lines for NULL model
    # mutate(fitted = if_else(text %in% c("NULL", ""), NA_real_, fitted),
    #        plo = if_else(text %in% c("NULL", ""), NA_real_, plo),
    #        phi = if_else(text %in% c("NULL", ""), NA_real_, phi)) |>
    ggplot(aes(x = Elevation_m, y = mean, colour = Gradient)) +
    geom_point(alpha = 0.5) +
    geom_line(aes(y = fitted, colour = Gradient)) +
    geom_ribbon(aes(ymin = plo, ymax = phi, fill = Gradient), alpha = 0.3, linetype = 0) +
    scale_colour_manual(name = "", values = c("green4", "grey"), labels = c("Nutrient input", "Reference")) +
    scale_fill_manual(name = "", values = c("green4", "grey"), labels = c("Nutrient input", "Reference")) +
    labs(x = "Elevation m a.s.l.", y = "Bootstrapped trait mean") +
    # add label
    geom_text(data = community_trait_model |>
                ungroup() |>
                distinct(trait_trans, trait_fancy, text),
              aes(x = Inf, y = Inf, label = text),
              size = 3, colour = "black", hjust = 1, vjust = 1) +
    facet_wrap(~ trait_fancy, scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "top")

}
