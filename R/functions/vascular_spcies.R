# traits_raw %>% distinct(Taxon)
#
# # choosing species and traits
# traits_raw %>%
#   distinct(Gradient, Site, Taxon) %>%
#   group_by(Gradient, Taxon) %>%
#   count() %>%
#   filter(n > 3)
#
# bryo_traits_raw %>%
#   distinct(Gradient, Site, Taxon) %>%
#   group_by(Gradient, Taxon) %>%
#   count()


make_ind_sp_plot <- function(traits_raw, bryo_traits_raw){

  # cobmine species vascular and bryophytes
  all_traits <- bind_rows(traits_raw %>%
                            filter(Taxon %in% c("luzula confusa", "salix polaris", "cerastium arcticum"),
                                   trait_trans %in% c("Plant_Height_cm_log", "Leaf_Area_cm2_log", "LDMC", "N_percent", "dN15_permil")),
                          bryo_traits_raw %>%
                            filter(Taxon %in% c("aulacomnium turgidum", "hylocomium splendens", "sanionia uncinata"),
                                   trait_trans %in% c("Shoot_Length_cm_log", "Shoot_Length_Green_cm_log", "Shoot_ratio", "WHC_g_g", "P_percent"))) %>%
    mutate(trait_trans = factor(trait_trans, levels = c("Plant_Height_cm_log", "Leaf_Area_cm2_log", "LDMC", "N_percent", "dN15_permil", "Shoot_Length_cm_log", "Shoot_Length_Green_cm_log", "Shoot_ratio", "WHC_g_g", "P_percent")),
           Taxon = factor(Taxon, levels = c("cerastium arcticum", "luzula confusa", "salix polaris", "aulacomnium turgidum", "hylocomium splendens", "sanionia uncinata")))


  v <- fancy_trait_name_dictionary(all_traits) %>%
    filter(Functional_group == "vascular") %>%
    ggplot(aes(x = Elevation_m, y = value_trans, colour = Gradient)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE) +
    scale_colour_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
    facet_grid(trait_fancy ~ Taxon, scales = "free_y") +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(aspect.ratio = 0.8,
          plot.margin = margin(0.5, 0.5, 0.5, 0.5))


  b <- fancy_trait_name_dictionary(all_traits) %>%
    filter(Functional_group == "bryophyte") %>%
    ggplot(aes(x = Elevation_m, y = value_trans, colour = Gradient)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE) +
    scale_colour_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
    facet_grid(trait_fancy ~ Taxon, scales = "free_y") +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(aspect.ratio = 0.8,
          plot.margin = margin(0.5, 0.5, 0.5, 0.5))

  ind_sp_traits <- v + b + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

  return(ind_sp_traits)

}
