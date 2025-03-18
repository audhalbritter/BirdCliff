# comm figure

make_community_figure <- function(comm_structure, comm_raw, family_list){

  dat <- comm_structure |>
    filter(!Variable %in% c("MedianHeight_cm", "MaxHeight_cm")) |>
    mutate(Variable = case_when(grepl("Lichen", Variable)  ~ "Lichen",
                                Variable == "Vascular" ~ "Vascular plants",
                                Variable == "BioCrust" ~ "Biocrust",
                                Variable == "BareGround" ~ "Bare ground",
                                TRUE ~ Variable),
           Variable = factor(Variable, levels = c("Rock", "Bare ground", "Biocrust", "Litter", "Lichen",  "Bryophytes", "Vascular plants")),
           Gradient = if_else(Gradient == "C", "Reference", "Nutrient"),
           Gradient = factor(Gradient, levels = c("Reference", "Nutrient")))
    # make reference gradient negative
    #mutate(Value = if_else(Gradient == "Reference", -1*Value, Value))

comm1 <- ggplot(dat, aes(x = Site, y = Value, fill = Variable)) +
    geom_col(position = "fill") +
    coord_flip() +
    scale_fill_manual(values = hcl.colors(palette = "Batlow", n = 7, rev = TRUE), name = "") +
    labs(x = "Site",
         y = "Relative cover",
         tag = "a)") +
    facet_wrap(~ Gradient) +
    theme_minimal()



  sums <- comm_raw |>
    tidylog::left_join(family_list, by = "Taxon") |>
    group_by(Gradient, FunctionalGroup, Site) |>
    summarise(sumofCover = sum(Cover)) |>
    mutate(Gradient = if_else(Gradient == "C", "Reference", "Nutrient"),
           Gradient = factor(Gradient, levels = c("Reference", "Nutrient")),
           FunctionalGroup = case_when(FunctionalGroup == "dshrub" ~ "Decidious shrubs",
                                       FunctionalGroup == "eshrub" ~ "Evergreen shrubs",
                                       FunctionalGroup == "forbs" ~ "Forbs",
                                       FunctionalGroup == "graminoid" ~ "Graminoids",
                                       is.na(FunctionalGroup) ~ "unknown",
                                       TRUE ~ FunctionalGroup),
           FunctionalGroup = factor(FunctionalGroup, levels = c("Graminoids", "Forbs", "Decidious shrubs", "Evergreen shrubs", "unknown"))) |>
    filter(!FunctionalGroup %in% c("unknown"))
    # make reference gradient negative
    #mutate(sumofCover = if_else(Gradient == "Reference", -1*sumofCover, sumofCover))

  comm2 <- ggplot(sums, aes(x = Site, y = sumofCover, fill = FunctionalGroup)) +
    geom_col(position = "fill") +
    coord_flip() +
    scale_fill_manual(values = hcl.colors(palette = "Hawaii", n = 4, rev = TRUE), name = "") +
    labs(x = "Sites",
         y = "Relative cover",
         tag = "b)") +
    facet_wrap(~ Gradient) +
    theme_minimal()

  comm1 / comm2

}


