# ORDINATIONS

## COMMUNITY (PCA)

make_community_pca <- function(comm_raw){

  set.seed(32)

  comm_wide <- comm_raw %>%
    #filter(!c(GS == "B3" & PlotID == "C")) |>
    group_by(Site) %>%
    mutate(Mean_elevation = mean(Elevation_m)) %>%
    ungroup() %>%
    mutate(Cover = sqrt(Cover)) |>
    select(Gradient, Site, GS, Mean_elevation, PlotID, Taxon, Cover) %>%
    mutate(Gradient = recode(Gradient, "B" = "Nutrient", "C" = "Reference")) %>%
    pivot_wider(names_from = Taxon,
                values_from = Cover,
                values_fill = 0)

  comm_sp <- comm_wide %>%
    select(-c(Gradient:PlotID))

  # meta data
  comm_info <- comm_wide %>%
    select(Gradient:PlotID)

  # make pca
  res <- rda(comm_sp)

  out <- bind_cols(comm_info, fortify(res) |>
                     filter(score == "sites"))

  sp <- fortify(res) |>
    filter(score == "species")

  # adonis test
  adonis <- adonis(comm_sp ~ Gradient * Mean_elevation , data = comm_info, permutations = 999, method = "euclidean")

  return(list(out, sp, res, adonis))

}


make_sp_pca_figure <- function(comm_pca){

  # bird <- grid::rasterGrob(png::readPNG("bird.png"), interpolate = TRUE)
  # ref <- grid::rasterGrob(png::readPNG("ref.png"), interpolate = TRUE)

  # elevational range
  range <- range(comm_pca[[1]]$Mean_elevation)

  # PC range
  PC1_min <- min(comm_pca[[1]]$PC1) - 1
  PC1_max <- max(comm_pca[[1]]$PC1) + 1
  PC2_min <- min(comm_pca[[1]]$PC2) - 1
  PC2_max <- max(comm_pca[[1]]$PC2) + 1

  e_B <- eigenvals(comm_pca[[3]])/sum(eigenvals(comm_pca[[3]]))

  important_species <- comm_pca[[2]] |>
    mutate(length = sqrt(PC1^2 + PC2^2),
           label = capitalize(label)) |>
    filter(length > 0.5) |>
    select(label, length, PC1, PC2)

  # make main figure
  # reference
  pca_C <- comm_pca[[1]] |>
    filter(Gradient == "Reference") |>
    ggplot(aes(x = PC1, y = PC2, colour = Mean_elevation)) +
    geom_point(size = 2) +
    #annotation_custom(ref, xmin = -3.5, xmax = -5.3, ymin = 1, ymax = 2.1) +
    coord_equal(xlim = c(PC1_min, PC1_max),
                ylim = c(PC2_min, PC2_max)) +
    labs(x = glue("PCA1 ({round(e_B[1] * 100, 1)}%)"),
         y = glue("PCA2 ({round(e_B[2] * 100, 1)}%)"),
         title = "Reference",
         tag = "(a)") +
    stat_ellipse(aes(colour = Mean_elevation, group = GS)) +
    scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1, name = "Elevation m a.s.l.", limits = c(range[1], range[2])) +
    scale_shape_manual(values = c(16, 1)) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    theme_bw() +
    theme(text = element_text(size = 13),
          legend.box="vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio = 1)

  # Nutrient gradient
  pca_B <- comm_pca[[1]] |>
    filter(Gradient == "Nutrient") |>
    ggplot(aes(x = PC1, y = PC2, colour = Mean_elevation)) +
    geom_point(size = 2) +
    #annotation_custom(bird, xmin = -3.5, xmax = -4.5, ymin = 1.3, ymax = 2.6) +
    coord_equal(xlim = c(PC1_min, PC1_max),
                ylim = c(PC2_min, PC2_max)) +
    labs(x = glue("PCA1 ({round(e_B[1] * 100, 1)}%)"),
         y = glue("PCA2 ({round(e_B[2] * 100, 1)}%)"),
         title = "Nutrient",
         tag = "(b)") +
    stat_ellipse(aes(colour = Mean_elevation, group = GS)) +
    scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1, name = "Elevation m a.s.l.", limits = c(range[1], range[2])) +
    scale_shape_manual(values = c(16, 1)) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    theme_bw() +
    theme(text = element_text(size = 13),
          legend.box="vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio = 1)

  # species names
  species <- comm_pca[[1]] |>
    ggplot(aes(x = PC1, y = PC2)) +
    coord_equal(xlim = c(-3.4, PC1_max),
                ylim = c(PC2_min, PC2_max)) +
    geom_segment(data = comm_pca[[2]] |>
                   mutate(length = sqrt(PC1^2 + PC2^2)),
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2,"cm")),
                 alpha = 0.75,
                 color = 'grey70') +
    geom_text(data = comm_pca[[2]] |>
                mutate(label = capitalize(label)) |>
                inner_join(important_species) |>
                mutate(PC1 = case_when(label == "Salix polaris" ~ -1.5,
                                       label == "Oxyria digyna" ~ 1.2,
                                       label == "Cerastium arcticum" ~ 0.5,
                                       label == "Bistorta vivipara" ~ 0.5,
                                       label == "Dryas octopetala" ~ 0.5,
                                        TRUE ~ PC1),
                       PC2 = case_when(label == "Salix polaris" ~ 0.2,
                                       label == "Oxyria digyna" ~ 0.4,
                                       label == "Cerastium arcticum" ~ 0.5,
                                       label == "Bistorta vivipara" ~ -1.3,
                                       label == "Dryas octopetala" ~ -2,
                                       TRUE ~ PC2)),
              aes(x = PC1, y = PC2, label = label),
              size = 2, col = 'black') +
    labs(x = "PC 1", y = "PC 2", tag = "(c)") +
    theme_bw() +
    theme(text = element_text(size = 13),
          legend.box="vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio = 1)

  plot_layout <- "
  12
  34
  "

  pca_C + pca_B + species + guide_area() + plot_layout(design = plot_layout, guides = "collect") &
    theme(legend.box = "horizontal")

}


## TRAITS (PCA)
make_trait_pca <- function(trait_mean){

  set.seed(32)

  # make wide trait table
  cwm_fat <- trait_mean %>%
    # remove nutrient ratio traits
    filter(!trait_trans %in% c("CN_ratio", "NP_ratio")) |>
    group_by(Site) %>%
    mutate(Mean_elevation = mean(Elevation_m),
           GS = paste0(Gradient, Site)) %>%
    select(Gradient:mean, Mean_elevation, GS, SoilMoisture, SoilTemperature) %>%
    pivot_wider(names_from = "trait_trans", values_from = "mean") %>%
    ungroup()

  # environemntal variables
  env <- cwm_fat |>
    select(SoilMoisture, SoilTemperature)

  pca_output <- cwm_fat %>%
    select(-(Gradient:SoilTemperature)) %>%
    rda(scale = TRUE, center = TRUE)

  # envfit
  env_fit <- envfit(pca_output, env)
  env_out <- scores(env_fit, "vectors") |>
    as_tibble() |>
    mutate(label = c("Moisture", "Temperature"),
           class = "Environment",
           figure_names = c("Moisture", "Temperature"))

  pca_sites <- bind_cols(
    cwm_fat %>%
      select(Gradient:SoilTemperature),
    fortify(pca_output, display = "sites")
  )

  # arrows
  pca_traits <- fortify(pca_output, display = "species") %>%
    mutate(trait_trans = label) %>%
    fancy_trait_name_dictionary() |>
    mutate(figure_names = str_remove(figure_names, "Size~-~|LES~-~|I~-~"),
           figure_names = str_extract(figure_names, "[^~]+"),
           figure_names = recode(figure_names,
                                     "δ^{13}" = "δ^{13}~C",
                                     "δ^{15}" = "δ^{15}~N")) |>
    # add environmental variables
    bind_rows(env_out) |>
    mutate(class = as.character(class),
           class = factor(class, levels = c("Size", "Leaf economics", "Isotopes", "Environment")))

  # permutation test
  # traits
  raw <- cwm_fat %>% select(-(Gradient:SoilTemperature))
  # meta data
  meta <- cwm_fat %>% select(Gradient:SoilTemperature) %>%
    mutate(Site = factor(Site))

  # adonis test
  if(meta %>% distinct(Gradient) %>% count() == 2){
    adonis_result <- adonis2(raw ~ Gradient*Mean_elevation , data = meta, permutations = 999, method = "euclidean")
  } else {
    adonis_result <- adonis2(raw ~ Mean_elevation, data = meta, permutations = 999, method = "euclidean")
  }

  outputList <- list(pca_sites, pca_traits, pca_output, adonis_result)

  return(outputList)
}



make_pca_plot <- function(trait_pca){

  # adding pics
  # bird <- grid::rasterGrob(png::readPNG("bird.png"), interpolate = TRUE)
  # ref <- grid::rasterGrob(png::readPNG("ref.png"), interpolate = TRUE)
  # add this to figure
  # annotation_custom(bird, xmin = -2.5, xmax = -3.8, ymin = 1.3, ymax = 2.6)

  # elevational range
  range <- range(trait_pca[[1]]$Mean_elevation)

  # PC range
  PC1_min <- min(trait_pca[[1]]$PC1) - 1
  PC1_max <- max(trait_pca[[1]]$PC1) + 1
  PC2_min <- min(trait_pca[[1]]$PC2) - 1.5
  PC2_max <- max(trait_pca[[1]]$PC2) + 1.5

  # both gradients
  e_B <- eigenvals(trait_pca[[3]])/sum(eigenvals(trait_pca[[3]]))

  # Reference
  pcaC <- trait_pca[[1]] %>%
    filter(Gradient == "C") |>
    ggplot(aes(x = PC1, y = PC2, colour = Mean_elevation, group = GS)) +
    geom_point(size = 2) +
    coord_equal(xlim = c(PC1_min, PC1_max), ylim = c(PC2_min, PC2_max)) +
    #coord_equal(clip = "off", xlim = c(PC1_min, PC1_max), ylim = c(PC2_min, PC2_max)) +
    stat_ellipse(aes(colour = Mean_elevation)) +
    scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1, name = "Elevation m a.s.l.", limits = c(range[1], range[2])) +
    labs(x = glue("PCA1 ({round(e_B[1] * 100, 1)}%)"),
         y = glue("PCA2 ({round(e_B[2] * 100, 1)}%)"),
         tag = "(a) Reference") +
    theme_minimal() +
    theme(aspect.ratio = 1,
          plot.tag.position = c(0, 0.9),
          plot.tag = element_text(vjust = -1.5, hjust = -0.5, size = 10))

  # Nutrient input
  pcaN <- trait_pca[[1]] %>%
    filter(Gradient == "B") |>
    ggplot(aes(x = PC1, y = PC2, colour = Mean_elevation, group = GS)) +
    geom_point(size = 2) +
    coord_equal(xlim = c(PC1_min, PC1_max), ylim = c(PC2_min, PC2_max)) +
    #coord_equal(clip = "off", xlim = c(PC1_min, PC1_max), ylim = c(PC2_min, PC2_max)) +
    stat_ellipse(aes(colour = Mean_elevation)) +
    scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1, name = "Elevation m a.s.l.", limits = c(range[1], range[2])) +
    labs(x = glue("PCA1 ({round(e_B[1] * 100, 1)}%)"),
         y = glue("PCA2 ({round(e_B[2] * 100, 1)}%)"),
         tag = "(b) Nutrient") +
    theme_minimal() +
    theme(aspect.ratio = 1,
          plot.tag.position = c(0, 0.9),
          plot.tag = element_text(vjust = -1.5, hjust = -0.5, size = 10))

  arrow <- trait_pca[[1]] %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_segment(data = trait_pca[[2]],
                 aes(x = 0, y = 0, xend = PC1, yend = PC2, colour = class, linetype = class),
                 arrow = arrow(length = unit(0.2, "cm")),
                 inherit.aes = FALSE) +
    geom_text(data = trait_pca[[2]] |>
                mutate(PC1 = case_when(label == "Thickness_mm_log" ~ 1.23,
                                       label == "LDMC" ~ -1,
                                       label == "Temperature" ~ -0.4,
                                       label == "Moisture" ~ -1,
                                       TRUE ~ PC1),
                       PC2 = case_when(label == "dC13_permil" ~ -1.25,
                                       label == "C_percent" ~ -0.3,
                                       label == "dN15_permil" ~ -0.3,
                                       label == "Temperature" ~ 0.4,
                                       label == "Moisture" ~ 0.06,
                                       TRUE ~ PC2)),
              aes(x = PC1+0.3, y = PC2, label = figure_names, colour = class),
              size = 2.5,
              inherit.aes = FALSE,
              show.legend = FALSE, parse = TRUE) +
    coord_equal(clip = "off", xlim = c(PC1_min+1, PC1_max-1.5), ylim = c(PC2_min+1.5, PC2_max-3)) +
    labs(x = "PCA1", y = "PCA2", tag = "(c)") +
    scale_x_continuous(expand = c(.2, 0)) +
    scale_linetype_manual(name = "", values = c("solid", "dashed", "solid", "solid")) +
    scale_colour_manual(name = "", values = c("black", "grey40", "grey70", "cornflowerblue")) +
    theme_minimal() +
    theme(aspect.ratio = 1,
          plot.tag.position = c(0, 1),
          plot.tag = element_text(vjust = 1.5, hjust = -2.85, size = 10))

  plot_layout <- "
  12
  34
  "

  trait_ordination_plot <- pcaC + pcaN + arrow + guide_area() + plot_layout(design = plot_layout, guides = "collect") &
    theme(legend.box = "horizontal")

  return(trait_ordination_plot)

}

make_pca_new <- function(trait_pca){

  # elevational range
  range <- range(trait_pca[[1]]$Mean_elevation)

  # both gradients
  e_B <- eigenvals(trait_pca[[3]])/sum(eigenvals(trait_pca[[3]]))

  # calculte centroids
  centroids <- trait_pca[[1]] |>
    left_join(trait_pca[[1]] |>
                group_by(Gradient, Site) |>
                summarise(centroid1 = mean(PC1),
                          centroid2 = mean(PC2)),
              by = c("Gradient", "Site"))

  # One PCA with centroids
  pca <- trait_pca[[1]] %>%
    ggplot(aes(x = PC1, y = PC2, colour = Mean_elevation, shape = Gradient)) +
    #geom_point(size = 2) +
    # centroid and lines
    geom_segment(data = centroids,
                 aes(x = centroid1, y = centroid2,
                     xend = PC1, yend = PC2, colour = Mean_elevation,
                     linetype = Gradient),
                 size = 0.6,
                 alpha = 0.5,
                 show.legend = FALSE) +
    geom_path(data = centroids |>
                distinct(Gradient, Site, GS, Mean_elevation, centroid1, centroid2),
              aes(x = centroid1, y = centroid2,
                  linetype = Gradient),
              colour = "grey20") +
    geom_point(data = centroids,
               aes(x = centroid1, y = centroid2,
                   colour = Mean_elevation),
               size = 2.5) +

    coord_equal() +

    #stat_ellipse(aes(colour = Mean_elevation, group = GS, colour = Mean_elevation)) +
    scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1, name = "Elevation m a.s.l.", limits = c(range[1], range[2])) +
    scale_shape_manual(values = c(1, 16), name = "", labels = c("Reference", "Nutrient")) +
    scale_linetype_manual(values = c("dashed", "solid"), name = "", labels = c("Reference", "Nutrient")) +
    labs(x = glue("PCA1 ({round(e_B[1] * 100, 1)}%)"),
         y = glue("PCA2 ({round(e_B[2] * 100, 1)}%)"),
         tag = "(a)") +
    theme_minimal() +
    theme(aspect.ratio = 1,
          plot.tag.position = c(0, 0.9),
          plot.tag = element_text(vjust = -1.5, hjust = -0.5, size = 10))


  arrow <- trait_pca[[1]] %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_segment(data = trait_pca[[2]],
                 aes(x = 0, y = 0, xend = PC1, yend = PC2, colour = class, linetype = class),
                 arrow = arrow(length = unit(0.2, "cm")),
                 inherit.aes = FALSE) +
    geom_text(data = trait_pca[[2]] |>
                mutate(PC1 = case_when(label == "Thickness_mm_log" ~ 1.23,
                                       label == "LDMC" ~ -1,
                                       label == "Temperature" ~ -0.4,
                                       label == "Moisture" ~ -1,
                                       TRUE ~ PC1),
                       PC2 = case_when(label == "dC13_permil" ~ -1.25,
                                       label == "C_percent" ~ -0.3,
                                       label == "dN15_permil" ~ -0.3,
                                       label == "Temperature" ~ 0.4,
                                       label == "Moisture" ~ 0.06,
                                       TRUE ~ PC2)),
              aes(x = PC1+0.3, y = PC2, label = figure_names, colour = class),
              size = 2.5,
              inherit.aes = FALSE,
              show.legend = FALSE, parse = TRUE) +
    coord_equal() +
    #coord_equal(clip = "off", xlim = c(PC1_min+1, PC1_max-1.5), ylim = c(PC2_min+1.5, PC2_max-3)) +
    labs(x = "PCA1", y = "PCA2", tag = "(b)") +
    scale_x_continuous(expand = c(.2, 0)) +
    scale_linetype_manual(name = "", values = c("solid", "dashed", "solid", "solid")) +
    scale_colour_manual(name = "", values = c("black", "grey40", "grey70", "cornflowerblue")) +
    theme_minimal() +
    theme(aspect.ratio = 1,
          plot.tag.position = c(0, 1),
          plot.tag = element_text(vjust = 1.5, hjust = -2.85, size = 10))

  var <- trait_pca[[4]] |>
    tidy() |>
    filter(term != "Total") |>
    mutate(term = case_when(term == "Gradient" ~ "N",
                            term == "Mean_elevation" ~ "E",
                            term == "Gradient:Mean_elevation" ~ "NxE",
                            TRUE ~ term),
           term = factor(term, levels = c("N", "E", "NxE", "Residual"))) |>
    ggplot(aes(x = term, y = R2)) +
    geom_col() +
    labs(x = "",
         y = "Variation explained (%)",
         tag = "(c)") +
    theme_minimal()

  plot_layout <- "
  1122
  1122
  3444
  "

  trait_ordination_plot <- pca + arrow + var + guide_area() + plot_layout(design = plot_layout, guides = "collect") &
    theme(legend.box = "horizontal")


}



make_pca_3_plot <- function(trait_pca){

  # elevational range
  range <- range(trait_pca[[1]]$Mean_elevation)

  # PC range
  PC3_min <- min(trait_pca[[1]]$PC3) - 1
  PC3_max <- max(trait_pca[[1]]$PC3) + 0.5
  PC4_min <- min(trait_pca[[1]]$PC4) - 0.5
  PC4_max <- max(trait_pca[[1]]$PC5) + 1

  # both gradients
  e_B <- eigenvals(trait_pca[[3]])/sum(eigenvals(trait_pca[[3]]))

  # Reference
  pcaC <- trait_pca[[1]] %>%
    filter(Gradient == "C") |>
    ggplot(aes(x = PC3, y = PC4, colour = Mean_elevation, group = GS)) +
    geom_point(size = 2) +
    coord_equal(xlim = c(PC3_min, PC3_max), ylim = c(PC4_min, PC4_max)) +
    stat_ellipse(aes(colour = Mean_elevation)) +
    scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1, name = "Elevation m a.s.l.", limits = c(range[1], range[2])) +
    labs(x = glue("PCA3 ({round(e_B[3] * 100, 1)}%)"),
         y = glue("PCA4 ({round(e_B[4] * 100, 1)}%)"),
         tag = "(a) Reference") +
    theme_minimal() +
    theme(aspect.ratio = 1,
          plot.tag.position = c(0, 0.9),
          plot.tag = element_text(vjust = -1.5, hjust = -0.5, size = 10))

  # Nutrient input
  pcaN <- trait_pca[[1]] %>%
    filter(Gradient == "B") |>
    ggplot(aes(x = PC3, y = PC4, colour = Mean_elevation, group = GS)) +
    geom_point(size = 2) +
    coord_equal(xlim = c(PC3_min, PC3_max), ylim = c(PC4_min, PC4_max)) +
    stat_ellipse(aes(colour = Mean_elevation)) +
    scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1, name = "Elevation m a.s.l.", limits = c(range[1], range[2])) +
    labs(x = glue("PCA3 ({round(e_B[3] * 100, 1)}%)"),
         y = glue("PCA4 ({round(e_B[4] * 100, 1)}%)"),
         tag = "(b) Nutrient") +
    theme_minimal() +
    theme(aspect.ratio = 1,
          plot.tag.position = c(0, 0.9),
          plot.tag = element_text(vjust = -1.5, hjust = -0.5, size = 10))

  arrow <- trait_pca[[1]] %>%
    ggplot(aes(x = PC3, y = PC4)) +
    geom_segment(data = trait_pca[[2]],
                 aes(x = 0, y = 0, xend = PC3, yend = PC4, colour = class, linetype = class),
                 arrow = arrow(length = unit(0.2, "cm")),
                 inherit.aes = FALSE) +
    geom_text(data = trait_pca[[2]] |>
                mutate(PC3 = case_when(label == "C_percent" ~ -0.1,
                                       TRUE ~ PC3),
                       PC4 = case_when(label == "Plant_Height_cm_log" ~ 0,
                                       label == "Thickness_mm_log" ~ -0.8,
                                       label == "C_percent" ~ 0.5,
                                       label == "dN15_permil" ~ 0.42,
                                       label == "P_percent" ~ -0.1,
                                       TRUE ~ PC4)),
              aes(x = PC3 + 0.2, y = PC4, label = figure_names, colour = class),
              size = 2.5,
              inherit.aes = FALSE,
              show.legend = FALSE, parse = TRUE) +
    coord_equal(xlim = c(PC3_min+1.5, PC3_max-1), ylim = c(PC4_min+1.5, PC4_max-1.5)) +
    labs(x = "PCA3", y = "PCA4", tag = "(c)") +
    scale_x_continuous(expand = c(.2, 0)) +
    scale_linetype_manual(name = "", values = c("solid", "dashed", "solid", "solid")) +
    scale_colour_manual(name = "", values = c("black", "grey40", "grey70", "cornflowerblue")) +
    theme_minimal() +
    theme(aspect.ratio = 1,
          plot.tag.position = c(0, 1),
          plot.tag = element_text(vjust = 1.5, hjust = -2.85, size = 10))

  plot_layout <- "
  12
  34
  "

  trait_ordination_plot <- pcaC + pcaN + arrow + guide_area() + plot_layout(design = plot_layout, guides = "collect") &
    theme(legend.box = "horizontal")

  return(trait_ordination_plot)

}



# make_full_trait_pca_plot <- function(trait_pca){
#
#   # elevational range
#   range <- range(trait_pca[[1]]$Mean_elevation)
#
#   # both gradients
#   e_B <- eigenvals(trait_pca[[3]])/sum(eigenvals(trait_pca[[3]]))
#
#   # PC1 vs. PC2
#   pca12 <- trait_pca[[1]] %>%
#     ggplot(aes(x = PC1, y = PC2, colour = Mean_elevation, linetype = Gradient, group = GS)) +
#     geom_point(size = 2) +
#     coord_equal(clip = "off", xlim = c(-1.5, 3)) +
#     stat_ellipse(aes(colour = Mean_elevation)) +
#     scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1, name = "Elevation m a.s.l.", limits = c(range[1], range[2])) +
#     scale_linetype_manual(values = c("dashed", "solid"), labels = c("Reference", "Nutrient input")) +
#     labs(x = glue("PCA1 ({round(e_B[1] * 100, 1)}%)"),
#          y = glue("PCA2 ({round(e_B[2] * 100, 1)}%)"),
#          tag = "(a)") +
#     theme_minimal() +
#     theme(aspect.ratio = 1,
#           plot.tag.position = c(0, 0.9),
#           plot.tag = element_text(vjust = -1.5, hjust = -0.5, size = 10))
#
#   # PC1 vs. PC3
#   pca13 <- trait_pca[[1]] %>%
#     ggplot(aes(x = PC1, y = PC3, colour = Mean_elevation, linetype = Gradient, group = GS)) +
#     geom_point(size = 2) +
#     coord_equal(clip = "off", xlim = c(-1.5, 3)) +
#     stat_ellipse(aes(colour = Mean_elevation)) +
#     scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1, name = "Elevation m a.s.l.", limits = c(range[1], range[2])) +
#     scale_linetype_manual(values = c("dashed", "solid"), labels = c("Reference", "Nutrient input")) +
#     labs(x = glue("PCA1 ({round(e_B[1] * 100, 1)}%)"),
#          y = glue("PCA3 ({round(e_B[3] * 100, 1)}%)"),
#          tag = "(c)") +
#     theme_minimal() +
#     theme(aspect.ratio = 1,
#           plot.tag.position = c(0, 0.9),
#           plot.tag = element_text(vjust = -1.5, hjust = -0.5, size = 10))
#
#   arrow12 <- trait_pca[[1]] %>%
#     ggplot(aes(x = PC1, y = PC2)) +
#     geom_segment(data = trait_pca[[2]],
#                  aes(x = 0, y = 0, xend = PC1, yend = PC2, colour = class, linetype = class),
#                  arrow = arrow(length = unit(0.2, "cm")),
#                  inherit.aes = FALSE) +
#     geom_text(data = trait_pca[[2]],
#               aes(x = PC1 * 1.1, y = PC2 * 1.1, label = trait_fancy, colour = class),
#               size = 2.5,
#               inherit.aes = FALSE,
#               show.legend = FALSE) +
#     labs(x = "PC 1", y = "PC 2", tag = "(b)") +
#     scale_x_continuous(expand = c(.2, 0)) +
#     scale_linetype_manual(name = "", values = c("solid", "dashed", "solid")) +
#     scale_colour_manual(name = "", values = c("black", "grey40", "grey70")) +
#     theme_minimal() +
#     theme(aspect.ratio = 1,
#           plot.tag.position = c(0, 1),
#           plot.tag = element_text(vjust = 1.5, hjust = -2.85, size = 10))
#
#
#   arrow13 <- trait_pca[[1]] %>%
#     ggplot(aes(x = PC1, y = PC3)) +
#     geom_segment(data = trait_pca[[2]],
#                  aes(x = 0, y = 0, xend = PC1, yend = PC3, colour = class, linetype = class),
#                  arrow = arrow(length = unit(0.2, "cm")),
#                  inherit.aes = FALSE) +
#     geom_text(data = trait_pca[[2]],
#               aes(x = PC1 * 1.1, y = PC3 * 1.1, label = trait_fancy, colour = class),
#               size = 2.5,
#               inherit.aes = FALSE,
#               show.legend = FALSE) +
#     labs(x = "PC 1", y = "PC 2", tag = "(d)") +
#     scale_x_continuous(expand = c(.2, 0)) +
#     scale_linetype_manual(name = "", values = c("solid", "dashed", "solid")) +
#     scale_colour_manual(name = "", values = c("black", "grey40", "grey70")) +
#     theme_minimal() +
#     theme(aspect.ratio = 1,
#           plot.tag.position = c(0, 1),
#           plot.tag = element_text(vjust = 1.5, hjust = -2.85, size = 10))
#
#
#   trait_ordination_plot <- wrap_plots(pca12, arrow12, pca13, arrow13, ncol = 2) +
#     plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.box="vertical")
#
#   return(trait_ordination_plot)
# }
