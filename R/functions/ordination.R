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
    mutate(Gradient = recode(Gradient, "B" = "Nutrient input", "C" = "Reference")) %>%
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
                     filter(Score == "sites"))

  sp <- fortify(res) |>
    filter(Score == "species")

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
           Label = capitalize(Label)) |>
    filter(length > 0.5) |>
    select(Label, length, PC1, PC2)

  # make main figure
  # reference
  pca_C <- comm_pca[[1]] |>
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
    filter(Gradient == "Nutrient input") |>
    ggplot(aes(x = PC1, y = PC2, colour = Mean_elevation)) +
    geom_point(size = 2) +
    #annotation_custom(bird, xmin = -3.5, xmax = -4.5, ymin = 1.3, ymax = 2.6) +
    coord_equal(xlim = c(PC1_min, PC1_max),
                ylim = c(PC2_min, PC2_max)) +
    labs(x = glue("PCA1 ({round(e_B[1] * 100, 1)}%)"),
         y = glue("PCA2 ({round(e_B[2] * 100, 1)}%)"),
         title = "Nutrient input",
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
                mutate(Label = capitalize(Label)) |>
                inner_join(important_species) |>
                mutate(PC1 = case_when(Label == "Salix polaris" ~ -1.5,
                                       Label == "Oxyria digyna" ~ 1.2,
                                       Label == "Cerastium arcticum" ~ 0.5,
                                       Label == "Bistorta vivipara" ~ 0.5,
                                       Label == "Dryas octopetala" ~ 0.5,
                                        TRUE ~ PC1),
                       PC2 = case_when(Label == "Salix polaris" ~ 0.2,
                                       Label == "Oxyria digyna" ~ 0.4,
                                       Label == "Cerastium arcticum" ~ 0.5,
                                       Label == "Bistorta vivipara" ~ -1.3,
                                       Label == "Dryas octopetala" ~ -2,
                                       TRUE ~ PC2)),
              aes(x = PC1, y = PC2, label = Label),
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

  pca_output <- cwm_fat %>%
    select(-(Gradient:SoilTemperature)) %>%
    rda(scale = TRUE, center = TRUE)

  pca_sites <- bind_cols(
    cwm_fat %>%
      select(Gradient:SoilTemperature),
    fortify(pca_output, display = "sites")
  )

  pca_traits <- fortify(pca_output, display = "species") %>%
    mutate(trait_trans = Label) %>%
    fancy_trait_name_dictionary()


  # permutation test
  # traits
  raw <- cwm_fat %>% select(-(Gradient:SoilTemperature))
  # meta data
  meta <- cwm_fat %>% select(Gradient:SoilTemperature) %>%
    mutate(Site = factor(Site))

  # adonis test
  if(meta %>% distinct(Gradient) %>% count() == 2){
    adonis_result <- adonis(raw ~ Gradient*Mean_elevation , data = meta, permutations = 999, method = "euclidean")
  } else {
    adonis_result <- adonis(raw ~ Mean_elevation, data = meta, permutations = 999, method = "euclidean")
  }

  outputList <- list(pca_sites, pca_traits, pca_output, adonis_result)

  return(outputList)
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



make_pca_plot <- function(trait_pca, PC1, PC2){

  # adding pics
  # bird <- grid::rasterGrob(png::readPNG("bird.png"), interpolate = TRUE)
  # ref <- grid::rasterGrob(png::readPNG("ref.png"), interpolate = TRUE)
  # add this to figure
  # annotation_custom(bird, xmin = -2.5, xmax = -3.8, ymin = 1.3, ymax = 2.6)

  # elevational range
  range <- range(trait_pca[[1]]$Mean_elevation)

  # PC range
  PC1_min <- min(trait_pca[[1]]$PC1) - 0.5
  PC1_max <- max(trait_pca[[1]]$PC1) + 1
  PC2_min <- min(trait_pca[[1]]$PC2) - 1
  PC2_max <- max(trait_pca[[1]]$PC2) + 0.5

  # both gradients
  e_B <- eigenvals(trait_pca[[3]])/sum(eigenvals(trait_pca[[3]]))

  # Reference
  pcaC <- trait_pca[[1]] %>%
    filter(Gradient == "C") |>
    ggplot(aes(x = PC1, y = PC2, colour = Mean_elevation, group = GS)) +
    geom_point(size = 2) +
    coord_equal(clip = "off", xlim = c(PC1_min, PC1_max), ylim = c(PC2_min, PC2_max)) +
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
    coord_equal(clip = "off", xlim = c(PC1_min, PC1_max), ylim = c(PC2_min, PC2_max)) +
    stat_ellipse(aes(colour = Mean_elevation)) +
    scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1, name = "Elevation m a.s.l.", limits = c(range[1], range[2])) +
    labs(x = glue("PCA1 ({round(e_B[1] * 100, 1)}%)"),
         y = glue("PCA2 ({round(e_B[2] * 100, 1)}%)"),
         tag = "(b) Nutrient input") +
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
                mutate(PC1 = case_when(Label == "N_percent" ~ 1,
                                       Label == "P_percent" ~ 1.3,
                                       TRUE ~ PC1),
                       PC2 = case_when(Label == "N_percent" ~ -0.4,
                                       Label == "P_percent" ~ -0.16,
                                       Label == "Thickness_mm_log" ~ 0.02,
                                       Label == "Leaf_Area_cm2_log" ~ 0.5,
                                       TRUE ~ PC2)),
              aes(x = PC1 * 1.1, y = PC2 * 1.1, label = trait_fancy, colour = class),
              size = 2.5,
              inherit.aes = FALSE,
              show.legend = FALSE) +
    coord_equal(clip = "off", xlim = c(PC1_min+0.5, PC1_max-1), ylim = c(PC2_min+1, PC2_max-2)) +
    labs(x = "PCA1", y = "PCA2", tag = "(c)") +
    scale_x_continuous(expand = c(.2, 0)) +
    scale_linetype_manual(name = "", values = c("solid", "dashed", "solid")) +
    scale_colour_manual(name = "", values = c("black", "grey40", "grey70")) +
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

