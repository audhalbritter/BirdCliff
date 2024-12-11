## COMMUNITY ORDINATION (PCA)

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
  adonis <- adonis2(comm_sp ~ Gradient * Mean_elevation , data = comm_info, permutations = 999, method = "euclidean")

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

  # calculte centroids
  centroids <- comm_pca[[1]] |>
    left_join(comm_pca[[1]] |>
                group_by(Gradient, Site) |>
                summarise(centroid1 = mean(PC1),
                          centroid2 = mean(PC2)),
              by = c("Gradient", "Site"))

  e_B <- eigenvals(comm_pca[[3]])/sum(eigenvals(comm_pca[[3]]))

  important_species <- comm_pca[[2]] |>
    mutate(length = sqrt(PC1^2 + PC2^2),
           label = capitalize(label)) |>
    filter(length > 0.5) |>
    select(label, length, PC1, PC2)

  # make main figure
  # reference
  pca <- comm_pca[[1]] |>
    ggplot(aes(x = PC1, y = PC2, colour = Mean_elevation, shape = Gradient)) +
    geom_point(size = 2) +
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
    labs(x = glue("PCA1 ({round(e_B[1] * 100, 1)}%)"),
         y = glue("PCA2 ({round(e_B[2] * 100, 1)}%)"),
         tag = "(a)") +
    stat_ellipse(aes(linetype = Gradient), alpha = 0.2) +
    scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1, name = "Elevation m a.s.l.", limits = c(range[1], range[2])) +
    scale_shape_manual(values = c(1, 16)) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    theme_bw() +
    theme(aspect.ratio = 1,
          legend.position = "top",
          plot.tag.position = c(0, 0.9),
          plot.tag = element_text(vjust = -1.5, hjust = -0.5, size = 10))

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
    labs(x = "PC 1", y = "PC 2", tag = "(b)") +
    theme_bw() +
    theme(aspect.ratio = 1,
          plot.tag.position = c(0, 0.9),
          plot.tag = element_text(vjust = -1.5, hjust = -0.5, size = 10))

  plot_layout <- "
  12
  12
  33
  "

  pca + species + guide_area() + plot_layout(design = plot_layout, guides = "collect") &
    theme(legend.box = "vertical")

}
