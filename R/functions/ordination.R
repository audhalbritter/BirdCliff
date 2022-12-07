# ORDINATIONS

# COMMUNITY
check_dimensions_NMDS <- function(comm_raw){

  set.seed(32)

  comm_fat <- comm_raw %>%
    group_by(Site) %>%
    mutate(Mean_elevation = mean(Elevation_m)) %>%
    ungroup() %>%
    select(Gradient, Site, GS, Mean_elevation, PlotID, Taxon) %>%
    mutate(presence = 1) %>%
    pivot_wider(names_from = Taxon,
                values_from = presence,
                values_fill = 0)

  comm_fat_spp <- comm_fat %>% select(-(Gradient:PlotID))

  NMDS_1 <-  metaMDS(comm_fat_spp, noshare = TRUE, try = 30, k = 1)
  NMDS_2 <-  metaMDS(comm_fat_spp, noshare = TRUE, try = 30, k = 2)
  NMDS_3 <-  metaMDS(comm_fat_spp, noshare = TRUE, try = 30, k = 3)
  NMDS_4 <-  metaMDS(comm_fat_spp, noshare = TRUE, try = 30, k = 4)
  NMDS_5 <-  metaMDS(comm_fat_spp, noshare = TRUE, try = 30, k = 5)
  NMDS_6 <-  metaMDS(comm_fat_spp, noshare = TRUE, try = 30, k = 6)

  stress_plot <- tibble(
    stress = c(NMDS_1$stress, NMDS_2$stress, NMDS_3$stress, NMDS_4$stress, NMDS_5$stress, NMDS_6$stress),
    dimensions = c(1:6)) %>%
    ggplot(aes(x = dimensions, y = stress)) +
    geom_point()

  return(stress_plot)
}

make_ordination <- function(comm_raw){

  set.seed(32)

  comm_fat <- comm_raw %>%
    group_by(Site) %>%
    mutate(Mean_elevation = mean(Elevation_m)) %>%
    ungroup() %>%
    select(Gradient, Site, GS, Mean_elevation, PlotID, Taxon) %>%
    mutate(presence = 1) %>%
    pivot_wider(names_from = Taxon,
                values_from = presence,
                values_fill = 0)

  comm_fat_spp <- comm_fat %>% select(-(Gradient:PlotID))

  NMDS <-  metaMDS(comm_fat_spp, noshare = TRUE, try = 30, k = 2)

  stress <- NMDS$stress

  #stressplot(NMDS)

  # fortify
  fNMDS <- fortify(NMDS) %>%
    filter(Score == "sites") %>%
    bind_cols(comm_fat %>% select(Gradient:PlotID))

  return(list(NMDS, fNMDS, stress))
}


### test model using adonis
test_ordination <- function(comm_raw){


  comm_fat <- comm_raw %>%
    group_by(Site) %>%
    mutate(Mean_elevation = mean(Elevation_m)) %>%
    ungroup() %>%
    select(Gradient, Site, GS, Mean_elevation, PlotID, Taxon) %>%
    mutate(presence = 1) %>%
    pivot_wider(names_from = Taxon,
                values_from = presence,
                values_fill = 0)

  comm_fat_spp <- comm_fat %>% select(-(Gradient:PlotID))
  meta <- comm_fat %>% select(Gradient:PlotID)

  adon.results <- adonis(comm_fat_spp ~ Gradient*Mean_elevation, data = meta, method = "bray", perm = 999)

  return(adon.results)

  }



# figure
make_ordination_plot <- function(comm_raw, NMDS, fNMDS){

  # env <- comm_raw %>%
  #   group_by(Site) %>%
  #   mutate(Mean_elevation = mean(Elevation_m)) %>%
  #   ungroup() %>%
  #   select(Gradient, Site, GS, Mean_elevation, PlotID, Taxon) %>%
  #   mutate(presence = 1) %>%
  #   pivot_wider(names_from = Taxon,
  #               values_from = presence,
  #               values_fill = 0) %>%
  #   select(Gradient:PlotID)
  #
  # # extract NMDS scores
  # sites <- as.data.frame(scores(NMDS, display = "sites"))
  # sites <- cbind(sites, GS = env$GS)
  # sites$GS <- factor(sites$GS)
  #
  # ### NMDS 1 and 2
  # #data for ellipse, in this case using the management factor
  # df_ell_12 <- data.frame() #sets up a data frame before running the function.
  # for(g in levels(sites$GS)){
  #   df_ell_12 <- rbind(df_ell_12, cbind(as.data.frame(with(sites [sites$GS==g,],
  #                                                    veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) ,GS=g))
  # }
  #
  # df_ell_12 <- df_ell_12 %>%
  #   left_join(env %>% distinct(GS, Gradient, Site, Elevation_m), by = "GS")

  NMDS_plot <- fNMDS %>%
    ggplot(aes(x = NMDS1, y = NMDS2, group = GS, colour = Mean_elevation, shape = Gradient, linetype = Gradient)) +
    geom_point(size = 2) +
    stat_ellipse() +
    coord_equal() +
    scale_colour_viridis_c(name = "Elevation in m a.s.l.",
                           end = 0.8,
                           option = "inferno",
                           direction = -1) +
    scale_linetype_manual(name = "", values = c(1, 2), , labels = c("Nutrient input", "Reference")) +
    scale_shape_manual(name = "", values = c(16, 2), , labels = c("Nutrient input", "Reference")) +
    labs(x = "NMDS axis 1", y = "NMDS axis 2") +
    theme_minimal()

  # ordi_plot12 <- g12 +
  #   geom_path(data = df_ell_12, aes(x = NMDS1, y = NMDS2, group = GS, colour = Elevation_m)) +
  #   scale_linetype_manual(values = c(1, 2), , labels = c("Bird cliff", "Reference"))

  return(NMDS_plot)
}

## COMMUNITY (PCA)

make_community_pca <- function(comm_raw){

  set.seed(32)

  comm_wide <- comm_raw %>%
    group_by(Site) %>%
    mutate(Mean_elevation = mean(Elevation_m)) %>%
    ungroup() %>%
    select(Gradient, Site, GS, Mean_elevation, PlotID, Taxon) %>%
    mutate(presence = 1,
           Gradient = recode(Gradient, "B" = "Nutrient input", "C" = "Reference")) %>%
    pivot_wider(names_from = Taxon,
                values_from = presence,
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
  adonis <- adonis(comm_sp ~ Mean_elevation , data = comm_info, permutations = 999, method = "euclidean")

  return(list(out, sp, res, adonis))

}


make_sp_pca_figure <- function(comm_pca_B, comm_pca_C){

  bird <- grid::rasterGrob(png::readPNG("bird.png"), interpolate = TRUE)
  ref <- grid::rasterGrob(png::readPNG("ref.png"), interpolate = TRUE)

  # elevational range
  range <- range(comm_pca_C[[1]]$Mean_elevation)

  # Nutrient gradient
  e_B <- eigenvals(comm_pca_B[[3]])/sum(eigenvals(comm_pca_B[[3]]))

  important_species_B <- comm_pca_B[[2]] |>
    mutate(length = sqrt(PC1^2 + PC2^2),
           Label = capitalize(Label)) |>
    filter(length > 0.5) |>
    select(Label, length)

  # make main figure
  pca_B <- comm_pca_B[[1]] |>
    ggplot(aes(x = PC1, y = PC2, colour = Mean_elevation)) +
    geom_point(size = 2) +
    coord_equal(clip = "off", xlim = c(-1.5, 3)) +
    labs(x = glue("PCA1 ({round(e_B[1] * 100, 1)}%)"),
         y = glue("PCA2 ({round(e_B[2] * 100, 1)}%)"),
         tag = "(a)") +
    annotate("text", x = -3.2, y = -1, angle = 90, size = 5, label = "Nutrient input") +
    annotation_custom(bird, xmin = -2.5, xmax = -3.8, ymin = 1.3, ymax = 2.6) +
    stat_ellipse(aes(colour = Mean_elevation, group = GS)) +
    scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1, name = "Elevation m a.s.l.", limits = c(range[1], range[2])) +
    scale_shape_manual(values = c(16, 1)) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    theme_bw() +
    theme(text = element_text(size = 13),
          legend.box="vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  # species names
  species_B <- comm_pca_B[[1]] |>
    ggplot(aes(x = PC1, y = PC2)) +
    coord_equal() +
    geom_segment(data = comm_pca_B[[2]] |>
                   mutate(length = sqrt(PC1^2 + PC2^2)),
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2,"cm")),
                 alpha = 0.75,
                 color = 'grey70') +
    geom_text(data = comm_pca_B[[2]] |>
                mutate(Label = capitalize(Label)) |>
                inner_join(important_species_B, by = "Label") |>
                mutate(PC1 = case_when(Label == "Salix polaris" ~ 0.43,
                                       Label == "Saxifraga cespitosa" ~ -0.38,
                                       Label == "Saxifraga cernua" ~ -0.45,
                                       Label == "Cochleria groenlandica" ~ -0.35,
                                       Label == "Oxyria digyna" ~ -0.45,
                                       TRUE ~ PC1),
                       PC2 = case_when(Label == "Salix polaris" ~ 0.12,
                                       Label == "Saxifraga cespitosa" ~ 0.05,
                                       Label == "Saxifraga cernua" ~ 0,
                                       Label == "Cochleria groenlandica" ~ -0.1,
                                       Label == "Oxyria digyna" ~ -0.18,
                                       TRUE ~ PC2)),
              aes(x = PC1*1.1, y = PC2*1.1, label = Label),
              size = 2, col = 'black') +
    labs(x = "PC 1", y = "PC 2", tag = "(b)") +
    theme_bw() +
    theme(text = element_text(size = 13),
          legend.box="vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  # Reference gradient
  e_C <- eigenvals(comm_pca_C[[3]])/sum(eigenvals(comm_pca_C[[3]]))

  important_species_C <- comm_pca_C[[2]] |>
    mutate(length = sqrt(PC1^2 + PC2^2),
           Label = capitalize(Label)) |>
    filter(length > 0.5) |>
    select(Label, length)

  # make main figure
  pca_C <- comm_pca_C[[1]] |>
    ggplot(aes(x = PC1*-1, y = PC2, colour = Mean_elevation)) +
    geom_point(size = 2) +
    coord_equal(clip = "off", xlim = c(-1.5, 3)) +
    annotate("text", x = -3.2, y = -1, angle = 90, size = 5, label = "Reference") +
    annotation_custom(ref, xmin = -2.5, xmax = -3.8, ymin = 1.3, ymax = 2.6) +
    labs(x = glue("PCA1 ({round(e_C[1] * 100, 1)}%)"),
         y = glue("PCA2 ({round(e_C[2] * 100, 1)}%)"),
         tag = "(c)") +
    stat_ellipse(aes(colour = Mean_elevation, group = GS)) +
    scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1, name = "Elevation m a.s.l.", limits = c(range[1], range[2])) +
    scale_shape_manual(values = c(16, 1)) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    theme_bw() +
    theme(text = element_text(size = 13),
          legend.box="vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  # species names
  species_C <- comm_pca_C[[1]] |>
    ggplot(aes(x = PC1*-1, y = PC2)) +
    coord_equal(xlim = c(-0.6, 0.4), ylim = c(-0.5, 0.75)) +
    geom_segment(data = comm_pca_C[[2]] |>
                   mutate(length = sqrt(PC1^2 + PC2^2)),
                 aes(x = 0, y = 0, xend = PC1*-1, yend = PC2),
                 arrow = arrow(length = unit(0.2,"cm")),
                 alpha = 0.75, color = 'grey70') +
    geom_text(data = comm_pca_C[[2]] |>
                mutate(Label = capitalize(Label)) |>
                inner_join(important_species_C, by = "Label") |>
                mutate(PC1 = case_when(Label == "Cerastium arcticum" ~ 0.4,
                                       Label == "Saxifraga cespitosa" ~ 0.42,
                                       Label == "Saxifraga cernua" ~ 0.36,
                                       Label == "Draba alpina" ~ 0.49,
                                       TRUE ~ PC1),
                       PC2 = case_when(Label == "Cerastium arcticum" ~ 0.13,
                                       Label == "Saxifraga cespitosa" ~ -0.16,
                                       Label == "Saxifraga cernua" ~ -0.08,
                                       Label == "Draba alpina" ~ 0.02,
                                       TRUE ~ PC2)),
              aes(x = PC1*-1.1, y = PC2*1.1, label = Label),
              size = 2, col = 'black') +
    labs(x = "PC 1", y = "PC 2", tag = "(d)") +
    theme_bw() +
    theme(text = element_text(size = 13),
          legend.box="vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  wrap_plots(pca_B, species_B, pca_C, species_C) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.box="vertical")

}


# comm_pca_B[[1]] |>
#   ggplot(aes(x = PC1, y = PC2, colour = Mean_elevation, shape = Gradient, linetype = Gradient, group = GS)) +
#   geom_point(size = 2) +
#   coord_equal() +
#   labs(x = glue("PCA1 ({round(e_B[1] * 100, 1)}%)"),
#        y = glue("PCA2 ({round(e_B[2] * 100, 1)}%)")) +
#   stat_ellipse(aes(colour = Mean_elevation)) +
#   scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1, name = "Elevation m a.s.l.", limits = c(range[1], range[2])) +
#   scale_shape_manual(values = c(16, 1)) +
#   scale_linetype_manual(values = c("solid", "dashed")) +
#   facet_wrap(~Gradient) +
#   theme_bw() +
#   theme(text = element_text(size = 13),
#         legend.box="vertical",
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())

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
    rda(scale = TRUE)

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


make_trait_pca_plot <- function(trait_pca_B, trait_pca_C){

  bird <- grid::rasterGrob(png::readPNG("bird.png"), interpolate = TRUE)
  ref <- grid::rasterGrob(png::readPNG("ref.png"), interpolate = TRUE)


  # elevational range
  range <- range(trait_pca_C[[1]]$Mean_elevation)

  # prop explained
  e_B <- eigenvals(trait_pca_B[[3]])/sum(eigenvals(trait_pca_B[[3]]))

  plot_B <- trait_pca_B[[1]] %>%
    ggplot(aes(x = PC1, y = PC2, colour = Mean_elevation, group = Site)) +
    geom_point(size = 2) +
    annotate("text", x = -3.2, y = -1, angle = 90, size = 5, label = "Nutrient input") +
    annotation_custom(bird, xmin = -2.5, xmax = -3.8, ymin = 1.3, ymax = 2.6) +
    coord_equal(clip = "off", xlim = c(-1.5, 3)) +
    stat_ellipse(aes(colour = Mean_elevation)) +
    scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1, name = "Elevation m a.s.l.", limits = c(range[1], range[2])) +
    labs(x = glue("PCA1 ({round(e_B[1] * 100, 1)}%)"),
         y = glue("PCA1 ({round(e_B[2] * 100, 1)}%)"),
         tag = "(a)") +
    theme_minimal() +
    theme(aspect.ratio = 1,
          plot.tag.position = c(0, 0.9),
          plot.tag = element_text(vjust = -1.5, hjust = -0.5, size = 10))

  arrow_B <- trait_pca_B[[1]] %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_segment(data = trait_pca_B[[2]],
                 aes(x = 0, y = 0, xend = PC1, yend = PC2, colour = class, linetype = class),
                 arrow = arrow(length = unit(0.2, "cm")),
                 inherit.aes = FALSE) +
    geom_text(data = trait_pca_B[[2]],
              aes(x = PC1 * 1.1, y = PC2 * 1.1, label = trait_fancy, colour = class),
              size = 2.5,
              inherit.aes = FALSE,
              show.legend = FALSE) +
    labs(x = "PC 1", y = "PC 2", tag = "(b)") +
    scale_x_continuous(expand = c(.2, 0)) +
    scale_linetype_manual(name = "", values = c("solid", "dashed", "solid")) +
    scale_colour_manual(name = "", values = c("black", "grey40", "grey40")) +
    theme_minimal() +
    theme(aspect.ratio = 1,
          plot.tag.position = c(0, 1),
          plot.tag = element_text(vjust = 1.5, hjust = -2.85, size = 10))


  # prop explained
  e_C <- eigenvals(trait_pca_C[[3]])/sum(eigenvals(trait_pca_C[[3]]))

  plot_C <- trait_pca_C[[1]] %>%
    ggplot(aes(x = PC1, y = PC2, colour = Mean_elevation, group = Site)) +
    geom_point(size = 2) +
    annotate("text", x = -4.5, y = -1, angle = 90, size = 5, label = "Reference") +
    annotation_custom(ref, xmin = -3.5, xmax = -5.3, ymin = 1.5, ymax = 2.8) +
    coord_equal(clip = "off", xlim = c(-2.5, 3)) +
    stat_ellipse(aes(colour = Mean_elevation)) +
    scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1, name = "Elevation m a.s.l.") +
    labs(x = glue("PCA1 ({round(e_C[1] * 100, 1)}%)"),
         y = glue("PCA1 ({round(e_C[2] * 100, 1)}%)"),
         tag = "(c)") +
    theme_minimal() +
    theme(aspect.ratio = 1,
          plot.tag.position = c(0, 0.9),
          plot.tag = element_text(vjust = -1, hjust = -0.5, size = 10))

  arrow_C <- trait_pca_C[[1]] %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_segment(data = trait_pca_C[[2]] %>%
                   mutate(PC1 = PC1, #*-1,
                          PC2 = PC2 *-1),
                 aes(x = 0, y = 0, xend = PC1, yend = PC2, colour = class, linetype = class),
                 arrow = arrow(length = unit(0.2, "cm")),
                 inherit.aes = FALSE) +
    geom_text(data = trait_pca_C[[2]] %>%
                mutate(PC1 = PC1, #*-1,
                       PC2 = PC2 *-1),
              aes(x = PC1 * 1.1, y = PC2 * 1.1, label = trait_fancy, colour = class),
              size = 2.5,
              inherit.aes = FALSE,
              show.legend = FALSE) +
    labs(x = "PC 1", y = "PC 2", tag = "(d)") +
    scale_x_continuous(expand = c(.2, 0)) +
    scale_linetype_manual(name = "", values = c("solid", "dashed", "solid")) +
    scale_colour_manual(name = "", values = c("black", "grey40", "grey40")) +
    theme_minimal() +
    theme(aspect.ratio = 1,
          plot.tag.position = c(0, 1),
          plot.tag = element_text(vjust = 1.5, hjust = -2.85, size = 10))

  trait_ordination_plot <- wrap_plots(plot_B, arrow_B, plot_C, arrow_C) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.box="vertical")

  return(trait_ordination_plot)

}



# Function to make ellipse
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100){
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
