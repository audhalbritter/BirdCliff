# NMDS

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
