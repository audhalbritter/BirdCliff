# Ordination
make_ordination <- function(comm_raw){

  set.seed(32)

  comm_fat <- comm_raw %>%
    select(Gradient, Site, GS, Elevation_m, PlotID, Taxon) %>%
    mutate(presence = 1) %>%
    pivot_wider(names_from = Taxon,
                values_from = presence,
                values_fill = 0)

  comm_fat_spp <- comm_fat %>% select(-(Gradient:PlotID))

  NMDS <-  metaMDS(comm_fat_spp, noshare = TRUE, try = 30)

  # fortify
  fNMDS <- fortify(NMDS) %>%
    filter(Score == "sites") %>%
    bind_cols(comm_fat %>% select(Gradient:PlotID))

  return(fNMDS)
}


### NEEDS A BETTER MODEL!!!
test_ordination <- function(comm_raw){
  set.seed(32)

  comm_fat <- comm_raw %>%
    select(Gradient, Site, GS, Elevation_m, PlotID, Taxon) %>%
    mutate(presence = 1) %>%
    pivot_wider(names_from = Taxon,
                values_from = presence,
                values_fill = 0)

  env <- comm_fat %>% select(Gradient:PlotID)
  attach(env)
  dis <- vegdist(comm_fat_spp)
  dune.ano <- anosim(dis, GS)
  summary(dune.ano)
  plot(dune.ano)

  }




make_ordination_plot <- function(comm_raw){

  set.seed(32)

  comm_fat <- comm_raw %>%
    select(Gradient, Site, GS, Elevation_m, PlotID, Taxon) %>%
    mutate(presence = 1) %>%
    pivot_wider(names_from = Taxon,
                values_from = presence,
                values_fill = 0)

  comm_fat_spp <- comm_fat %>% select(-(Gradient:PlotID))

  NMDS <-  metaMDS(comm_fat_spp, noshare = TRUE, try = 30)

  # extract important data
  fNMDS <- fortify(NMDS) %>%
    filter(Score == "sites") %>%
    bind_cols(comm_fat %>% select(Gradient:PlotID))

  env <- comm_fat %>% select(Gradient:Elevation_m)

  # extract NMDS scores
  sites <- as.data.frame(scores(NMDS, display = "sites"))
  sites <- cbind(sites, GS = env$GS)
  sites$GS <- factor(sites$GS)

  #data for ellipse, in this case using the management factor
  df_ell <- data.frame() #sets up a data frame before running the function.
  for(g in levels(sites$GS)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(sites [sites$GS==g,],
                                                     veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) ,GS=g))
  }

  df_ell <- df_ell %>%
    left_join(env %>% distinct(GS, Gradient, Site, Elevation_m), by = "GS")

  g0 <- fNMDS %>%
    ggplot(aes(x = NMDS1, y = NMDS2, group = Gradient, colour = Elevation_m, shape = Gradient)) +
    geom_point(size = 2) +
    coord_equal() +
    scale_colour_viridis_c(name = "Elevation in m a.s.l.",
                           end = 0.8,
                           option = "inferno",
                           direction = -1) +
    scale_fill_viridis_c(name = "Elevation", end = 0.8, option = "inferno", direction = -1) +
    scale_shape_manual(values = c(16, 2), , labels = c("Birdcliff", "Reference")) +
    labs(x = "NMDS axis 1", y = "NMDS axis 2") +
    theme_minimal()

  ordi_plot <- g0 +
    geom_path(data = df_ell, aes(x = NMDS1, y = NMDS2, group = GS, colour = Elevation_m, shape = Gradient, linetype = Gradient)) +
    scale_linetype_manual(values = c(1, 2), , labels = c("Bird cliff", "Reference"))

  return(ordi_plot)
}



## trait ordinations
make_trait_pca <- function(trait_mean){

  # make wide trait table
  cwm_fat <- trait_mean %>%
    mutate(GS = paste0(Gradient, Site)) %>%
    select(Gradient:mean, Elevation_m, GS) %>%
    pivot_wider(names_from = "trait_trans", values_from = "mean")

  pca_output <- cwm_fat %>%
    select(-(Gradient:GS)) %>%
    rda(scale = TRUE)

  pca_sites <- bind_cols(
    cwm_fat %>%
      select(Gradient:GS),
    fortify(pca_output, display = "sites")
  )

  pca_traits <- fortify(pca_output, display = "species") %>%
    mutate(trait_trans = Label) %>%
    fancy_trait_name_dictionary()


  # permutation test
  # traits
  raw <- cwm_fat %>% select(-(Gradient:GS))
  # meta data
  meta <- cwm_fat %>% select(Gradient:GS) %>%
    mutate(Site = factor(Site))

  # test
  adonis_result <- adonis2(raw ~ Site, data = meta, permutations = 999, method = "euclidean")

  outputList <- list(pca_sites, pca_traits, pca_output, adonis_result)

  return(outputList)
}





make_trait_pca_plot <- function(trait_pca_B, trait_pca_C){

  plot_B <- trait_pca_B[[1]] %>%
    ggplot(aes(x = PC1, y = PC2, colour = Elevation_m, group = Site)) +
    geom_point(size = 2) +
    coord_equal() +
    stat_ellipse() +
    scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1, name = "Elevation m a.s.l.") +
    scale_shape_manual(values = c(16)) +
    labs(x = "PC 1 (47.2%)", y = "PC 2 (17.2%)", title = "Bird cliff") +
    theme_minimal() +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5),
          aspect.ratio = 1)

  arrow_B <- trait_pca_B[[1]] %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_segment(data = trait_pca_B[[2]],
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "cm")),
                 colour = "grey50",
                 inherit.aes = FALSE) +
    geom_text(data = trait_pca_B[[2]],
              aes(x = PC1 * 1.1,y = PC2 * 1.1, label = trait_fancy),
              size = 3,
              inherit.aes = FALSE, colour = "black") +
    labs(x = "PC 1", y = "PC 2") +
    scale_x_continuous(expand = c(.2, 0)) +
    theme_minimal() +
    theme(aspect.ratio = 1)

  plot_C <- trait_pca_C[[1]] %>%
    ggplot(aes(x = PC1, y = PC2, colour = Elevation_m, group = Site)) +
    geom_point(size = 2, shape = 2) +
    coord_equal() +
    stat_ellipse(linetype = 2) +
    scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1, name = "Elevation m a.s.l.") +
    labs(x = "PC 1 (35.5%)", y = "PC 2 (21.6%)", title = "Reference") +
    theme_minimal() +
    theme(aspect.ratio = 1)

  arrow_C <- trait_pca_C[[1]] %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_segment(data = trait_pca_C[[2]] %>%
                   mutate(PC1 = PC1*-1,
                          PC2 = PC2*-1),
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "cm")),
                 colour = "grey50",
                 inherit.aes = FALSE) +
    geom_text(data = trait_pca_C[[2]] %>%
                mutate(PC1 = PC1*-1,
                       PC2 = PC2*-1),
              aes(x = PC1 * 1.1, y = PC2 * 1.1, label = trait_fancy),
              size = 3,
              inherit.aes = FALSE, colour = "black") +
    labs(x = "PC 1", y = "PC 2") +
    scale_x_continuous(expand = c(.2, 0)) +
    theme_minimal() +
    theme(aspect.ratio = 1)

  layout <- "
  AAB
  CCD
"
  trait_ordination_plot <- wrap_plots(plot_B, arrow_B, plot_C, arrow_C) + plot_layout(design = layout)

}



# Function to make ellipse
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100){
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
