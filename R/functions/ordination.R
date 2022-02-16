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
    scale_linetype_manual(values = c(1, 2), , labels = c("Birdcliff", "Reference"))

  return(ordi_plot)
}



## trait ordinations
make_trait_pca <- function(trait_mean){

  # make wide trait table
  cwm_fat <- trait_mean %>%
    select(Gradient:mean) %>%
    pivot_wider(names_from = "trait_trans", values_from = "mean")

  pca_output <- cwm_fat %>%
    select(-(Gradient:PlotID)) %>%
    rda(scale = TRUE)

  pca_sites <- bind_cols(
    cwm_fat %>%
      select(Gradient:PlotID),
    fortify(pca_output, display = "sites")
  )

  pca_traits <- fortify(pca_output, display = "species")

  outputList <- list(pca_sites, pca_traits)

  return(outputList)
}

# Function to make ellipse
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100){
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
