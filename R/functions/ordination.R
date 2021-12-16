# Ordination
make_ordination <- function(comm_raw){

  set.seed(32)

  comm_fat <- comm_raw %>%
    select(Gradient, Site, PlotID, Taxon) %>%
    mutate(presence = 1) %>%
    pivot_wider(names_from = Taxon,
                values_from = presence,
                values_fill = 0)

  comm_fat_spp <- comm_fat %>% select(-(Gradient:PlotID))

  NMDS <-  metaMDS(comm_fat_spp, noshare = TRUE, try = 30)

  fNMDS <- fortify(NMDS) %>%
    filter(Score == "sites") %>%
    bind_cols(comm_fat %>% select(Gradient:PlotID))

  return(fNMDS)
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
