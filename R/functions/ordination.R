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
    scale_linetype_manual(values = c(1, 2), , labels = c("Bird cliff", "Reference")) +
    scale_shape_manual(values = c(16, 2), , labels = c("Bird cliff", "Reference")) +
    labs(x = "NMDS axis 1", y = "NMDS axis 2") +
    theme_minimal()

  # ordi_plot12 <- g12 +
  #   geom_path(data = df_ell_12, aes(x = NMDS1, y = NMDS2, group = GS, colour = Elevation_m)) +
  #   scale_linetype_manual(values = c(1, 2), , labels = c("Bird cliff", "Reference"))

  return(NMDS_plot)
}



## TRAITS (PCA)
make_trait_pca <- function(trait_mean){

  # make wide trait table
  cwm_fat <- trait_mean %>%
    group_by(Site) %>%
    mutate(Mean_elevation = mean(Elevation_m),
           GS = paste0(Gradient, Site)) %>%
    select(Gradient:mean, Mean_elevation, GS) %>%
    pivot_wider(names_from = "trait_trans", values_from = "mean") %>%
    ungroup()

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

  # adonis test
  if(meta %>% distinct(Gradient) %>% count() == 2){
    adonis_result <- adonis2(raw ~ Site*Gradient, data = meta, permutations = 999, method = "euclidean")
  } else {
    adonis_result <- adonis2(raw ~ Site, data = meta, permutations = 999, method = "euclidean")
  }


  outputList <- list(pca_sites, pca_traits, pca_output, adonis_result)

  return(outputList)
}





make_trait_pca_plot <- function(trait_pca_B, trait_pca_C){

  # elevational range
  range <- range(trait_pca_C[[1]]$Mean_elevation)

  # prop explained
  e_B <- eigenvals(trait_pca_B[[3]])/sum(eigenvals(trait_pca_B[[3]]))

  plot_B <- trait_pca_B[[1]] %>%
    ggplot(aes(x = PC1, y = PC2, colour = Mean_elevation, group = Site)) +
    geom_point(size = 2) +
    coord_equal() +
    stat_ellipse(aes(colour = Mean_elevation)) +
    scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1, name = "Elevation m a.s.l.", limits = c(range[1], range[2])) +
    labs(x = glue("PCA1 ({round(e_B[1] * 100, 1)}%)"),
         y = glue("PCA1 ({round(e_B[2] * 100, 1)}%)"),
         tag = "(a) Bird cliff") +
    theme_minimal() +
    theme(aspect.ratio = 1,
          plot.tag.position = c(0, 0.9),
          plot.tag = element_text(vjust = -1.5, hjust = -0.5, size = 10))

  arrow_B <- trait_pca_B[[1]] %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_segment(data = trait_pca_B[[2]],
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "cm")),
                 colour = "grey50",
                 inherit.aes = FALSE) +
    geom_text(data = trait_pca_B[[2]],
              aes(x = PC1 * 1.1,y = PC2 * 1.1, label = trait_fancy),
              size = 2.5,
              inherit.aes = FALSE, colour = "black") +
    labs(x = "PC 1", y = "PC 2", tag = "(b)") +
    scale_x_continuous(expand = c(.2, 0)) +
    theme_minimal() +
    theme(aspect.ratio = 1,
          plot.tag.position = c(0, 1),
          plot.tag = element_text(vjust = 1.5, hjust = -2.85, size = 10))


  # prop explained
  e_C <- eigenvals(trait_pca_C[[3]])/sum(eigenvals(trait_pca_C[[3]]))

  plot_C <- trait_pca_C[[1]] %>%
    ggplot(aes(x = PC1, y = PC2, colour = Mean_elevation, group = Site)) +
    geom_point(size = 2, shape = 2) +
    coord_equal() +
    stat_ellipse(aes(colour = Mean_elevation), linetype = 2) +
    scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1, name = "Elevation m a.s.l.") +
    labs(x = glue("PCA1 ({round(e_C[1] * 100, 1)}%)"),
         y = glue("PCA1 ({round(e_C[2] * 100, 1)}%)"),
         tag = "(c) Reference") +
    theme_minimal() +
    theme(aspect.ratio = 1,
          plot.tag.position = c(0, 0.9),
          plot.tag = element_text(vjust = -1, hjust = -0.5, size = 10))

  arrow_C <- trait_pca_C[[1]] %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_segment(data = trait_pca_C[[2]] %>%
                   mutate(PC1 = PC1, #*-1,
                          PC2 = PC2), #*-1),
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "cm")),
                 colour = "grey50",
                 inherit.aes = FALSE) +
    geom_text(data = trait_pca_C[[2]] %>%
                mutate(PC1 = PC1, #*-1,
                       PC2 = PC2), #*-1),
              aes(x = PC1 * 1.1, y = PC2 * 1.1, label = trait_fancy),
              size = 2.5,
              inherit.aes = FALSE, colour = "black") +
    labs(x = "PC 1", y = "PC 2", tag = "(d)") +
    scale_x_continuous(expand = c(.2, 0)) +
    theme_minimal() +
    theme(aspect.ratio = 1,
          plot.tag.position = c(0, 1),
          plot.tag = element_text(vjust = 1.5, hjust = -2.85, size = 10))



  trait_ordination_plot <- wrap_plots(plot_B, arrow_B, plot_C, arrow_C) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

  return(trait_ordination_plot)

}



# Function to make ellipse
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100){
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
