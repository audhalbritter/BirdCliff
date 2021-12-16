#### TRAIT BOOTSTRAPPING ####


make_bootstrapping <- function(comm_raw, trait_raw){

  #prepare community data
  comm <- comm_raw

  #prepare trait data
  trait <- trait_raw %>%
    # remove bryophyte traits
    filter(! Trait %in% c("Shoot_Length_cm", "Shoot_Length_Green_cm")) %>%
    select(Gradient, Site, PlotID, Taxon, Trait, Value) %>%
    filter(Trait != "Wet_Mass_g")

  #prepare trait data without intraspecific variation
  # trait.null <- trait_raw %>%
  #   # remove bryophyte traits
  #   filter(! Trait %in% c("Shoot_Length_cm", "Shoot_Length_Green_cm")) %>%
  #   select(Treatment, Site, PlotID, Taxon, Trait, Value) %>%
  #   filter(Trait != "Wet_Mass_g") %>%
  #   group_by(Taxon, Trait) %>%
  #   summarize(Value = mean(as.numeric(Value), na.rm = T)) %>%
  #   right_join(trait, by = c("Taxon", "Trait")) %>%
  #   select(-Value.y, "Value" = Value.x)

  #set seed for bootstrapping repeatability
  set.seed(2525)
  trait_imp <- trait_impute(comm = comm,
                            traits = trait,
                            scale_hierarchy = c("Gradient", "Site", "PlotID"),
                            global = F,
                            taxon_col = "Taxon",
                            trait_col = "Trait",
                            value_col = "Value",
                            abundance_col = "Cover",
                            min_n_in_sample = 2
  )

  # trait_imp_null <- trait_impute(comm = comm,
  #                                traits = trait.null,
  #                                scale_hierarchy = c("Site", "Site_trt", "PlotID"),
  #                                global = F,
  #                                taxon_col = "Taxon",
  #                                trait_col = "Trait",
  #                                value_col = "Value",
  #                                abundance_col = "Cover",
  #                                min_n_in_sample = 2
  # )

  #check trait coverage
  trait_imp %>%
    #filter(Trait == "C_percent") %>%
    autoplot(.) +
    theme(axis.text.x = element_text(angle = 90))

  # trait_imp_null %>%
  #   autoplot(.) +
  #   theme(axis.text.x = element_text(angle = 90))

  #do the bootstrapping
  CWM <- trait_np_bootstrap(trait_imp, nrep = 100, sample_size = 200)
  #CWM_notiv <- trait_np_bootstrap(trait_imp_null, nrep = 100, sample_size = 200)

  CWM_mean <- trait_summarise_boot_moments(CWM) %>%
    select(Site:mean)

  # CWM_notiv_mean <- trait_summarise_boot_moments(CWM_notiv) %>%
  #   select(Site:mean) %>%
  #   rename("mean_noitv" = "mean")

  # traitMean <- CWM_mean %>%
  #   left_join(CWM_notiv_mean) %>%
  #   select(-n)

  #prepare bootstrapped trait data for analyses
  traitMean <- CWM_mean %>%
    ungroup() %>%
    select(-n) %>%
    mutate(Trait = factor(Trait, levels = c("Plant_Height_cm", "Dry_Mass_g", "Leaf_Area_cm2", "Leaf_Thickness_mm", "SLA_cm2_g", "LDMC", "C_percent", "N_percent", "CN_ratio", "P_percent", "NP_ratio", "dC13_permil", "dN15_permil")))

  return(traitMean)

}
