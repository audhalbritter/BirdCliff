#### TRAIT BOOTSTRAPPING ####


make_bootstrapping <- function(comm_raw, traits_raw){

  #prepare community data
  comm <- comm_raw

  #prepare trait data
  trait <- traits_raw %>%
    select(Gradient, Site, PlotID, Taxon, trait_trans, value_trans)

  #prepare trait data without intraspecific variation
  trait.null <- traits_raw %>%
    select(Gradient, Site, PlotID, Taxon, trait_trans, value_trans) %>%
    # do we want to have fixed traits per gradient or across all the gradients???
    group_by(Taxon, trait_trans) %>%
    summarise(value_trans = mean(as.numeric(value_trans), na.rm = TRUE)) %>%
    right_join(trait, by = c("Taxon", "trait_trans")) %>%
    select(-value_trans.y, "value_trans" = value_trans.x)

  #set seed for bootstrapping repeatability
  set.seed(2525)
  trait_imp <- trait_impute(comm = comm,
                            traits = trait,
                            scale_hierarchy = c("Gradient", "Site", "PlotID"),
                            global = F,
                            taxon_col = "Taxon",
                            trait_col = "trait_trans",
                            value_col = "value_trans",
                            abundance_col = "Cover",
                            min_n_in_sample = 2
  )

  trait_imp_null <- trait_impute(comm = comm,
                                 traits = trait.null,
                                 scale_hierarchy = c("Gradient", "Site", "PlotID"),
                                 global = F,
                                 taxon_col = "Taxon",
                                 trait_col = "trait_trans",
                                 value_col = "value_trans",
                                 abundance_col = "Cover",
                                 min_n_in_sample = 2
  )

  #check trait coverage
  trait_imp %>%
    autoplot(.) +
    theme(axis.text.x = element_text(angle = 90))

  trait_imp_null %>%
    autoplot(.) +
    theme(axis.text.x = element_text(angle = 90))

  #do the bootstrapping
  CWM <- trait_np_bootstrap(trait_imp, nrep = 100, sample_size = 200)
  CWM_notiv <- trait_np_bootstrap(trait_imp_null, nrep = 100, sample_size = 200)

  CWM_mean <- trait_summarise_boot_moments(CWM) %>%
    select(Gradient:mean, -n)

  CWM_notiv_mean <- trait_summarise_boot_moments(CWM_notiv) %>%
    select(Gradient:mean, -n) %>%
    rename("mean_noitv" = "mean")

  traitMean <- CWM_mean %>%
    left_join(CWM_notiv_mean, by = c("Gradient", "Site", "PlotID", "trait_trans"))

  #prepare bootstrapped trait data for analyses
  traitMean <- traitMean %>%
    ungroup() %>%
    mutate(trait_trans = factor(trait_trans, levels = c("Plant_Height_cm_log", "Dry_Mass_g_log", "Leaf_Area_cm2_log", "Thickness_mm_log", "SLA_cm2_g", "LDMC", "C_percent", "N_percent", "CN_ratio", "P_percent", "NP_ratio", "dC13_permil", "dN15_permil")))

  return(traitMean)

}
