#rename traits to fancy names for figures
#rename traits to fancy names for figures

fancy_trait_name_dictionary <- function(dat){

  dat <- dat %>%
    # make fancy names for figures
    mutate(trait_fancy = recode(trait_trans,
                              C_percent = "C %",
                              CN_ratio = "CN",
                              dC13_permil = "δC13 ‰",
                              dN15_permil = "δN15 ‰",
                              Plant_Height_cm_log = "Height cm",
                              Dry_Mass_g_log = "Dry mass g",
                              LDMC = "LDMC",
                              Leaf_Area_cm2_log = "Area cm2",
                              N_percent = "N %",
                              NP_ratio = "NP",
                              P_percent = "P %",
                              SLA_cm2_g = "SLA cm2/g",
                              Thickness_mm_log = "Thickness mm",
                              Shoot_Length_cm_log = "Shoot length cm",
                              Shoot_Length_Green_cm_log = "Green length cm",
                              Shoot_ratio = "Shoot ratio",
                              SSL_cm_g = "SSL cm/g",
                              WHC_g_g = "WHC g/g")) |>
    mutate(figure_names = case_match(trait_trans,
                                     "Plant_Height_cm_log" ~ "Size~-~Height~cm",
                                     "Dry_Mass_g_log" ~ "Size~-~Dry~mass~g",
                                     "Leaf_Area_cm2_log" ~ "Size~-~Area~cm^2",
                                     "Thickness_mm_log" ~ "Size~-~Thickness~mm",
                                     "LDMC" ~ "LES~-~LDMC",
                                     "SLA_cm2_g" ~ "LES~-~SLA~cm^2*g^{-1}",
                                     "C_percent" ~ "LES~-~C~'%'",
                                     "N_percent" ~ "LES~-~N~'%'",
                                     "CN_ratio" ~ "LES~-~CN",
                                     "P_percent" ~ "LES~-~P~'%'",
                                     "NP_ratio" ~ "LES~-~NP",
                                     "dC13_permil" ~ "I~-~δ^{13}~C~'‰'",
                                     "dN15_permil" ~ "I~-~δ^{15}~N~'‰'")) |>

    # add class
    mutate(class = case_when(trait_trans %in% c("Plant_Height_cm_log", "Dry_Mass_g_log", "Leaf_Area_cm2_log", "Thickness_mm_log", "Shoot_Length_cm_log", "Shoot_Length_Green_cm_log") ~ "Size",
                             trait_trans %in% c("dC13_permil", "dN15_permil") ~ "Isotopes",
                             TRUE ~ "Leaf economics"),
           class = factor(class, levels = c("Size", "Leaf economics", "Isotopes")))

  return(dat)
}

