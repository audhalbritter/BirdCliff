#rename traits to fancy names for figures
#rename traits to fancy names for figures

fancy_trait_name_dictionary <- function(dat){

  dat <- dat %>%
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
                              WHC_g_g = "WHC g/g"))

  return(dat)
}
