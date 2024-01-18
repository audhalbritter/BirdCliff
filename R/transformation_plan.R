# try to make as few targets as possible as each target is cached.
# With many intermediate steps, it uses a lot of disk space.

transformation_plan <- list(


  # import coordinates
  tar_target(
    name = coordinates,
    command = read_csv(coords) %>%
      mutate(Site = as.numeric(Site),
             Gradient = factor(Gradient, levels = c("C", "B")))
  ),

  # soil
  tar_target(
    name = cn_data,
    command = read_csv(soil_cn) |>
      mutate(GS = paste0(Gradient, Site),
             Gradient = factor(Gradient, levels = c("C", "B"))) |>
      select(-Unit, -Weight_mg)
  ),

  # import community
  tar_target(
    name = comm_raw,
    command = read_csv(community) %>%
      # remove duplicates
      distinct() %>%
      mutate(GS = paste0(Gradient, Site),
             Gradient = factor(Gradient, levels = c("C", "B")))
  ),

  # calculate diversity indices
  tar_target(
    name = diversity_grad,
    command = comm_raw %>%
      group_by(Gradient, Site, PlotID, Elevation_m) %>%
      summarise(Richness = n(),
                Diversity = diversity(Cover),
                Evenness = Diversity/log(Richness),
                sumAbundance = sum(Cover)) %>%
      pivot_longer(cols = Richness:sumAbundance, names_to = "DiversityIndex", values_to = "Value")
  ),

  # import traits
  tar_target(
    name = traits_raw,
    command = read_csv(traits) %>%
      # remove bryo
      filter(Project != "Bryophytes") %>%
      # remove wet mass, correlated with dry mass
      filter(Trait != "Wet_Mass_g") %>%
      #log transform size and area traits
      mutate(
        value_trans = if_else(
          Trait %in% c(
            "Plant_Height_cm",
            "Dry_Mass_g",
            "Leaf_Area_cm2",
            "Leaf_Thickness_mm"
          ),
          true = suppressWarnings(log(Value)),# suppress warnings from log(-value) in isotopes (these are calculated but not kept)
          false = Value
        ),
        trait_trans = recode(
          Trait,
          "Plant_Height_cm" = "Plant_Height_cm_log",
          "Dry_Mass_g" = "Dry_Mass_g_log",
          "Leaf_Area_cm2" = "Leaf_Area_cm2_log",
          "Leaf_Thickness_mm" = "Thickness_mm_log"
        )) |>
      # order traits
      mutate(trait_trans = factor(trait_trans, levels = c("Plant_Height_cm_log", "Dry_Mass_g_log", "Leaf_Area_cm2_log", "Thickness_mm_log", "LDMC", "SLA_cm2_g", "C_percent", "N_percent", "CN_ratio", "P_percent", "NP_ratio", "dC13_permil", "dN15_permil")),
             Gradient = factor(Gradient, levels = c("C", "B")))
  ),

  # import bryophyte traits
  tar_target(
    name = bryo_traits_raw,
    command = read_csv(traits) %>%
      # filter for bryo
      filter(Project == "Bryophytes") %>%
      # filter for important traits
      filter(!Trait %in% c("Wet_Mass_g", "Dry_Mass_g")) %>%
      #log transform shoot length traits
      mutate(
        value_trans = if_else(
          Trait %in% c(
            "Shoot_Length_cm",
            "Shoot_Length_Green_cm"
          ),
          true = suppressWarnings(log(Value)),# suppress warnings from log(-value) in isotopes (these are calculated but not kept)
          false = Value
        ),
        trait_trans = recode(
          Trait,
          "Shoot_Length_cm" = "Shoot_Length_cm_log",
          "Shoot_Length_Green_cm" = "Shoot_Length_Green_cm_log"
        )) %>%

      # merge niphotrichum and polytrichum species
      mutate(Taxon = case_when(Taxon == "niphotrichum canescens" ~ "niphotrichum sp",
                               Taxon == "polytrichum piliferum" ~ "polytrichum sp",
                               TRUE ~ Taxon)) |>
      select(-c(Elevation_m:Longitude_E)) %>%
      left_join(coordinates %>%
                  group_by(Gradient, Site) %>%
                  summarise(Elevation_m = mean(Elevation_m)), by = c("Gradient", "Site")) |>
      # order traits
      mutate(trait_trans = factor(trait_trans, levels = c("Shoot_Length_cm_log", "Shoot_Length_Green_cm_log", "Shoot_ratio", "WHC_g_g", "SSL_cm_g", "C_percent", "N_percent", "CN_ratio", "P_percent", "NP_ratio", "dC13_permil", "dN15_permil")),
             Gradient = factor(Gradient, levels = c("C", "B")))
  ),

  # import climate
  tar_target(
    name = environment,
    command = read_csv(climate) |>
      mutate(Gradient = factor(Gradient, levels = c("C", "B")),
             GS = paste0(Gradient, Site)) |>
      select(Gradient:GS) |>
      left_join(coordinates, by = c("Gradient", "Site", "PlotID")) |>
      bind_rows(cn_data)
  ),

  # BOOTSTRAPPING
  # make trait impute
  tar_target(
    name = trait_impute,
    command = make_trait_impute(comm_raw, traits_raw)

  ),

  # trait coverage
  tar_target(
    name = trait_coverage,
    command = fortify_filled_trait(trait_impute) |>
      ungroup() |>
      #filter(Trait == "CN_ratio") |>
      complete(.id, level, trait_trans, fill = list(s = 0)) |> # .id does not exist anymore...
      complete(level, trait_trans, fill = list(s = 0)) |>
      filter(level == "PlotID") |>
      group_by(Gradient, trait_trans) |>
      # prob = 0.25 gives 75% of the plots
      # also run prob = 0.5 for 50% of the plots
      summarise(q = quantile(s, prob = 0.25))

  ),

  tar_target(
    name = trait_null_impute,
    command = make_trait_null_impute(comm_raw, traits_raw)

  ),

  # make bootstrap
  tar_target(
    name = trait_mean,
    command = make_bootstrapping(trait_impute, trait_null_impute) %>%
      left_join(coordinates, by = c("Gradient", "Site", "PlotID")) %>%
      # join climate data
      left_join(environment |>
                  pivot_wider(names_from = Variable, values_from = Value),
                by = c("Gradient", "Site", "PlotID", "Elevation_m", "Longitude_E", "Latitude_N"))
    )
)
