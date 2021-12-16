# try to make as few targets as possible as each target is cached.
# With many intermediate steps, it uses a lot of disk space.

transformation_plan <- list(

  # import coordinates
  tar_target(
    name = coordinates,
    command = read_excel(path = "clean_data/PFTC4_Svalbard_Coordinates.xlsx") %>%
      filter(Project == "T", Treatment %in% c("C", "B"), !Site %in% c("CAS", "BIS", "DRY")) %>%
      mutate(Site = as.numeric(Site)) %>%
      select(Gradient = Treatment, Site:Longitude_E)
  ),

  # import traits
  tar_target(
    name = traits_raw,
    command = read_csv(file = "clean_data/traits/PFTC4_Svalbard_2018_Gradient_Traits.csv") %>%
      #log transform size and area traits
      mutate(
        value_trans = if_else(
          Trait %in% c(
            "Plant_Height_cm",
            "Dry_Mass_g",
            "Leaf_Area_cm2",
            "Thickness_mm"
          ),
          true = suppressWarnings(log(Value)),# suppress warnings from log(-value) in isotopes (these are calculated but not kept)
          false = Value
        ),
        trait_trans = recode(
          Trait,
          "Plant_Height_cm" = "Plant_Height_cm_log",
          "Dry_Mass_g" = "Dry_Mass_g_log",
          "Leaf_Area_cm2" = "Leaf_Area_cm2_log",
          "Thickness_mm" = "Thickness_mm_log"
        )) %>% distinct(trait_trans)
  ),

  # import community
  tar_target(
    name = comm_raw,
    command = read_csv(file = "clean_data/community/PFTC4_Svalbard_2018_Community_Gradient.csv") %>%
      # remove duplicates
      distinct()
  ),

  # import climate
  tar_target(
    name = climate_data,
    command = read_csv(file = "clean_data/climate/PFTC4_Svalbard_2018_Gradient_Climate.csv")
  ),

  # make bootstrap
  tar_target(
    name = trait_mean,
    command = make_bootstrapping(comm_raw, traits_raw)  %>%
      left_join(coordinates, by = c("Gradient", "Site"))

  )

)