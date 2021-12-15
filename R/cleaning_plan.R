# try to make as few targets as possible as each target is cached.
# With many intermediate steps, it uses a lot of disk space.

cleaning_plan <- list(
  # import traits
  tar_target(
    name = traits_gradient,
    command = read_csv(file = "clean_data/traits/PFTC4_Svalbard_2018_Gradient_Traits.csv")
  ),

  # import community
  tar_target(
    name = comm_gradient,
    command = read_csv(file = "clean_data/community/PFTC4_Svalbard_2018_Community_Gradient.csv")
  ),

  # import climate
  tar_target(
    name = climate_data,
    command = read_csv(file = "clean_data/climate/PFTC4_Svalbard_2018_Gradient_Climate.csv")
  )

)
