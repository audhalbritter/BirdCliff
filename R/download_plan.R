#

download_plan <- list(

  # meta data
  tar_target(
    name = coords,
    command = get_file(node = "smbqh",
                       file = "PFTC4_Svalbard_Coordinates.xlsx",
                       path = "clean_data",
                       remote_path = "MetaData")
  ),

  # trait data
  tar_target(
    name = traits,
    command = get_file(node = "smbqh",
                       file = "PFTC4_Svalbard_2018_Gradient_Traits.csv",
                       path = "clean_data/traits",
                       remote_path = "Traits")
    ),

  # community data
  tar_target(
    name = community,
    command = get_file(node = "smbqh",
                       file = "PFTC4_Svalbard_2018_Community_Gradient.csv",
                       path = "clean_data/community",
                       remote_path = "Community")
  ),
  #PFTC4_Svalbard_2018_Community_Structure_Gradient.csv

  # climate data
  tar_target(
    name = climate,
    command = get_file(node = "smbqh",
                       file = "PFTC4_Svalbard_2018_Gradient_Climate.csv",
                       path = "clean_data/climate",
                       remote_path = "Climate")
  )
)
