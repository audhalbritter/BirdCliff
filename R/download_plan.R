#

download_plan <- list(

  # meta data
  tar_target(
    name = coords,
    command = get_file(node = "smbqh",
                       file = "PFTC4_Svalbard_Coordinates_Gradient.csv",
                       path = "clean_data",
                       remote_path = "MetaData"),
    format = "file"
  ),

  # trait data
  tar_target(
    name = traits,
    command = get_file(node = "smbqh",
                       file = "PFTC4_Svalbard_2018_Gradient_Traits.csv",
                       path = "clean_data/traits/",
                       remote_path = "Traits"),
    format = "file"
    ),

  # community data
  tar_target(
    name = community,
    command = get_file(node = "smbqh",
                       file = "PFTC4_Svalbard_2018_Community_Gradient.csv",
                       path = "clean_data/community",
                       remote_path = "Community"),
    format = "file"
  ),

  # climate data
  tar_target(
    name = climate,
    command = get_file(node = "smbqh",
                       file = "PFTC4_Svalbard_2018_Gradient_Climate.csv",
                       path = "clean_data/climate",
                       remote_path = "Climate"),
    format = "file"
  ),

  # soil data
  tar_target(
    name = soil_cn,
    command = get_file(node = "smbqh",
                       file = "PFTC4_Svalbard_2022_Gradient_Clean_Soil_CN.csv",
                       path = "clean_data/soil",
                       remote_path = "Soil"),
    format = "file"
  )
)
