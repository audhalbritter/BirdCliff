#

download_plan <- list(
  #select datasets to download
  # Leaf trait data
  tar_target(
    name = traits,
    command = download_PFTC_data(country = "Svalbard",
                                 datatype = "trait",
                                 path = "clean_data/traits")
    ),

  # community data
  tar_target(
    name = community,
    command = download_PFTC_data(country = "Svalbard",
                                 datatype = "community",
                                 path = "clean_data/community")
  )
)
