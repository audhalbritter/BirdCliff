library("targets")
library("tarchetypes")
#remotes::install_github("Between-the-Fjords/dataDownloader")

tar_option_set(packages = c("dataDownloader", "tidyverse", "readxl", "traitstrap", "vegan", "ggvegan", "kableExtra", "viridis", "patchwork", "broom", "nlme", "lme4", "broom.mixed", "glue", "MuMIn", "Hmisc"))
# "ape", "egg", "FD"
#"performance",

# source target plans and functions
tar_source()


#Combine target plans
combined_plan <- c(
  download_plan,
  transformation_plan,
  analysis_plan,
  figures_plan,
  si_figures_plan#,
  #manuscript_plan
)

#clean up subplans
# rm(
#   neotoma_download_plan,
#   cleaning_plan,
#   analysis_plan,
#   figures_plan,
#   manuscript_plan
# )

#return combined plan
combined_plan
