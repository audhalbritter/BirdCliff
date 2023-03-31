library("targets")
library("tarchetypes")
#remotes::install_github("Between-the-Fjords/dataDownloader")

tar_option_set(packages = c("dataDownloader", "tidyverse", "readxl", "traitstrap", "vegan", "ggvegan", "kableExtra", "viridis", "patchwork", "broom", "ape", "nlme", "lme4", "broom.mixed", "glue", "MuMIn", "Hmisc", "egg", "FD"))
#"performance",

# source target plans - can also construct plans directly in this file.
source("R/download_plan.R")
source("R/transformation_plan.R")
source("R/analysis_plan.R")
source("R/figures_plan.R")
source("R/si_figures_plan.R")
# source("R/manuscript_plan.R")

# source functions
source("R/functions/fancy_traits.R")
source("R/functions/bootstrapping.R")
source("R/functions/correlation_plot.R")
source("R/functions/ordination.R")
source("R/functions/trait-analysis.R")
source("R/functions/trait_figures.R")
source("R/functions/itv_analysis.R")
source("R/functions/individual_species.R")
source("R/functions/microclimate.R")


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
