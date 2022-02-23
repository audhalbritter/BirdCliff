library("targets")
library("tarchetypes")
#remotes::install_github("Between-the-Fjords/dataDownloader")

tar_option_set(packages = c("dataDownloader", "tidyverse", "readxl", "traitstrap", "vegan", "ggvegan", "kableExtra", "viridis", "patchwork", "broom", "ape", "nlme", "Hmisc", "lme4", "broom.mixed", "MuMIn"))
#"performance",

# source target plans - can also construct plans directly in this file.
source("R/download_plan.R")
source("R/transformation_plan.R")
source("R/analysis_plan.R")
source("R/si_figures_plan.R")
source("R/manuscript_plan.R")

# source functions
source("R/functions/bootstrapping.R")
source("R/functions/fancy_traits.R")
source("R/functions/ordination.R")
source("R/functions/intra_vs_inter.R")
source("R/functions/inter_intra_anova.R")
source("R/functions/correlation_plot.R")
source("R/functions/trait-analysis.R")
source("R/functions/trait_output.R")
source("R/functions/ITV.R")
source("R/functions/individual_species.R")


#Combine target plans
combined_plan <- c(
  download_plan,
  transformation_plan,
  analysis_plan,
  figures_plan,
  manuscript_plan
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
