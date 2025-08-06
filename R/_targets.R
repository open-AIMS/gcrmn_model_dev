
# Load the targets package
library(targets)
library(tarchetypes)

# Set global options
tar_option_set(
  packages = c("tidyverse", "sf", "synthos",
    "glmmTMB", "emmeans", "DHARMa", "patchwork",
    "brms", "rstan", "bayesplot", "tidybayes",
    "posterior", "gbm", "dbarts", "HDInterval"),  # Load required packages
  format = "rds"                    # Default storage format
)
## lapply(packages, library, character.only = TRUE)
 
source("helper_functions.R")    # Load the modelling of incomplete spatial script
source("synthetic_data.R")      # Load the synthetic data generation script
source("site_replacement.R")    # Load the modelling of site replacement script
source("missing_years.R")       # Load the modelling of missing years script
source("incomplete_spatial.R")  # Load the modelling of incomplete spatial script

list(
  ## tar_target(synthetic_data_output,
  synthetic_data(),
  ## ),
  ## tar_target(site_replacement_output,
  site_replacement(),
  ## )
  missing_years(),
  helper_functions(),
  incomplete_spatial()
)
