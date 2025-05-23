
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

 
source("synthetic_data.R")  # Load the synthetic data generation script
source("site_replacement.R")  # Load the modelling of site replacement script

list(
  synthetic_data(),
  site_replacement()
)
