source("model_functions.R")
## source("helper_functions.R")
missing_years <- function() {
  targets <- list(
    tar_target(missing_years_link,
      lnk <- synthetic_covariates_),
    # Target: Load raw data
    tar_target(
      missing_years_libraries_,
      {
        ## ---- site replacement libraries
        library(tidyverse) # for data manipulation and visualisation
        library(sf) # for spatial data handling and visualisation
        library(glmmTMB) # for frequentist GLMM
        library(emmeans) # for post-hoc analysis
        library(DHARMa) # for model diagnostics
        library(patchwork)
        library(brms) # for Bayesian GLMM
        library(rstan) # for Bayesian modelling
        library(bayesplot) # for Bayesian plotting
        library(tidybayes) # for tidy bayesian analysis
        library(posterior) # for posterior analysis
        library(gbm) # for gradient boosted regression trees
        library(dbarts) # for Bayesian additive regression trees
        library(HDInterval) # for Bayesian credible intervals
        ## ----end
      }
    ),
    tar_target(
      missing_years_extra_functions_,
      {
        ## ---- site replacement functions
        source("model_functions.R")
        ## source("helper_functions.R")
        ## ----end
      }
    ),
    tar_target(
      missing_years_global_parameters_,
      {
        ## ---- site replacement global parameters
        assign(x = "data_path", value = "../data/", envir = .GlobalEnv)
        assign(x = "output_path", value = "../output/", envir = .GlobalEnv)
        paths <- list(
          data_path = data_path,
          synthetic_path = paste0(data_path, "synthetic/"),
          output_path = output_path,
          fig_path = paste0(output_path, "figures/")
        )
        lapply(paths, function(x) {
          if (!dir.exists(x)) {
            dir.create(x)
          }
        })
        ## ----end
        paths
      }
    ),
    ## Import data =====================================================  
    ## The processing of all data and all sampled data occurs in the
    ## site_replacement.R script, so rather than repeat those steps here,
    ## we will instead just load the processed data

    ## Temporal gap.  Years 3-4 are missing for some sites
    tar_target(
      read_sampled_reefs_data_3_,
      {
        data_path <- missing_years_global_parameters_$data_path
        benthos_reefs_sf <- read_all_reefs_data_
        tmp <- synthetic_temporal_gap_
        ## ---- read sampled reefs data 3
        benthos_fixed_locs_obs_3 <- readRDS(
          file = paste0(
            data_path,
            "synthetic/benthos_fixed_locs_obs_3.rds"
          )
        )
        benthos_fixed_locs_obs_3 <- benthos_fixed_locs_obs_3 |>
          left_join(
            benthos_reefs_sf |>
              st_drop_geometry() |>
              dplyr::select(Year, Reef, CYC, DHW, OTHER) |>
              group_by(Year, Reef) |>
              summarise(across(c(CYC, DHW, OTHER), mean)),
            by = c("Year", "Reef")
          )
        ## ----end
        benthos_fixed_locs_obs_3 
      }
    ),
    tar_target(sampled_reefs_data_3_plot_, {
      benthos_fixed_locs_obs_3 <- read_sampled_reefs_data_3_
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      ## ---- sampled reefs data 3 plot
      g <- benthos_fixed_locs_obs_3 |>
        dplyr::select(-MA, -SC) |>
        pivot_wider(
          id_cols = c(Reef, Longitude, Latitude, Site, Transect),
          names_from = Year,
          values_from = HCC
        ) |>
        pivot_longer(
          cols = -c(Reef, Longitude, Latitude, Site, Transect),
          names_to = "Year",
          values_to = "HCC"
        ) |>
        mutate(Year = as.numeric(Year)) |>
        ggplot() +
        geom_line(aes(
          y = HCC, x = Year, colour = Site,
          group = interaction(Site, Transect)
        )) +
        facet_wrap(~Reef) +
        theme_bw()

      ggsave(
        filename = paste0(
          fig_path, "R_sampled_reefs_3_plot.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
      ## ----end
    }),
    ## Temporal gap.  Years 3-4 are missing for all sites
    tar_target(
      read_sampled_reefs_data_4_,
      {
        data_path <- missing_years_global_parameters_$data_path
        benthos_reefs_sf <- read_all_reefs_data_
        tmp <- synthetic_temporal_gap_2_
        ## ---- read sampled reefs data 4
        benthos_fixed_locs_obs_4 <- readRDS(
          file = paste0(
            data_path,
            "synthetic/benthos_fixed_locs_obs_4.rds"
          )
        )
        benthos_fixed_locs_obs_4 <- benthos_fixed_locs_obs_4 |>
          left_join(
            benthos_reefs_sf |>
              st_drop_geometry() |>
              dplyr::select(Year, Reef, CYC, DHW, OTHER) |>
              group_by(Year, Reef) |>
              summarise(across(c(CYC, DHW, OTHER), mean)),
            by = c("Year", "Reef")
          )
        ## ----end
        benthos_fixed_locs_obs_4 
      }
    ),
    tar_target(sampled_reefs_data_4_plot_, {
      benthos_fixed_locs_obs_4 <- read_sampled_reefs_data_4_
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      ## ---- sampled reefs data 4 plot
      yrs <- benthos_fixed_locs_obs_4 |>
        pull(Year) |>
        full_seq(period = 1)
      g <- benthos_fixed_locs_obs_4 |>
        complete(Year = yrs, nesting(Reef, Longitude, Latitude, Site, Transect)) |>
        ggplot() +
        geom_line(aes(
          y = HCC, x = Year, colour = Site,
          group = interaction(Site, Transect)
        )) +
        facet_wrap(~Reef) +
        theme_bw()

      ggsave(
        filename = paste0(
          fig_path, "R_sampled_reefs_4_plot.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
      ## ----end
    })
    
  )
  return(targets)
}
