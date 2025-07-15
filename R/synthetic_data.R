synthetic_data <- function() {
  targets <- list(
    # Target: Load raw data
    tar_target(
      synthetic_libraries_,
      {
        ## ---- synthetic libraries
        library(tidyverse)  # for data manipulation and visualisation
        library(synthos)    # for synthetic data generation
        library(sf)         # for spatial data handling and visualisation
        ## ----end
      }
    ),
    tar_target(
      synthetic_global_parameters_,
      {
        ## ---- synthetic global parameters
        assign(x = "data_path", value = "../data/", envir = .GlobalEnv)
        paths <- list(
          data_path = data_path,
          synthetic_path = paste0(data_path, "synthetic/")
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
    tar_target(synthetic_landscape_config_, {
      ## ---- synthetic landscape config
      config <- list(
        seed = 1,
        crs = 4326,
        model = "Exp",
        psill = 1,
        range = 15,
        nugget = 0,
        alpha = 2,
        kappa = 1,
        variance = 1,
        patch_threshold = 1.75,
        reef_width = 0.01,
        years = 1:12,
        ## years = 2000:2011,
        dhw_weight = 0.5,
        cyc_weight = 0.4,
        other_weight = 0.1,
        hcc_cover_range = c(0.1, 0.4),
        hcc_growth = 0.3,
        sc_cover_range = c(0.01, 0.1),
        sc_growth =  0.3
      )
      ## ----end
      config
    }),
    tar_target(synthetic_landscape_spatial_domain_, {
      config <- synthetic_landscape_config_
      ## ---- synthetic landscape spatial domain
      spatial_domain <- st_geometry(
        st_multipoint(
          x = rbind(
            c(0, -10),
            c(3, -10),
            c(10, -20),
            c(1, -21),
            c(2, -16),
            c(0, -10)
          )
        )
      ) |>
        st_set_crs(config$crs) |>
        st_cast("POLYGON")
      set.seed(config$seed)
      spatial_grid <- spatial_domain |>
        st_set_crs(NA) |>
        st_sample(size = 10000, type = "regular") |>
        st_set_crs(config$crs)
      ## ----end
    }),
    tar_target(synthetic_landscape_benthos_reefs_, {
      config <- synthetic_landscape_config_
      spatial_grid <- synthetic_landscape_spatial_domain_
      data_path <- synthetic_global_parameters_$data_path
      ## ---- synthetic landscape benthos reefs
      sf_use_s2(FALSE)
      benthos_reefs_pts <- create_synthetic_reef_landscape(spatial_grid,
        config,
        include_disturbances = TRUE,
        verbose = FALSE
      )
      saveRDS(benthos_reefs_pts, file = paste0(data_path, "synthetic/benthos_reefs_pts.rds"))
      ## ----end
      benthos_reefs_pts
    }),
    tar_target(synthetic_landscape_benthos_reefs_response_scale_, {
      benthos_reefs_pts <- synthetic_landscape_benthos_reefs_
      data_path <- synthetic_global_parameters_$data_path
      ## ---- synthetic landscape benthos reefs response scale
      benthos_reefs_sf <- benthos_reefs_pts |>
        mutate(across(c(HCC, SC, MA), ~ plogis(.))) |>
        dplyr::select(Year, Reef, geometry, everything())
      saveRDS(benthos_reefs_sf, file = paste0(data_path, "synthetic/benthos_reefs_sf.rds"))
      write_csv(benthos_reefs_sf, paste0(data_path, "synthetic/benthos_reefs_sf.csv"))
      ## ----end
      benthos_reefs_sf
    }),
    tar_target(synthetic_landscape_benthos_reefs_locs_, {
      benthos_reefs_pts <- synthetic_landscape_benthos_reefs_
      config <- synthetic_landscape_config_
      data_path <- synthetic_global_parameters_$data_path
      ## ---- synthetic landscape benthos reefs locs
      config <- list(n_locs = 25, n_sites = 2, seed = 123)
      benthos_fixed_locs_sf <- sampling_design_large_scale_fixed(benthos_reefs_pts, config)
      saveRDS(benthos_fixed_locs_sf, file = paste0(
        data_path,
        "synthetic/benthos_fixed_locs_sf.rds"
      ))
      ## ----end
      benthos_fixed_locs_sf 
    }),

    tar_target(synthetic_landscape_benthos_reefs_locs_obs_, {
      benthos_fixed_locs_sf <- synthetic_landscape_benthos_reefs_locs_
      config <- synthetic_landscape_config_
      data_path <- synthetic_global_parameters_$data_path
      ## ---- synthetic landscape benthos reefs locs obs
      config <- list(
        Number_of_transects_per_site = 5,
        Depths = 2,
        Number_of_frames_per_transect = 100,
        Points_per_frame = 5,
        ## Note, the following are on the link scale
        hcc_site_sigma = 0.5, # variability in Sites within Locations
        hcc_transect_sigma = 0.2, # variability in Transects within Sites
        hcc_sigma = 0.1, # random noise

        sc_site_sigma = 0.05, # variability in Sites within Locations
        sc_transect_sigma = 0.02, # variability in Transects within Sites
        sc_sigma = 0.01, # random noise

        ma_site_sigma = 0.5, # variability in Sites within Locations
        ma_transect_sigma = 0.2, # variability in Transects within Sites
        ma_sigma = 0.1 # random noise
      )

      benthos_fixed_locs_obs <- sampling_design_fine_scale_fixed(
        benthos_fixed_locs_sf,
        config
      )
      saveRDS(benthos_fixed_locs_obs, file = paste0(
        data_path,
        "synthetic/benthos_fixed_locs_obs.rds"
      ))
      ## save as csv
      write_csv(
        benthos_fixed_locs_obs,
        paste0(data_path, "synthetic/benthos_fixed_locs_obs.csv")
      )
      ## ----end
      benthos_fixed_locs_obs
    }),
    tar_target(synthetic_landscape_benthos_reefs_locs_obs_disturb_, {
      benthos_fixed_locs_obs <- synthetic_landscape_benthos_reefs_locs_obs_
      benthos_fixed_locs_sf <- synthetic_landscape_benthos_reefs_locs_
      data_path <- synthetic_global_parameters_$data_path
      ## ---- synthetic landscape benthos reefs locs obs disturb
      benthos_fixed_locs_obs_disturb <- benthos_fixed_locs_obs |>
        left_join(
          benthos_fixed_locs_sf |>
            dplyr::select(Reef, Year, CYC, DHW, OTHER) |>
            distinct(),
          by = c("Reef", "Year"),
          relationship = "many-to-many"
        ) 
      saveRDS(benthos_fixed_locs_obs_disturb, file = paste0(
        data_path,
        "synthetic/benthos_fixed_locs_obs_disturb.rds"
      ))
      write_csv(
        benthos_fixed_locs_obs_disturb,
        paste0(data_path, "synthetic/benthos_fixed_locs_obs_disturb.csv")
      )
      ## ----end
      benthos_fixed_locs_obs_disturb 
    }),

    ## Versions of the data with specific issues/challenges

    tar_target(synthetic_replace_reefs_, {
      benthos_fixed_locs_obs <- synthetic_landscape_benthos_reefs_locs_obs_
      data_path <- synthetic_global_parameters_$data_path
      ## ---- synthetic replace reefs
      reef_names <- benthos_fixed_locs_obs |>
        group_by(Reef) |>
        summarise(HCC = mean(HCC)) |>
        filter(HCC == max(HCC) | HCC == min(HCC)) |>
        arrange(HCC) |>
        pull(Reef)
      yrs <- benthos_fixed_locs_obs |> pull(Year) |> unique()
      benthos_fixed_locs_obs_1 <-
        benthos_fixed_locs_obs |>
        filter(
          !(Reef == reef_names[1] & Year > yrs[floor(length(yrs)*3/4)] |
              Reef == reef_names[2] & Year < yrs[floor(length(yrs)*3/4) + 1])
        )
      saveRDS(benthos_fixed_locs_obs_1,
        file = paste0(data_path, "synthetic/benthos_fixed_locs_obs_1.rds")
      )
      write_csv(
        benthos_fixed_locs_obs_1,
        paste0(data_path, "synthetic/benthos_fixed_locs_obs_1.csv")
      )
      ## ----end
    }),
    ## also have a version of this where there are only five reefs
    ## (including the swapped reefs from above)
    tar_target(synthetic_replace_reefs_2_, {
      benthos_fixed_locs_obs_1 <- synthetic_replace_reefs_
      benthos_fixed_locs_obs <- synthetic_landscape_benthos_reefs_locs_obs_
      data_path <- synthetic_global_parameters_$data_path
      ## ---- synthetic replace reefs 2
      set.seed(123)
      reef_names <- benthos_fixed_locs_obs |>
        group_by(Reef) |>
        summarise(HCC = mean(HCC)) |>
        filter(HCC == max(HCC) | HCC == min(HCC)) |>
        arrange(HCC) |>
        pull(Reef)
      reef_names_rnd <- benthos_fixed_locs_obs |>
        group_by(Reef) |>
        summarise(HCC = mean(HCC)) |>
        arrange(HCC) |>
        filter(!Reef %in% reef_names) |>
        pull(Reef) |>
        sample(3) |>
        c(reef_names)
      benthos_fixed_locs_obs_2 <-
        benthos_fixed_locs_obs_1 |>
        filter(Reef %in% reef_names_rnd)

      saveRDS(benthos_fixed_locs_obs_2,
        file = paste0(data_path, "synthetic/benthos_fixed_locs_obs_2.rds")
      )
      write_csv(
        benthos_fixed_locs_obs_2,
        paste0(data_path, "synthetic/benthos_fixed_locs_obs_2.csv")
      )
      ## ----end
    }),
    ## Temporal gap.  Years 3-4 are missing
    tar_target(synthetic_temporal_gap_, {
      benthos_fixed_locs_obs <- synthetic_landscape_benthos_reefs_locs_obs_
      data_path <- synthetic_global_parameters_$data_path
      ## ---- synthetic temporal gap
      set.seed(123)
      reef_names <- benthos_fixed_locs_obs |>
        pull(Reef) |>
        unique() |> 
        sample(5)
      benthos_fixed_locs_obs_3 <-
        benthos_fixed_locs_obs |>
        filter(
          !(Reef %in% reef_names & Year %in% 3:5)
        )

      saveRDS(benthos_fixed_locs_obs_3,
        file = paste0(data_path, "synthetic/benthos_fixed_locs_obs_3.rds")
      )
      write_csv(
        benthos_fixed_locs_obs_3,
        paste0(data_path, "synthetic/benthos_fixed_locs_obs_3.csv")
      )
      ## ----end
    }),
    tar_target(synthetic_temporal_gap_2_, {
      benthos_fixed_locs_obs <- synthetic_landscape_benthos_reefs_locs_obs_
      data_path <- synthetic_global_parameters_$data_path
      ## ---- synthetic temporal gap 2
      benthos_fixed_locs_obs_4 <-
        benthos_fixed_locs_obs |>
        filter(
          !(Year %in% 3:5)
        )
      saveRDS(benthos_fixed_locs_obs_4,
        file = paste0(data_path, "synthetic/benthos_fixed_locs_obs_4.rds")
      )
      write_csv(
        benthos_fixed_locs_obs_4,
        paste0(data_path, "synthetic/benthos_fixed_locs_obs_4.csv")
      )
      ## ----end
    }),
    ## covariates
    tar_target(synthetic_covariates_, {
      benthos_fixed_locs_obs_disturb <- synthetic_landscape_benthos_reefs_locs_obs_disturb_
      benthos_reefs_sf <- synthetic_landscape_benthos_reefs_response_scale_
      data_path <- synthetic_global_parameters_$data_path
      ## ---- synthetic covariates
      benthos_fixed_locs_obs_disturb_list <- list(
        benthos_fixed_locs_obs_disturb = benthos_fixed_locs_obs_disturb,
        benthos_reefs_sf = benthos_reefs_sf
      )
      saveRDS(benthos_fixed_locs_obs_disturb_list,
        file = paste0(data_path, "synthetic/benthos_fixed_locs_obs_disturb_list.rds")
      )
      ## ----end
    })
  )
  return(targets)
}
