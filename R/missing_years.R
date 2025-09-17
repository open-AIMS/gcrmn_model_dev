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
    }),
    
    ## Data preparations ==============================================
    ## All sampled reefs ----------------------------------------------
    
    ## Temporal gap.  Years 3-4 are missing for some sites ------------
    tar_target(missing_years_data_prep_3_, {
      benthos_fixed_locs_obs_3 <- read_sampled_reefs_data_3_
      ## ---- sampled data prep 3
      benthos_fixed_locs_obs_3 <- benthos_fixed_locs_obs_3 |>
        mutate(
          fYear = as.factor(Year),
          Reef = as.factor(Reef),
          Site = interaction(Reef, Site),
          Transect = interaction(Site, Transect),
          cover = HCC/100
        )
      benthos_fixed_locs_obs_3 
      ## ----end
      benthos_fixed_locs_obs_3 
    }),
    tar_target(missing_years_newdata_3_, {
      benthos_fixed_locs_obs_3 <- missing_years_data_prep_3_
      ## ---- newdata 3
      newdata_3 <-
        benthos_fixed_locs_obs_3 |>
        tidyr::expand(fYear)
      newdata_3
      ## ----end
      newdata_3
    }),
    tar_target(missing_years_newdata_3b_, {
      benthos_fixed_locs_obs_3 <- missing_years_data_prep_3_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- newdata 3b
      newdata_3b <-
        benthos_fixed_locs_obs_3 |> 
        ## mutate(Longitude = st_coordinates(geometry)[, 1],
        ##   Latitude = st_coordinates(geometry)[, 2]) |>
        ## st_drop_geometry() |>
        group_by(Year, Reef) |>
        summarise(
          HCC = mean(HCC),
          CYC = mean(CYC),
          DHW = mean(DHW),
          OTHER = mean(OTHER),
          Longitude = mean(Longitude),
          Latitude = mean(Latitude),
          HCC = HCC/100
          )
      newdata_3b
      saveRDS(newdata_3b,
        file = paste0(
          data_path,
          "synthetic/newdata_3b.rds"
        )
      )
      ## ----end
      newdata_3b
    }),
    ## Temporal gap.  Years 3-4 are missing for all sites -------------
    tar_target(missing_years_data_prep_4_, {
      benthos_fixed_locs_obs_4 <- read_sampled_reefs_data_4_
      ## ---- sampled data prep 4
      benthos_fixed_locs_obs_4 <- benthos_fixed_locs_obs_4 |>
        mutate(
          fYear = as.factor(Year),
          Reef = as.factor(Reef),
          Site = interaction(Reef, Site),
          Transect = interaction(Site, Transect),
          cover = HCC/100
        )
      ## ----end
      benthos_fixed_locs_obs_4 
    }),
    tar_target(missing_years_newdata_4_, {
      benthos_fixed_locs_obs_4 <- missing_years_data_prep_4_
      ## ---- newdata 4
      newdata_4 <-
        benthos_fixed_locs_obs_4 |>
        tidyr::expand(fYear)
      ## ----end
      newdata_4
    }),
    tar_target(missing_years_newdata_4b_, {
      benthos_fixed_locs_obs_4 <- missing_years_data_prep_4_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- newdata 4b
      newdata_4b <-
        benthos_fixed_locs_obs_4 |> 
        ## mutate(Longitude = st_coordinates(geometry)[, 1],
        ##   Latitude = st_coordinates(geometry)[, 2]) |>
        ## st_drop_geometry() |>
        group_by(Year, Reef) |>
        summarise(
          HCC = mean(HCC),
          CYC = mean(CYC),
          DHW = mean(DHW),
          OTHER = mean(OTHER),
          Longitude = mean(Longitude),
          Latitude = mean(Latitude),
          HCC = HCC/100
          )
      newdata_4b
      saveRDS(newdata_4b,
        file = paste0(
          data_path,
          "synthetic/newdata_4b.rds"
        )
      )
      ## ----end
      newdata_4b
    }),

    ## Models =========================================================

    ## Temporal gap.  Years 3-4 are missing for some sites ------------
    ## simple ---------------------------------------------------------
    tar_target(mod_simple_3_, {
      benthos_fixed_locs_obs_3 <- missing_years_data_prep_3_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- simple_3
      mod_simple_3 <- benthos_fixed_locs_obs_3 |>
        group_by(Year, Site) |>
        summarise(
          Mean = mean(cover),
          Median = median(cover)
        ) |>
        ungroup() |>
        group_by(Year) |> 
        summarise(
          Mean = mean(Mean),
          Median = median(Median)
        ) 
      saveRDS(mod_simple_3,
        file = paste0(data_path, "synthetic/mod_simple_3.rds")
      ) 
      ## ----end
      mod_simple_3
    }),

    tar_target(mod_glmmTMB_3_sample_data_, {
      benthos_fixed_locs_obs_3 <- missing_years_data_prep_3_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_3 sample data
      benthos_fixed_locs_obs_3 
      benthos_fixed_locs_obs_3 |>
        group_by(Year) |>
        summarise(
          Mean = mean(cover),
          Median = median(cover)
          )
      ## ----end
      benthos_fixed_locs_obs_3
    }),
    tar_target(mod_glmmTMB_3_newdata_3_, {
      newdata_3 <- missing_years_newdata_3_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_3 newdata 3
      newdata_3 
      ## ----end
      newdata_3
    }),
    tar_target(mod_glmmTMB_3_newdata_3b_, {
      newdata_3b <- missing_years_newdata_3b_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_3 newdata 3b
      newdata_3b 
      unique(newdata_3b$Year)
      ## ----end
      newdata_3b
    }),
    tar_target(mod_glmmTMB_3_newdata_, {
      newdata <- site_replacements_newdata_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_3 newdata
      newdata 
      ## ----end
      newdata
    }),
    tar_target(mod_glmmTMB_3_newdata_b_, {
      newdata_b <- site_replacements_newdata_b_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_3 newdata_b
      newdata_b 
      unique(newdata_b$Year)
      ## ----end
      newdata_b
    }),

    
    ## glmmTMB --------------------------------------------------------
    tar_target(mod_glmmTMB_3_, {
      benthos_fixed_locs_obs_3 <- missing_years_data_prep_3_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_3
      mod_glmmTMB_3 <- glmmTMB(cover ~ fYear + (1 | Site) + (1 | Transect),
        data = benthos_fixed_locs_obs_3,
        family = "beta_family"
      ) 
      saveRDS(mod_glmmTMB_3,
        file = paste0(data_path, "synthetic/mod_glmmTMB_3.rds")
      ) 
      ## ----end
      mod_glmmTMB_3
    }),
    tar_target(dharma_mod_glmmTMB_3_, {
      mod_glmmTMB_3 <- mod_glmmTMB_3_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_3_dharma
      glmmTMB_3_dharma <- mod_glmmTMB_3 |> 
        simulateResiduals(n = 1000)
      saveRDS(glmmTMB_3_dharma,
        file = paste0(data_path, "synthetic/glmmTMB_3_dharma.rds")
      ) 
      ## ----end
      glmmTMB_3_dharma
    }),
    tar_target(dharma_mod_glmmTMB_3_plot_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      glmmTMB_3_dharma <- dharma_mod_glmmTMB_3_$glmmTMB_3_dharma
      ## ---- glmmTMB_3_dharma plot
      glmmTMB_3_dharma <- readRDS(
        file = paste0(data_path, "synthetic/glmmTMB_3_dharma.rds")
      )
      g <- wrap_elements(~testUniformity(glmmTMB_3_dharma)) +
        wrap_elements(~plotResiduals(glmmTMB_3_dharma)) +
        wrap_elements(~testDispersion(glmmTMB_3_dharma))
      ggsave(
        filename = paste0(
          fig_path, "R_dharma_mod_glmmTMB_3.png"
        ),
        g,
        width = 10, height = 4, dpi = 72
      )
      ## ----end
    }),
    tar_target(emmeans_mod_glmmTMB_3, {
      mod_glmmTMB_3 <- mod_glmmTMB_3_
      newdata_3 <- missing_years_newdata_3_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_3_emmeans
      glmmTMB_3_sum <- 
        mod_glmmTMB_3 |> emmeans(~fYear, at = newdata_3, type = "response") |>
        as.data.frame() |> 
        rename(median = response, lower = asymp.LCL, upper = asymp.UCL) |>
        mutate(Year = as.numeric(as.character(fYear))) |>
        mutate(type = "glmmTMB")
      saveRDS(glmmTMB_3_sum,
        file = paste0(data_path, "synthetic/glmmTMB_3_sum.rds")
      ) 
      ## ----end
      glmmTMB_3_sum 
    }),
    tar_target(emmeans_mod_glmmTMB_3_plot_, {
      glmmTMB_3_sum <- emmeans_mod_glmmTMB_3
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_simple_3 <- mod_simple_3_
      ## ---- glmmTMB_3_emmeans plot
      g1 <- glmmTMB_3_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "glmmTMB")) +
        geom_line(data = mod_simple_3,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_3,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        ## geom_line(data = all_sampled_sum,
        ##   aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      g2 <- glmmTMB_3_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "glmmTMB")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "true mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "true median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_glmmTMB_3.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    
    tar_target(pred_1_mod_glmmTMB_3_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_3 <- mod_glmmTMB_3_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_3b <- missing_years_newdata_3b_
      mod_simple_3 <- mod_simple_3_ #sampled_simple_raw_means_
      ## ---- glmmTMB_3_pred_1
      newdata <- newdata_3b
      true_sum <- mod_simple_3
      glmmTMB_3_pred <- pred_glmmTMB(mod_glmmTMB_3,
        type = 1,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_3"
      )
      ## ----end
      glmmTMB_3_pred
    }),
    tar_target(mse_1_mod_glmmTMB_3_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_3 <- mod_glmmTMB_3_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_3b <- missing_years_newdata_3b_
      ## ---- glmmTMB_3_mse_1
      glmmTMB_3_mse <- mse_glmmTMB(mod_glmmTMB_3,
        newdata = newdata_3b, type = 1, model_type = ""
      )
      ## ----end
      glmmTMB_3_mse
    }),
    tar_target(pred_2_mod_glmmTMB_3_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_3 <- mod_glmmTMB_3_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata <- site_replacements_newdata_
      benthos_reefs_temporal_summary_0 <- read_all_temporal_summary_
      ## ---- glmmTMB_3_pred_2
      newdata <- newdata
      true_sum <- benthos_reefs_temporal_summary_0
      glmmTMB_3_pred <- pred_glmmTMB(mod_glmmTMB_3,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_3"
      )
      ## ----end
      glmmTMB_3_pred
    }),
    tar_target(mse_2_mod_glmmTMB_3_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_3 <- mod_glmmTMB_3_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      ## ---- glmmTMB_3_mse_2
      glmmTMB_3_mse <- mse_glmmTMB(mod_glmmTMB_3,
        newdata = newdata_b, type = 2, model_type = ""
      )
      ## ----end
      glmmTMB_3_mse
    }),
    
    ## glmmTMB + covariates -------------------------------------------
    tar_target(mod_glmmTMB_3b_, {
      benthos_fixed_locs_obs_3 <- missing_years_data_prep_3_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_3b
      mod_glmmTMB_3b <- glmmTMB(cover ~ fYear + CYC + DHW + OTHER + Latitude + Longitude +
                                   (1 | Site) + (1 | Transect),
        data = benthos_fixed_locs_obs_3,
        family = "beta_family"
      )
      saveRDS(mod_glmmTMB_3b,
        file = paste0(data_path, "synthetic/mod_glmmTMB_3b.rds")
      ) 
      ## ----end
      mod_glmmTMB_3b
    }),
    tar_target(dharma_mod_glmmTMB_3b_, {
      DHARMa_glmmTMB <- DHARMa_glmmTMB_
      mod_glmmTMB_3b <- mod_glmmTMB_3b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- glmmTMB_3b_dharma
      glmmTMB_3b_dharma <- DHARMa_glmmTMB(mod_glmmTMB_3b,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_3b")
      ## ----end
      glmmTMB_3b_dharma
    }),
    tar_target(pred_1_mod_glmmTMB_3b_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_3b <- mod_glmmTMB_3b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_3b <- missing_years_newdata_3b_
      mod_simple_3 <- mod_simple_3_ #sampled_simple_raw_means_
      ## ---- glmmTMB_3b_pred_1
      true_sum <- mod_simple_3
      newdata <- newdata_3b
      glmmTMB_3b_pred <- pred_glmmTMB(mod_glmmTMB_3b,
        type = 1,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_3"
      )
      ## ----end
      glmmTMB_3b_pred
    }),
    tar_target(mse_1_mod_glmmTMB_3b_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_3 <- mod_glmmTMB_3b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_3b <- missing_years_newdata_3b_
      ## ---- glmmTMB_3b_mse_1
      glmmTMB_3_mse <- mse_glmmTMB(mod_glmmTMB_3,
        newdata = newdata_3b, type = 1, model_type = "covariates"
      )
      ## ----end
      glmmTMB_3_mse
    }),
    tar_target(pred_2_mod_glmmTMB_3b_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_3b <- mod_glmmTMB_3b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      benthos_reefs_temporal_summary_0 <- read_all_temporal_summary_
      ## ---- glmmTMB_3b_pred_2
      newdata <- newdata_b
      true_sum <- benthos_reefs_temporal_summary_0
      glmmTMB_3b_pred <- pred_glmmTMB(mod_glmmTMB_3b,
        type = 2,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_3"
      )
      ## ----end
      glmmTMB_3b_pred
    }),
    tar_target(mse_2_mod_glmmTMB_3b_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_3b <- mod_glmmTMB_3b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      ## ---- glmmTMB_3b_mse_2
      glmmTMB_3b_mse <- mse_glmmTMB(mod_glmmTMB_3b,
        newdata = newdata_b, type = 2, model_type = "covariates"
      )
      ## ----end
      glmmTMB_3b_mse
    }),

    ## brms -----------------------------------------------------------
    tar_target(mod_brms_3_, {
      benthos_fixed_locs_obs_3 <- missing_years_data_prep_3_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- brms_pre_3
      benthos_fixed_locs_obs_3 |>
        group_by(fYear) |>
        summarise(mean_abundance = qlogis(mean(cover)),
          median_abundance = qlogis(median(cover)),
          sd_abundance = sd(qlogis(cover)) ,
          mad_abundance = mad(qlogis(cover)) 
        ) 
      priors <- c(
        prior(normal(0, 0.5), class = "b"),
        prior(normal(-1.5, 0.5), class = "Intercept"),
        prior(student_t(3, 0, 0.5), class = "sd")
      )
      mod_form <- bf(cover ~ fYear + (1 | Site) + (1 | Transect),
        family = "Beta"
      )
      ## ----end
      ## ---- brms_3
      mod_brms_3 <- brm(mod_form,
        data = benthos_fixed_locs_obs_3,
        iter = 5000,
        warmup = 1000,
        chains = 3,
        cores = 3,
        prior = priors,
        thin =  5,
        control = list(adapt_delta = 0.99),
        backend = "cmdstanr"
      )
      saveRDS(mod_brms_3,
        file = paste0(data_path, "synthetic/mod_brms_3.rds")
      ) 
      ## ----end
      mod_brms_3
    }),
    tar_target(brms_trace_3_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_brms_3 <- mod_brms_3_
      ## ---- brms_trace_3
      mod_brms_3 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_3.rds")
      )
      vars <- mod_brms_3 |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_3$fit |> stan_trace(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_trace_3.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ac_3_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_brms_3 <- mod_brms_3_
      ## ---- brms_ac_3
      mod_brms_3 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_3.rds")
      )
      vars <- mod_brms_3 |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_3$fit |> stan_ac(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ac_3.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_rhat_3_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_brms_3 <- mod_brms_3_
      ## ---- brms_rhat_3
      mod_brms_3 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_3.rds")
      )
      g <- mod_brms_3$fit |> stan_rhat()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_rhat_3.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ess_3_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_brms_3 <- mod_brms_3_
      ## ---- brms_ess_3
      mod_brms_3 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_3.rds")
      )
      g <- mod_brms_3$fit |> stan_ess()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ess_3.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ppc_3_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_brms_3 <- mod_brms_3_
      ## ---- brms_ppc_3
      mod_brms_3 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_3.rds")
      )
      g <- mod_brms_3 |> pp_check( type='dens_overlay', ndraws=100)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ppc_3.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(pred_1_mod_brms_3_, {
      pred_brms <- pred_brms_
      mod_brms_3 <- mod_brms_3_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_3b <- missing_years_newdata_3b_
      mod_simple_3 <- mod_simple_3_ #sampled_simple_raw_means_
      ## ---- brms_3_pred_1
      newdata <- newdata_3b
      true_sum <- mod_simple_3
      brms_3_pred <- pred_brms(mod_brms_3,
        type = 1,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_3"
      )
      ## ----end
      brms_3_pred
    }),
    tar_target(mse_1_mod_brms_3_, {
      mse_brms <- mse_brms_
      mod_brms_3 <- mod_brms_3_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_3b <- missing_years_newdata_3b_
      ## ---- brms_3_mse_1
      brms_3_mse <- mse_brms(mod_brms_3,
        newdata = newdata_3b, type = 1, model_type = ""
      )
      ## ----end
      brms_3_mse
    }),
    tar_target(pred_2_mod_brms_3_, {
      pred_brms <- pred_brms_
      mod_brms_3 <- mod_brms_3_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata <- site_replacements_newdata_
      benthos_reefs_temporal_summary_0 <- read_all_temporal_summary_
      ## ---- brms_3_pred_2
      newdata <- newdata
      true_sum <- benthos_reefs_temporal_summary_0
      brms_3_pred <- pred_brms(mod_brms_3,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_3"
      )
      ## ----end
      brms_3_pred
    }),
    tar_target(mse_2_mod_brms_3_, {
      mse_brms <- mse_brms_
      mod_brms_3 <- mod_brms_3_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      ## ---- brms_3_mse_2
      brms_3_mse <- mse_brms(mod_brms_3,
        newdata = newdata_b, type = 2, model_type = ""
      )
      ## ----end
      brms_3_mse
    }),
    
    ## brms + covariates ----------------------------------------------
    tar_target(mod_brms_3b_, {
      benthos_fixed_locs_obs_3 <- missing_years_data_prep_3_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- brms_pre_3b
      benthos_fixed_locs_obs_3 |>
        group_by(fYear) |>
        summarise(mean_abundance = qlogis(mean(cover)),
          median_abundance = qlogis(median(cover)),
          sd_abundance = sd(qlogis(cover)) ,
          mad_abundance = mad(qlogis(cover)) 
        ) 
      priors <- c(
        prior(normal(0, 0.5), class = "b"),
        prior(normal(-1.5, 0.5), class = "Intercept"),
        prior(student_t(3, 0, 0.5), class = "sd")
      )
      mod_form <- bf(cover ~ fYear + scale(CYC) + scale(DHW) +
                       scale(OTHER) + scale(Latitude) + scale(Longitude) +
                       (1 | Site) + (1 | Transect),
        family = "Beta"
      )
      ## ----end
      ## ---- brms_3b
      mod_brms_3b <- brm(mod_form,
        data = benthos_fixed_locs_obs_3,
        iter = 5000,
        warmup = 1000,
        chains = 3,
        cores = 3,
        prior = priors,
        thin =  5,
        control = list(adapt_delta = 0.99),
        backend = "cmdstanr"
      )
      saveRDS(mod_brms_3b,
        file = paste0(data_path, "synthetic/mod_brms_3b.rds")
      ) 
      ## ----end
      mod_brms_3b
    }),
    tar_target(brms_trace_3b_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_brms_3b <- mod_brms_3b_
      ## ---- brms_trace_3b
      mod_brms_3b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_3b.rds")
      )
      vars <- mod_brms_3b |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_3b$fit |> stan_trace(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_trace_3b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ac_3b_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_brms_3b <- mod_brms_3b_
      ## ---- brms_ac_3b
      mod_brms_3b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_3b.rds")
      )
      vars <- mod_brms_3b |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_3b$fit |> stan_ac(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ac_3b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_rhat_3b_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_brms_3b <- mod_brms_3b_
      ## ---- brms_rhat_3b
      mod_brms_3b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_3b.rds")
      )
      g <- mod_brms_3b$fit |> stan_rhat()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_rhat_3b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ess_3b_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_brms_3b <- mod_brms_3b_
      ## ---- brms_ess_3b
      mod_brms_3b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_3b.rds")
      )
      g <- mod_brms_3b$fit |> stan_ess()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ess_3b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ppc_3b_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_brms_3b <- mod_brms_3b_
      ## ---- brms_ppc_3b
      mod_brms_3b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_3b.rds")
      )
      g <- mod_brms_3b |> pp_check( type='dens_overlay', ndraws=100)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ppc_3b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(pred_1_mod_brms_3b_, {
      pred_brms <- pred_brms_
      mod_brms_3b <- mod_brms_3b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_3b <- missing_years_newdata_3b_
      mod_simple_3 <- mod_simple_3_ #sampled_simple_raw_means_
      ## ---- brms_3b_pred_1
      newdata <- newdata_3b
      true_sum <- mod_simple_3
      brms_3b_pred <- pred_brms(mod_brms_3b,
        type = 1,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_3"
      )
      ## ----end
      brms_3b_pred
    }),
    tar_target(mse_1_mod_brms_3b_, {
      mse_brms <- mse_brms_
      mod_brms_3b <- mod_brms_3b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_3b <- missing_years_newdata_3b_
      ## ---- brms_3b_mse_1
      brms_3b_mse <- mse_brms(mod_brms_3b,
        newdata = newdata_3b, type = 1, model_type = "covariates"
      )
      ## ----end
      brms_3b_mse
    }),
    tar_target(pred_2_mod_brms_3b_, {
      pred_brms <- pred_brms_
      mod_brms_3b <- mod_brms_3b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      benthos_reefs_temporal_summary_0 <- read_all_temporal_summary_
      ## ---- brms_3b_pred_2
      newdata <- newdata_b
      true_sum <- benthos_reefs_temporal_summary_0
      brms_3b_pred <- pred_brms(mod_brms_3b,
        type = 2,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_3"
      )
      ## ----end
      brms_3b_pred
    }),
    tar_target(mse_2_mod_brms_3b_, {
      mse_brms <- mse_brms_
      mod_brms_3b <- mod_brms_3b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      ## ---- brms_3b_mse_2
      brms_3b_mse <- mse_brms(mod_brms_3b,
        newdata = newdata_b, type = 2, model_type = "covariates"
      )
      ## ----end
      brms_3b_mse
    }),

    ## stan -----------------------------------------------------------
    tar_target(mod_stan_3_, {
      source("model_functions.R")
      benthos_fixed_locs_obs_3 <- missing_years_data_prep_3_
      data_path <- missing_years_global_parameters_$data_path
      site_extra_functions_
      ## ---- stan_pre_3
      benthos_fixed_locs_obs_3 <-
        benthos_fixed_locs_obs_3 |>
        mutate(
          cYear = fYear,
          grid_id = factor(Reef),
          cSite = factor(interaction(Reef, Site)),
          cReplicate = ifelse(is.na(Transect),
            interaction(Site, Year),
            interaction(cSite, Transect)),
          Cover = cover,
          area = 1,
          sum = 1
        )
      saveRDS(benthos_fixed_locs_obs_3,
        file = paste0(data_path, "synthetic/saveRDS(benthos_fixed_locs_obs_3_forstan.rds")
      ) 
      stan_data <- prepare_data_for_stan(benthos_fixed_locs_obs_3, yrs = NULL)
      model_stan <- cmdstanr::cmdstan_model(stan_file = "model1.stan")
      ## ----end
      ## ---- stan_3
      mod_stan_3 <- model_stan$sample(
        data = stan_data,
        seed = 123,
        iter_sampling = 5000,
        iter_warmup = 1000,
        thin = 5,
        chains = 3,
        parallel_chains = 3,
        adapt_delta = 0.99,
        output_dir = paste0(data_path, "synthetic/"),
      )
      saveRDS(mod_stan_3,
        file = paste0(data_path, "synthetic/mod_stan_3.rds")
      ) 
      ## ----end
      mod_stan_3
    }),
    tar_target(stan_trace_3_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_stan_3 <- mod_stan_3_
      ## ---- stan_trace_3
      mod_stan_3 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_3.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_3$draws(variables = c("beta", "phi", "sd_1", "sd_2", "sd_3")) |>
        mcmc_trace() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_trace_3.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ac_3_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_stan_3 <- mod_stan_3_
      ## ---- stan_ac_3
      mod_stan_3 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_3.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_3$draws(variables = c("beta", "phi", "sd_1", "sd_2", "sd_3")) |>
        mcmc_acf() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ac_3.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_rhat_3_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_stan_3 <- mod_stan_3_
      ## ---- stan_rhat_3
      mod_stan_3 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_3.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_3 |> bayesplot::rhat() |> 
        mcmc_rhat_hist() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_rhat_3.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ess_3_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_stan_3 <- mod_stan_3_
      ## ---- stan_ess_3
      mod_stan_3 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_3.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_3 |> bayesplot::neff_ratio() |> 
        mcmc_neff_hist() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ess_3.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ppc_3_, {
      benthos_fixed_locs_obs_3 <- missing_years_data_prep_3_
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_stan_3 <- mod_stan_3_
      ## ---- stan_ppc_3
      mod_stan_3 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_3.rds")
      )
      g <- 
        bayesplot::pp_check(
          benthos_fixed_locs_obs_3$cover,
          mod_stan_3$draws("ypred", format = "matrix")[1:100, ],
          ppc_dens_overlay
        ) +
        theme_classic()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ppc_3.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(pred_1_mod_stan_3_, {
      pred_stan <- pred_stan_
      mod_stan_3 <- mod_stan_3_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_3 <- missing_years_newdata_3_
      mod_simple_3 <- mod_simple_3_ #sampled_simple_raw_means_
      ## ---- stan_3_pred_1
      newdata <- newdata_3
      true_sum <- mod_simple_3
      stan_3_pred <- pred_stan(mod_stan_3,
        type = 1,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_stan_3"
      )
      ## ----end
      stan_3_pred
    }),
    tar_target(mse_1_mod_stan_3_, {
      mse_stan <- mse_stan_
      mod_stan_3 <- mod_stan_3_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_3b <- missing_years_newdata_3b_
      ## ---- stan_3_mse_1
      stan_3_mse <- mse_stan(mod_stan_3,
        newdata = newdata_3b, type = 1, model_type = ""
      )
      ## ----end
      stan_3_mse
    }),
    tar_target(pred_2_mod_stan_3_, {
      pred_stan <- pred_stan_
      mod_stan_3 <- mod_stan_3_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata <- site_replacements_newdata_
      benthos_reefs_temporal_summary_0 <- read_all_temporal_summary_
      ## ---- stan_3_pred_2
      newdata <- newdata
      true_sum <- benthos_reefs_temporal_summary_0
      stan_3_pred <- pred_stan(mod_stan_3,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_stan_3"
      )
      ## ----end
      stan_3_pred
    }),
    tar_target(mse_2_mod_stan_3_, {
      mse_stan <- mse_stan_
      mod_stan_3 <- mod_stan_3_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      ## ---- stan_3_mse_2
      stan_3_mse <- mse_stan(mod_stan_3,
        newdata = newdata_b, type = 2, model_type = ""
      )
      ## ----end
      stan_3_mse
    }),
    
    ## gbm ------------------------------------------------------------
    tar_target(mod_gbm_3_, {
      benthos_fixed_locs_obs_3 <- missing_years_data_prep_3_
      data_path <- missing_years_global_parameters_$data_path
      ## ---- gbm_3
      mod_gbm_3 <- gbm(cover ~ fYear,
        data =  benthos_fixed_locs_obs_3,
        distribution = "gaussian",
        n.trees = 10000,
        interaction.depth = 5,
        shrinkage = 0.001,
        bag.fraction = 0.5,
        cv.folds = 5, n.cores = 1,
        verbose = TRUE
      )
      saveRDS(mod_gbm_3,
        file = paste0(data_path, "synthetic/mod_gbm_3.rds")
      ) 
      ## ----end
      ## ---- gbm_post_3
      n.trees <- gbm.perf(mod_gbm_3, method = "cv")
      ## ----end
      list(mod_gbm_3 = mod_gbm_3, n.trees = n.trees)
    }),
    tar_target(pdp_gbm_3_, {
      mod_gbm_3 <- mod_gbm_3_$mod_gbm_3
      n.trees <- mod_gbm_3_$n.trees
      newdata_3 <- missing_years_newdata_3_
      data_path <- missing_years_global_parameters_$data_path
      ## ---- gbm_pdp_3
      gbm_3_sum <- newdata_3 |>
        mutate(median = predict(mod_gbm_3, newdata_3, n.trees = n.trees, type = "response")) |> 
        mutate(Year = as.numeric(as.character(fYear))) |>
        mutate(type = "gbm")
      saveRDS(gbm_3_sum,
        file = paste0(data_path, "synthetic/gbm_3_sum.rds")
      ) 
      ## ----end
      gbm_3_sum
    }),
    tar_target(pred_1_mod_gbm_3_, {
      pred_gbm <- pred_gbm_
      mod_gbm_3 <- mod_gbm_3_$mod_gbm_3
      n.trees <- mod_gbm_3_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_3 <- missing_years_newdata_3_
      mod_simple_3 <- mod_simple_3_ #sampled_simple_raw_means_
      ## ---- gbm_3_pred_1
      newdata <- newdata_3
      true_sum <- mod_simple_3
      gbm_3_pred <- pred_gbm(mod_gbm_3,
        n.trees = n.trees,
        type = 1,
        model_type = "",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_3"
      )
      ## ----end
      gbm_3_pred
    }),
    tar_target(mse_1_mod_gbm_3_, {
      mse_gbm <- mse_gbm_
      mod_gbm_3 <- mod_gbm_3_$mod_gbm_3
      n.trees <- mod_gbm_3_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_3b <- missing_years_newdata_3b_
      ## ---- gbm_3_mse_1
      gbm_3_mse <- mse_gbm(mod_gbm_3,
        n.trees = n.trees,
        newdata = newdata_3b, type = 1, model_type = ""
      )
      ## ----end
      gbm_3_mse
    }),
    tar_target(pred_2_mod_gbm_3_, {
      pred_gbm <- pred_gbm_
      mod_gbm_3 <- mod_gbm_3_$mod_gbm_3
      n.trees <- mod_gbm_3_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata <- site_replacements_newdata_
      benthos_reefs_temporal_summary_0 <- read_all_temporal_summary_
      ## ---- gbm_3_pred_2
      newdata <- newdata
      true_sum <- benthos_reefs_temporal_summary_0
      gbm_3_pred <- pred_gbm(mod_gbm_3,
        n.trees =  n.trees,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_3"
      )
      ## ----end
      gbm_3_pred
    }),
    tar_target(mse_2_mod_gbm_3_, {
      mse_gbm <- mse_gbm_
      mod_gbm_3 <- mod_gbm_3_$mod_gbm_3
      n.trees <- mod_gbm_3_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      ## ---- gbm_3_mse_2
      gbm_3_mse <- mse_gbm(mod_gbm_3,
        n.trees = n.trees,
        newdata = newdata_b, type = 2, model_type = ""
      )
      ## ----end
      gbm_3_mse
    }),
    
    ## gbm + covartiates ----------------------------------------------
    tar_target(mod_gbm_3b_, {
      benthos_fixed_locs_obs_3 <- missing_years_data_prep_3_
      data_path <- missing_years_global_parameters_$data_path
      ## ---- gbm_3b
      mod_gbm_3b <- gbm(cover ~ fYear + Latitude + Longitude + CYC + DHW + OTHER,
        data =  benthos_fixed_locs_obs_3,
        distribution = "gaussian",
        var.monotone = c(0, 0, 0, -1, -1, -1),
        n.trees = 10000,
        interaction.depth = 5,
        shrinkage = 0.001,
        bag.fraction = 0.5,
        cv.folds = 5, n.cores = 1,
        verbose = TRUE
      )
      saveRDS(mod_gbm_3b,
        file = paste0(data_path, "synthetic/mod_gbm_3b.rds")
      )
      ## ----end
      ## ---- gbm_post_3b
      n.trees <- gbm.perf(mod_gbm_3b, method = "cv")
      ## ----end
      list(mod_gbm_3b = mod_gbm_3b, n.trees = n.trees)
    }),
    tar_target(pred_1_mod_gbm_3b_, {
      pred_gbm <- pred_gbm_
      mod_gbm_3b <- mod_gbm_3b_$mod_gbm_3b
      n.trees <- mod_gbm_3b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_3b <- missing_years_newdata_3b_
      mod_simple_3 <- mod_simple_3_ #sampled_simple_raw_means_
      ## ---- gbm_3b_pred_1
      newdata <- newdata_3b
      true_sum <- mod_simple_3
      gbm_3b_pred <- pred_gbm(mod_gbm_3b,
        n.trees = n.trees,
        type = 1,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_3"
      )
      ## ----end
      gbm_3b_pred
    }),
    tar_target(mse_1_mod_gbm_3b_, {
      mse_gbm <- mse_gbm_
      mod_gbm_3b <- mod_gbm_3b_$mod_gbm_3b
      n.trees <- mod_gbm_3b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_3b <- missing_years_newdata_3b_
      ## ---- gbm_3b_mse_1
      gbm_3b_mse <- mse_gbm(mod_gbm_3b,
        n.trees = n.trees,
        newdata = newdata_3b, type = 1, model_type = "covariates"
      )
      ## ----end
      gbm_3b_mse
    }),
    tar_target(pred_2_mod_gbm_3b_, {
      pred_gbm <- pred_gbm_
      mod_gbm_3b <- mod_gbm_3b_$mod_gbm_3b
      n.trees <- mod_gbm_3b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      benthos_reefs_temporal_summary_0 <- read_all_temporal_summary_
      ## ---- gbm_3b_pred_2
      newdata <- newdata_b
      true_sum <- benthos_reefs_temporal_summary_0
      gbm_3b_pred <- pred_gbm(mod_gbm_3b,
        n.trees =  n.trees,
        type = 2,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_3"
      )
      ## ----end
      gbm_3b_pred
    }),
    tar_target(mse_2_mod_gbm_3b_, {
      mse_gbm <- mse_gbm_
      mod_gbm_3b <- mod_gbm_3b_$mod_gbm_3b
      n.trees <- mod_gbm_3b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      ## ---- gbm_3b_mse_2
      gbm_3b_mse <- mse_gbm(mod_gbm_3b,
        n.trees = n.trees,
        newdata = newdata_b, type = 2, model_type = "covariates"
      )
      ## ----end
      gbm_3b_mse
    }),

    ## dbarts --------------------------------------------------------
    tar_target(mod_dbarts_3_, {
      benthos_fixed_locs_obs_3 <- missing_years_data_prep_3_
      data_path <- missing_years_global_parameters_$data_path
      newdata_3 <- missing_years_newdata_3_
      newdata_3b <- missing_years_newdata_3b_
      newdata_3_a <- missing_years_data_prep_3_
      newdata <- site_replacements_newdata_
      newdata_a <- site_replacements_newdata_b_
      ## ---- dbarts_3
      print(head(benthos_fixed_locs_obs_3))
      mod_dbarts_3 <- bart2(log(cover) ~ fYear,
        data =   benthos_fixed_locs_obs_3,
        keepTrees = TRUE
      )
      saveRDS(mod_dbarts_3,
        file = paste0(data_path, "synthetic/mod_dbarts_3.rds")
      ) 
      ## ----end
      ## Unfortunately, the next part must be in the same tar_target
      ## due to the way dbarts stores pointers - they cannot be stored
      ## ---- dbarts_pred_3
      preds <- predict(mod_dbarts_3, newdata_3, type = "ev") |>
        exp() 
      saveRDS(preds,
        file = paste0(data_path, "synthetic/mod_dbarts_3_preds.rds")
      ) 
      preds_sum <- preds |>
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_3_preds_sum.rds")
      ) 

      newdata.1 <- newdata_3_a |>
        group_by(Year, Reef) |>
        summarise(
          HCC = mean(HCC, na.rm = TRUE),
          Latitude = mean(Latitude),
          Longitude = mean(Longitude),
          CYC = mean(CYC, na.rm = TRUE),
          DHW = mean(DHW, na.rm = TRUE),
          OTHER = mean(OTHER, na.rm = TRUE)
        ) |>
        mutate(fYear = factor(Year)) |> 
        ungroup()
      preds_a <- predict(mod_dbarts_3, newdata.1, type = "ev") |>
        exp() 
      saveRDS(preds_a,
        file = paste0(data_path, "synthetic/mod_dbarts_3_a_preds.rds")
      ) 
      
      preds_0 <- predict(mod_dbarts_3, newdata, type = "ev") |>
        exp()
      saveRDS(preds_0,
        file = paste0(data_path, "synthetic/mod_dbarts_0_preds.rds")
      ) 
      preds_0_sum <- preds_0 |>
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_0_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_0_preds_sum.rds")
      ) 
      newdata_0_a <- newdata_a |>
        ungroup() |>
        mutate(fYear = factor(Year))
      preds_0_a <- predict(mod_dbarts_3, newdata_0_a, type = "ev") |>
        exp()
      saveRDS(preds_0_a,
        file = paste0(data_path, "synthetic/mod_dbarts_0_a_preds.rds")
      ) 
      ## ----end
      list(
        mod_dbarts_3 = mod_dbarts_3,
        preds = preds,
        preds_sum = preds_sum,
        preds_a = preds_a,
        preds_0 = preds_0,
        preds_0_sum = preds_0_sum,
        preds_0_a = preds_0_a
      )
    }),
    tar_target(pred_1_mod_dbarts_3_, {
      pred_dbarts <- pred_dbarts_
      preds <- mod_dbarts_3_$preds_sum
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_3 <- missing_years_newdata_3_
      mod_simple_3 <- mod_simple_3_ #sampled_simple_raw_means_
      ## ---- dbarts_3_pred_1
      newdata <- newdata_3
      true_sum <- mod_simple_3
      dbarts_3_pred <- pred_dbarts(
        preds,
        type = 1,
        model_type = "",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_3"
      )
      ## ----end
      dbarts_3_pred
    }),
    tar_target(mse_1_mod_dbarts_3_, {
      mse_dbarts <- mse_dbarts_
      preds <- mod_dbarts_3_$preds_a
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      newdata_3b <- missing_years_newdata_3b_
      ## ---- dbarts_3b_mse_1
      dbarts_3_mse <- mse_dbarts(preds,
        newdata = newdata_3b, type = 1, model_type = ""
      )
      ## ----end
      dbarts_3_mse
    }),
    tar_target(pred_2_mod_dbarts_3_, {
      pred_dbarts <- pred_dbarts_
      preds <- mod_dbarts_3_$preds_0_sum
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata <- site_replacements_newdata_
      benthos_reefs_temporal_summary_0 <- read_all_temporal_summary_
      ## ---- dbarts_3_pred_2
      newdata <- newdata
      true_sum <- benthos_reefs_temporal_summary_0
      dbarts_3_pred <- pred_dbarts(
        preds,
        type = 2,
        model_type = "",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_3"
      )
      ## ----end
      dbarts_3_pred
    }),
    tar_target(mse_2_mod_dbarts_3_, {
      mse_dbarts <- mse_dbarts_
      preds <- mod_dbarts_3_$preds_a
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_ 
      ## ---- dbarts_3_mse_2
      dbarts_3_mse <- mse_dbarts(preds,
        newdata = newdata_b, type = 2, model_type = ""
      )
      ## ----end
      dbarts_3_mse
    }),
    
    ## dbarts + covariates --------------------------------------------
    tar_target(mod_dbarts_3b_, {
      benthos_fixed_locs_obs_3 <- missing_years_data_prep_3_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_3 <- missing_years_newdata_3_
      newdata_3b <- missing_years_newdata_3b_
      newdata_b <- site_replacements_newdata_b_
      ## ---- dbarts_3b
      print(head(benthos_fixed_locs_obs_3))
      mod_dbarts_3b <- bart2(log(cover) ~ fYear + Latitude + Longitude + CYC + DHW + OTHER,
        data =   benthos_fixed_locs_obs_3,
        keepTrees = TRUE
      )
      saveRDS(mod_dbarts_3b,
        file = paste0(data_path, "synthetic/mod_dbarts_3b.rds")
      ) 
      ## ----end
      ## Unfortunately, the next part must be in the same tar_target
      ## due to the way dbarts stores pointers - they cannot be stored
      ## ---- dbarts_pred_3b
      newdata_3b <- newdata_3b |>
        mutate(fYear = factor(Year))
      preds_3b <- predict(mod_dbarts_3b, newdata_3b, type = "ev") |>
        exp() 
      saveRDS(preds_3b,
        file = paste0(data_path, "synthetic/mod_dbarts_3b_preds_3b.rds")
      ) 
      preds_3b_sum <- preds_3b |> summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_3b_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_3b_preds_3b_sum.rds")
      )
      
      newdata_b <- newdata_b |>
        mutate(fYear = factor(Year))
      preds_b <- predict(mod_dbarts_3b, newdata_b, type = "ev") |>
        exp() 
      saveRDS(preds_b,
        file = paste0(data_path, "synthetic/mod_dbarts_3b_preds_b.rds")
      ) 
      preds_b_sum <- preds_b |> 
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_b_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_3b_preds_b_sum.rds")
      ) 

      ## ----end
      list(
        mod_dbarts_3 = mod_dbarts_3b,
        preds_3b = preds_3b,
        preds_3b_sum = preds_3b_sum,
        preds_b = preds_b,
        preds_b_sum = preds_b_sum
      )
    }),
    tar_target(pred_1_mod_dbarts_3b_, {
      pred_dbarts <- pred_dbarts_
      preds_3b <- mod_dbarts_3b_$preds_3b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_3b <- missing_years_newdata_3b_
      mod_simple_3 <- mod_simple_3_ #sampled_simple_raw_means_
      ## ---- dbarts_3b_pred_1
      newdata <- newdata_3b
      true_sum <- mod_simple_3
      dbarts_3b_pred <- pred_dbarts(
        preds_3b,
        type = 1,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_3"
      )
      ## ----end
      dbarts_3b_pred
    }),
    tar_target(mse_1_mod_dbarts_3b_, {
      mse_dbarts <- mse_dbarts_
      preds_3b <- mod_dbarts_3b_$preds_3b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_3b <- missing_years_newdata_3b_
      ## ---- dbarts_3b_mse_1
      dbarts_3_mse <- mse_dbarts(preds_3b,
        newdata = newdata_3b, type = 1, model_type = "covariates"
      )
      ## ----end
      dbarts_3_mse
    }),
    tar_target(pred_2_mod_dbarts_3b_, {
      pred_dbarts <- pred_dbarts_
      preds_b <- mod_dbarts_3b_$preds_b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      benthos_reefs_temporal_summary_0 <- read_all_temporal_summary_
      ## ---- dbarts_3b_pred_2
      newdata <- newdata_b
      true_sum <- benthos_reefs_temporal_summary_0
      dbarts_3_pred <- pred_dbarts(
        preds_b, 
        type = 2, 
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_3"
      )
      ## ----end
      dbarts_3_pred
    }),
    tar_target(mse_2_mod_dbarts_3b_, {
      mse_dbarts <- mse_dbarts_
      preds_b <- mod_dbarts_3b_$preds_b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      ## ---- dbarts_3b_mse_2
      dbarts_3_mse <- mse_dbarts(preds_b,
        newdata = newdata_b, type = 2, model_type = "covariates"
      )
      ## ----end
      dbarts_3_mse
    }),

    ## xgboost + covariates -------------------------------------------
    tar_target(mod_xgboost_3b_prep_, {
      benthos_fixed_locs_obs_3 <- missing_years_data_prep_3_
      data_path <- missing_years_global_parameters_$data_path
      ## ---- xgboost_3_prep
      data_train_3b <- benthos_fixed_locs_obs_3 |>
        dplyr::select(cover, Year, Latitude, Longitude, CYC, DHW, OTHER) 
      ## ----end
      data_train_3b
    }),
    tar_target(mod_xgboost_3b_tune_, {
      data_train_3b <- mod_xgboost_3b_prep_
      data_path <- missing_years_global_parameters_$data_path
      ## ---- xgboost_3_tune
      ## Define the recipe
      tune_recipe <- recipe(cover ~ ., data = data_train_3b) |> 
        step_dummy(all_nominal_predictors())
      ##Define the model
      tune_model <- boost_tree(learn_rate = tune(),
        trees = tune(), 
        min_n = tune(), 
        tree_depth = tune()) |>    # Model type
        set_engine("xgboost") |>   # Model engine
        set_mode("regression")     # Model mode
      ## Define the workflow
      tune_workflow <- workflow() |> 
        add_recipe(tune_recipe) |> 
        add_model(tune_model)
      tune_grid_values <-
        grid_space_filling(learn_rate(),
          trees(),
          tree_depth(),
          min_n(),
          size = 20,
          type = "max_entropy")
      ##Run the hyper parameters tuning
      tuned_results <- tune_grid(tune_workflow,
        resamples = vfold_cv(data_train_3b, v = 5),
        grid = tune_grid_values)
      ## Get best set of parameters
      model_hyperparams <-
        select_best(tuned_results, metric = "rmse") |> 
        select(-".config") |> 
        as_tibble() |>
        mutate(nb_training = nrow(data_train_3b),
          grid_size = nrow(tune_grid_values)) 
      ## ----end
     list(
        tune_recipe = tune_recipe,
        tune_model = tune_model,
        tune_workflow = tune_workflow,
        tune_grid_values = tune_grid_values,
        tuned_results = tuned_results,
        model_hyperparams = model_hyperparams
     ) 
    }),
    tar_target(mod_xgboost_3b_fit_, {
      data_train_3b <- mod_xgboost_3b_prep_
      data_path <- missing_years_global_parameters_$data_path
      tune_recipe = mod_xgboost_3b_tune_$tune_recipe
      tune_model = mod_xgboost_3b_tune_$tune_model
      tune_workflow = mod_xgboost_3b_tune_$tune_workflow
      tune_grid_values = mod_xgboost_3b_tune_$tune_grid_values
      tuned_results = mod_xgboost_3b_tune_$tuned_results
      model_hyperparams = mod_xgboost_3b_tune_$model_hyperparams
      ## ---- xgboost_3_fit
      ## Redefine the model (with hyper parameters)
      tune_model <-
        boost_tree(learn_rate = model_hyperparams$learn_rate,
          trees = model_hyperparams$trees, 
          min_n = model_hyperparams$min_n, 
          tree_depth = model_hyperparams$tree_depth) |> # Model type
        set_engine("xgboost") |> # Model engine
        set_mode("regression") # Model mode
      ## Redefine the workflow
      tune_workflow <- workflow() |>
        add_recipe(tune_recipe) |> 
        add_model(tune_model)
      ## Fit the final model
      final_fitted_3b <- tune_workflow |>
        fit(data_train_3b)
      ## ----end
      final_fitted_3b
    }),
    tar_target(pred_1_mod_xgboost_3b_, {
      pred_xgboost <- pred_xgboost_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_3b <- missing_years_newdata_3b_
      final_fitted_3b <- mod_xgboost_3b_fit_
      mod_simple_3 <- mod_simple_3_ #sampled_simple_raw_means_
      ## ---- xgboost_3b_pred_1
      true_sum <- mod_simple_3
      newdata <- newdata_3b
      xgboost_3b_pred <- pred_xgboost(
        final_fitted_3b,
        type = 1,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_xgboost_3"
      )
      ## ----end
      xgboost_3b_pred
    }),
    tar_target(mse_1_mod_xgboost_3b_, {
      mse_xgboost <- mse_xgboost_
      final_fitted_3b <- mod_xgboost_3b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_3b <- missing_years_newdata_3b_
      ## ---- xgboost_3b_mse_1
      xgboost_3_mse <- mse_xgboost(final_fitted_3b,
        newdata = newdata_3b, type = 1, model_type = "covariates"
      )
      ## ----end
      xgboost_3_mse
    }),
    tar_target(pred_2_mod_xgboost_3b_, {
      pred_xgboost <- pred_xgboost_
      final_fitted_3b <- mod_xgboost_3b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      benthos_reefs_temporal_summary_0 <- read_all_temporal_summary_
      ## ---- xgboost_3b_pred_2
      newdata <- newdata_b
      true_sum <- benthos_reefs_temporal_summary_0
      xgboost_3_pred <- pred_xgboost(final_fitted_3b,
        type = 2,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_xgboost_3"
      )
      ## ----end
      xgboost_3_pred
    }),
    tar_target(mse_2_mod_xgboost_3b_, {
      mse_xgboost <- mse_xgboost_
      final_fitted_3b <- mod_xgboost_3b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      ## ---- xgboost_3b_mse_2
      xgboost_3_mse <- mse_xgboost(final_fitted_3b,
        newdata = newdata_b, type = 2, model_type = "covariates"
      )
      ## ----end
      xgboost_3_mse
    }),

    ## Comparisons ----------------------------------------------------

    tar_target(mse_1_mod_pymc_barts_3b_file_, {
      fl <- paste0(
        missing_years_global_parameters_$data_path,
        "modelled/pymc_bart_3b_mse_1.csv"
      )
      if (!file.exists(fl)) return(NULL)
      fl
    }, format = "file"),
    tar_target(mse_1_mod_pymc_barts_3b_, {
      if(is.null(mse_1_mod_pymc_barts_3b_file_)) return(NULL)
      read_csv(file = mse_1_mod_pymc_barts_3b_file_) 
    }, cue = tar_cue(mode = "always")),

    tar_target(mse_2_mod_pymc_barts_3b_file_, {
      fl <- paste0(
        missing_years_global_parameters_$data_path,
        "modelled/pymc_bart_3b_mse_2.csv"
      )
      if (!file.exists(fl)) return(NULL)
      fl
    }, format = "file"),
    tar_target(mse_2_mod_pymc_barts_3b_, {
      if(is.null(mse_2_mod_pymc_barts_3b_file_)) return(NULL)
      read_csv(file = mse_2_mod_pymc_barts_3b_file_) 
    }, cue = tar_cue(mode = "always")),

    tar_target(comparisons_3_, {
      data_path <- missing_years_global_parameters_$data_path
      mse_1_mod_glmmTMB_3 <- mse_1_mod_glmmTMB_3_
      mse_2_mod_glmmTMB_3 <- mse_2_mod_glmmTMB_3_
      mse_1_mod_glmmTMB_3b <- mse_1_mod_glmmTMB_3b_
      mse_2_mod_glmmTMB_3b <- mse_2_mod_glmmTMB_3b_

      mse_1_mod_brms_3 <- mse_1_mod_brms_3_
      mse_2_mod_brms_3 <- mse_2_mod_brms_3_
      mse_1_mod_brms_3b <- mse_1_mod_brms_3b_
      mse_2_mod_brms_3b <- mse_2_mod_brms_3b_

      mse_1_mod_stan_3 <- mse_1_mod_stan_3_
      mse_2_mod_stan_3 <- mse_2_mod_stan_3_

      mse_1_mod_gbm_3 <- mse_1_mod_gbm_3_
      mse_2_mod_gbm_3 <- mse_2_mod_gbm_3_
      mse_1_mod_gbm_3b <- mse_1_mod_gbm_3b_
      mse_2_mod_gbm_3b <- mse_2_mod_gbm_3b_

      mse_1_mod_dbarts_3 <- mse_1_mod_dbarts_3_
      mse_2_mod_dbarts_3 <- mse_2_mod_dbarts_3_
      mse_1_mod_dbarts_3b <- mse_1_mod_dbarts_3b_
      mse_2_mod_dbarts_3b <- mse_2_mod_dbarts_3b_

      mse_1_mod_xgboost_3b <- mse_1_mod_xgboost_3b_
      mse_2_mod_xgboost_3b <- mse_2_mod_xgboost_3b_

      mse_1_mod_pymc_barts_3b <- mse_1_mod_pymc_barts_3b_
      mse_2_mod_pymc_barts_3b <- mse_2_mod_pymc_barts_3b_
      ## ---- comparisons_3
      mse_1_mod_pymc_barts_3b <- read_csv(
        file = paste0(data_path, "modelled/pymc_bart_3b_mse_1.csv")
      )   
      mse_2_mod_pymc_barts_3b <- read_csv(
        file = paste0(data_path, "modelled/pymc_bart_3b_mse_2.csv")
      )   
      comparisons_3 <- bind_rows(
        mse_1_mod_glmmTMB_3,
        mse_2_mod_glmmTMB_3,
        mse_1_mod_glmmTMB_3b,
        mse_2_mod_glmmTMB_3b,

        mse_1_mod_brms_3,
        mse_2_mod_brms_3,
        mse_1_mod_brms_3b,
        mse_2_mod_brms_3b,
        
        mse_1_mod_stan_3,
        mse_2_mod_stan_3,
        
        mse_1_mod_gbm_3,
        mse_2_mod_gbm_3,
        mse_1_mod_gbm_3b,
        mse_2_mod_gbm_3b,

        mse_1_mod_dbarts_3,
        mse_2_mod_dbarts_3,
        mse_1_mod_dbarts_3b,
        mse_2_mod_dbarts_3b,

        mse_1_mod_xgboost_3b,
        mse_2_mod_xgboost_3b,

        mse_1_mod_pymc_barts_3b,
        mse_2_mod_pymc_barts_3b
        ## mse_3_mod_pymc_barts_3b
        )
      saveRDS(comparisons_3,
        file = paste0(data_path, "synthetic/comparisons_3.rds")
      ) 
       
      comps_3 <- 
        comparisons_3 |>
        dplyr::select(-lower, -upper) |> 
        pivot_longer(
          cols = c("mse_mean", "mse_median", "acc_mean", "acc_median"),
          names_to = "metric",
          values_to = "value"
        ) |>
        mutate(model_type = ifelse(model_type == "", "Without covariates", "With covariates")) |>
        separate(metric, into = c("metric", "stat"), sep = "_") |>
        mutate(type = case_when(
          type == 1 ~ "Predicting the sample data",
          type == 2 ~ "Predicting All reefs"
        )) 
      saveRDS(comps_3,
        file = paste0(data_path, "synthetic/comps_3.rds")
      )  
      ## ----end
      comps_3
    }),
    tar_target(comparisons_3_plots_, {
      comps_3 <- comparisons_3_
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- comparisons_3_tab
      g <-
        comps_3 |> 
        filter(stat == "mean", metric == "mse") |>
        ggplot(aes(x = value, y = model)) +
        geom_segment(aes(xend = 0, yend = model), color = "black") +
        geom_point() +
        facet_grid(model_type~type, scales = "free_y") +
        scale_y_discrete("") +
        scale_x_continuous("MSE") +
        theme_bw() +
        ggtitle("Mean Square Error of models built on Northern reefs data (25 reefs)")
 
      ggsave(
        filename = paste0(
          fig_path, "mse_3_1.png"
        ),
        g,
        width = 10, height = 6, dpi = 100
      )
      g <- 
        comps_3 |> 
        filter(stat == "mean", metric == "acc") |>
        ggplot(aes(x = value, y = model)) +
        geom_segment(aes(xend = 0, yend = model), color = "black") +
        geom_point() +
        facet_grid(model_type~type, scales = "free_y") +
        scale_y_discrete("") +
        scale_x_continuous("Mean inaccuracy (%)", labels = function(x) sprintf("%0.1f%%", x*100)) +
        theme_bw() +
        ggtitle("Mean accuracy of models built on Northern reefs data (25 reefs)")

      ggsave(
        filename = paste0(
          fig_path, "mse_3_2.png"
        ),
        g,
        width = 10, height = 6, dpi = 100
      )
      ## ----end
      comps_3
    }),
    
    

    ## Temporal gap.  Years 3-4 are missing for all sites ------------
    ## simple ---------------------------------------------------------
    tar_target(mod_simple_4_, {
      benthos_fixed_locs_obs_4 <- missing_years_data_prep_4_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- simple_4
      mod_simple_4 <- benthos_fixed_locs_obs_4 |>
        group_by(Year, Site) |>
        summarise(
          Mean = mean(cover),
          Median = median(cover)
        ) |>
        ungroup() |>
        group_by(Year) |> 
        summarise(
          Mean = mean(Mean),
          Median = median(Median)
        ) 
      saveRDS(mod_simple_4,
        file = paste0(data_path, "synthetic/mod_simple_4.rds")
      ) 
      ## ----end
      mod_simple_4
    }),

    tar_target(mod_glmmTMB_4_sample_data_, {
      benthos_fixed_locs_obs_4 <- missing_years_data_prep_4_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_4 sample data
      benthos_fixed_locs_obs_4 
      benthos_fixed_locs_obs_4 |>
        group_by(Year) |>
        summarise(
          Mean = mean(cover),
          Median = median(cover)
          )
      ## ----end
      benthos_fixed_locs_obs_4
    }),
    tar_target(mod_glmmTMB_4_newdata_4_, {
      newdata_4 <- missing_years_newdata_4_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_4 newdata 4
      newdata_4 
      ## ----end
      newdata_4
    }),
    tar_target(mod_glmmTMB_4_newdata_4b_, {
      newdata_4b <- missing_years_newdata_4b_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_4 newdata 4b
      newdata_4b 
      unique(newdata_4b$Year)
      ## ----end
      newdata_4b
    }),
 
    ## glmmTMB --------------------------------------------------------
    tar_target(mod_glmmTMB_4_, {
      benthos_fixed_locs_obs_4 <- missing_years_data_prep_4_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_4
      mod_glmmTMB_4 <- glmmTMB(cover ~ fYear + (1 | Site) + (1 | Transect),
        data = benthos_fixed_locs_obs_4,
        family = "beta_family"
      )
      saveRDS(mod_glmmTMB_4,
        file = paste0(data_path, "synthetic/mod_glmmTMB_4.rds")
      ) 
      ## ----end
      mod_glmmTMB_4
    }),
    tar_target(dharma_mod_glmmTMB_4_, {
      mod_glmmTMB_4 <- mod_glmmTMB_4_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_4_dharma
      glmmTMB_4_dharma <- mod_glmmTMB_4 |> 
        simulateResiduals(n = 1000)
      saveRDS(glmmTMB_4_dharma,
        file = paste0(data_path, "synthetic/glmmTMB_4_dharma.rds")
      ) 
      ## ----end
      glmmTMB_4_dharma
    }),
    tar_target(dharma_mod_glmmTMB_4_plot_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      glmmTMB_4_dharma <- dharma_mod_glmmTMB_4_$glmmTMB_4_dharma
      ## ---- glmmTMB_4_dharma plot
      glmmTMB_4_dharma <- readRDS(
        file = paste0(data_path, "synthetic/glmmTMB_4_dharma.rds")
      )
      g <- wrap_elements(~testUniformity(glmmTMB_4_dharma)) +
        wrap_elements(~plotResiduals(glmmTMB_4_dharma)) +
        wrap_elements(~testDispersion(glmmTMB_4_dharma))
      ggsave(
        filename = paste0(
          fig_path, "R_dharma_mod_glmmTMB_4.png"
        ),
        g,
        width = 10, height = 4, dpi = 72
      )
      ## ----end
    }),
    tar_target(pred_1_mod_glmmTMB_4_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_4 <- mod_glmmTMB_4_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_4b <- missing_years_newdata_4b_
      mod_simple_4 <- mod_simple_4_ #sampled_simple_raw_means_
      ## ---- glmmTMB_4_pred_1
      newdata <- newdata_4b
      true_sum <- mod_simple_4
      glmmTMB_4_pred <- pred_glmmTMB(mod_glmmTMB_4,
        type = 1,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_4"
      )
      ## ----end
      glmmTMB_4_pred
    }),
    tar_target(mse_1_mod_glmmTMB_4_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_4 <- mod_glmmTMB_4_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_4b <- missing_years_newdata_4b_
      ## ---- glmmTMB_4_mse_1
      glmmTMB_4_mse <- mse_glmmTMB(mod_glmmTMB_4,
        newdata = newdata_4b, type = 1, model_type = ""
      )
      ## ----end
      glmmTMB_4_mse
    }),
    tar_target(pred_2_mod_glmmTMB_4_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_4 <- mod_glmmTMB_4_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata <- site_replacements_newdata_
      benthos_reefs_temporal_summary_0 <- read_all_temporal_summary_
      ## ---- glmmTMB_4_pred_2
      newdata <- newdata
      true_sum <- benthos_reefs_temporal_summary_0
      glmmTMB_4_pred <- pred_glmmTMB(mod_glmmTMB_4,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_4"
      )
      ## ----end
      glmmTMB_4_pred
    }),
    tar_target(mse_2_mod_glmmTMB_4_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_4 <- mod_glmmTMB_4_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      ## ---- glmmTMB_4_mse_2
      glmmTMB_4_mse <- mse_glmmTMB(mod_glmmTMB_4,
        newdata = newdata_b, type = 2, model_type = ""
      )
      ## ----end
      glmmTMB_4_mse
    }),

    ## glmmTMB + covariates -------------------------------------------
    tar_target(mod_glmmTMB_4b_, {
      benthos_fixed_locs_obs_4 <- missing_years_data_prep_4_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_4b
      mod_glmmTMB_4b <- glmmTMB(cover ~ fYear + CYC + DHW + OTHER + Latitude + Longitude +
                                   (1 | Site) + (1 | Transect),
        data = benthos_fixed_locs_obs_4,
        family = "beta_family"
      )
      saveRDS(mod_glmmTMB_4b,
        file = paste0(data_path, "synthetic/mod_glmmTMB_4b.rds")
      ) 
      ## ----end
      mod_glmmTMB_4b
    }),
    tar_target(dharma_mod_glmmTMB_4b_, {
      DHARMa_glmmTMB <- DHARMa_glmmTMB_
      mod_glmmTMB_4b <- mod_glmmTMB_4b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- glmmTMB_4b_dharma
      glmmTMB_4b_dharma <- DHARMa_glmmTMB(mod_glmmTMB_4b,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_4b")
      ## ----end
      glmmTMB_4b_dharma
    }),
    tar_target(pred_1_mod_glmmTMB_4b_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_4b <- mod_glmmTMB_4b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_4b <- missing_years_newdata_4b_
      mod_simple_4 <- mod_simple_4_ #sampled_simple_raw_means_
      ## ---- glmmTMB_4b_pred_1
      true_sum <- mod_simple_4
      newdata <- newdata_4b
      glmmTMB_4b_pred <- pred_glmmTMB(mod_glmmTMB_4b,
        type = 1,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_4"
      )
      ## ----end
      glmmTMB_4b_pred
    }),
    tar_target(mse_1_mod_glmmTMB_4b_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_4 <- mod_glmmTMB_4b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_4b <- missing_years_newdata_4b_
      ## ---- glmmTMB_4b_mse_1
      glmmTMB_4_mse <- mse_glmmTMB(mod_glmmTMB_4,
        newdata = newdata_4b, type = 1, model_type = "covariates"
      )
      ## ----end
      glmmTMB_4_mse
    }),
    tar_target(pred_2_mod_glmmTMB_4b_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_4b <- mod_glmmTMB_4b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      benthos_reefs_temporal_summary_0 <- read_all_temporal_summary_
      ## ---- glmmTMB_4b_pred_2
      newdata <- newdata_b
      true_sum <- benthos_reefs_temporal_summary_0
      glmmTMB_4b_pred <- pred_glmmTMB(mod_glmmTMB_4b,
        type = 2,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_4"
      )
      ## ----end
      glmmTMB_4b_pred
    }),
    tar_target(mse_2_mod_glmmTMB_4b_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_4b <- mod_glmmTMB_4b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      ## ---- glmmTMB_4b_mse_2
      glmmTMB_4b_mse <- mse_glmmTMB(mod_glmmTMB_4b,
        newdata = newdata_b, type = 2, model_type = "covariates"
      )
      ## ----end
      glmmTMB_4b_mse
    }),
    
    ## brms -----------------------------------------------------------
    tar_target(mod_brms_4_, {
      benthos_fixed_locs_obs_4 <- missing_years_data_prep_4_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- brms_pre_4
      benthos_fixed_locs_obs_4 |>
        group_by(fYear) |>
        summarise(mean_abundance = qlogis(mean(cover)),
          median_abundance = qlogis(median(cover)),
          sd_abundance = sd(qlogis(cover)) ,
          mad_abundance = mad(qlogis(cover)) 
        ) 
      priors <- c(
        prior(normal(0, 0.5), class = "b"),
        prior(normal(-1.5, 0.5), class = "Intercept"),
        prior(student_t(3, 0, 0.5), class = "sd")
      )
      mod_form <- bf(cover ~ fYear + (1 | Site) + (1 | Transect),
        family = "Beta"
      )
      ## ----end
      ## ---- brms_4
      mod_brms_4 <- brm(mod_form,
        data = benthos_fixed_locs_obs_4,
        iter = 5000,
        warmup = 1000,
        chains = 3,
        cores = 3,
        prior = priors,
        thin =  5,
        control = list(adapt_delta = 0.99),
        backend = "cmdstanr"
      )
      saveRDS(mod_brms_4,
        file = paste0(data_path, "synthetic/mod_brms_4.rds")
      ) 
      ## ----end
      mod_brms_4
    }),
    tar_target(brms_trace_4_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_brms_4 <- mod_brms_4_
      ## ---- brms_trace_4
      mod_brms_4 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_4.rds")
      )
      vars <- mod_brms_4 |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_4$fit |> stan_trace(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_trace_4.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ac_4_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_brms_4 <- mod_brms_4_
      ## ---- brms_ac_4
      mod_brms_4 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_4.rds")
      )
      vars <- mod_brms_4 |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_4$fit |> stan_ac(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ac_4.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_rhat_4_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_brms_4 <- mod_brms_4_
      ## ---- brms_rhat_4
      mod_brms_4 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_4.rds")
      )
      g <- mod_brms_4$fit |> stan_rhat()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_rhat_4.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ess_4_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_brms_4 <- mod_brms_4_
      ## ---- brms_ess_4
      mod_brms_4 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_4.rds")
      )
      g <- mod_brms_4$fit |> stan_ess()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ess_4.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ppc_4_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_brms_4 <- mod_brms_4_
      ## ---- brms_ppc_4
      mod_brms_4 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_4.rds")
      )
      g <- mod_brms_4 |> pp_check( type='dens_overlay', ndraws=100)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ppc_4.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(pred_1_mod_brms_4_, {
      pred_brms <- pred_brms_
      mod_brms_4 <- mod_brms_4_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_4b <- missing_years_newdata_4b_
      mod_simple_4 <- mod_simple_4_ #sampled_simple_raw_means_
      ## ---- brms_4_pred_1
      newdata <- newdata_4b
      true_sum <- mod_simple_4
      brms_4_pred <- pred_brms(mod_brms_4,
        type = 1,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_4"
      )
      ## ----end
      brms_4_pred
    }),
    tar_target(mse_1_mod_brms_4_, {
      mse_brms <- mse_brms_
      mod_brms_4 <- mod_brms_4_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_4b <- missing_years_newdata_4b_
      ## ---- brms_4_mse_1
      brms_4_mse <- mse_brms(mod_brms_4,
        newdata = newdata_4b, type = 1, model_type = ""
      )
      ## ----end
      brms_4_mse
    }),
    tar_target(pred_2_mod_brms_4_, {
      pred_brms <- pred_brms_
      mod_brms_4 <- mod_brms_4_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata <- site_replacements_newdata_
      benthos_reefs_temporal_summary_0 <- read_all_temporal_summary_
      ## ---- brms_4_pred_2
      newdata <- newdata
      true_sum <- benthos_reefs_temporal_summary_0
      brms_4_pred <- pred_brms(mod_brms_4,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_4"
      )
      ## ----end
      brms_4_pred
    }),
    tar_target(mse_2_mod_brms_4_, {
      mse_brms <- mse_brms_
      mod_brms_4 <- mod_brms_4_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      ## ---- brms_4_mse_2
      brms_4_mse <- mse_brms(mod_brms_4,
        newdata = newdata_b, type = 2, model_type = ""
      )
      ## ----end
      brms_4_mse
    }),
    
    ## brms + covariates ----------------------------------------------
    tar_target(mod_brms_4b_, {
      benthos_fixed_locs_obs_4 <- missing_years_data_prep_4_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- brms_pre_4b
      benthos_fixed_locs_obs_4 |>
        group_by(fYear) |>
        summarise(mean_abundance = qlogis(mean(cover)),
          median_abundance = qlogis(median(cover)),
          sd_abundance = sd(qlogis(cover)) ,
          mad_abundance = mad(qlogis(cover)) 
        ) 
      priors <- c(
        prior(normal(0, 0.5), class = "b"),
        prior(normal(-1.5, 0.5), class = "Intercept"),
        prior(student_t(3, 0, 0.5), class = "sd")
      )
      mod_form <- bf(cover ~ fYear + scale(CYC) + scale(DHW) +
                       scale(OTHER) + scale(Latitude) + scale(Longitude) +
                       (1 | Site) + (1 | Transect),
        family = "Beta"
      )
      ## ----end
      ## ---- brms_4b
      mod_brms_4b <- brm(mod_form,
        data = benthos_fixed_locs_obs_4,
        iter = 5000,
        warmup = 1000,
        chains = 3,
        cores = 3,
        prior = priors,
        thin =  5,
        control = list(adapt_delta = 0.99),
        backend = "cmdstanr"
      )
      saveRDS(mod_brms_4b,
        file = paste0(data_path, "synthetic/mod_brms_4b.rds")
      ) 
      ## ----end
      mod_brms_4b
    }),
    tar_target(brms_trace_4b_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_brms_4b <- mod_brms_4b_
      ## ---- brms_trace_4b
      mod_brms_4b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_4b.rds")
      )
      vars <- mod_brms_4b |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_4b$fit |> stan_trace(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_trace_4b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ac_4b_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_brms_4b <- mod_brms_4b_
      ## ---- brms_ac_4b
      mod_brms_4b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_4b.rds")
      )
      vars <- mod_brms_4b |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_4b$fit |> stan_ac(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ac_4b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_rhat_4b_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_brms_4b <- mod_brms_4b_
      ## ---- brms_rhat_4b
      mod_brms_4b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_4b.rds")
      )
      g <- mod_brms_4b$fit |> stan_rhat()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_rhat_4b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ess_4b_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_brms_4b <- mod_brms_4b_
      ## ---- brms_ess_4b
      mod_brms_4b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_4b.rds")
      )
      g <- mod_brms_4b$fit |> stan_ess()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ess_4b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ppc_4b_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_brms_4b <- mod_brms_4b_
      ## ---- brms_ppc_4b
      mod_brms_4b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_4b.rds")
      )
      g <- mod_brms_4b |> pp_check( type='dens_overlay', ndraws=100)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ppc_4b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(pred_1_mod_brms_4b_, {
      pred_brms <- pred_brms_
      mod_brms_4b <- mod_brms_4b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_4b <- missing_years_newdata_4b_
      mod_simple_4 <- mod_simple_4_ #sampled_simple_raw_means_
      ## ---- brms_4b_pred_1
      newdata <- newdata_4b
      true_sum <- mod_simple_4
      brms_4b_pred <- pred_brms(mod_brms_4b,
        type = 1,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_4"
      )
      ## ----end
      brms_4b_pred
    }),
    tar_target(mse_1_mod_brms_4b_, {
      mse_brms <- mse_brms_
      mod_brms_4b <- mod_brms_4b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_4b <- missing_years_newdata_4b_
      ## ---- brms_4b_mse_1
      brms_4b_mse <- mse_brms(mod_brms_4b,
        newdata = newdata_4b, type = 1, model_type = "covariates"
      )
      ## ----end
      brms_4b_mse
    }),
    tar_target(pred_2_mod_brms_4b_, {
      pred_brms <- pred_brms_
      mod_brms_4b <- mod_brms_4b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      benthos_reefs_temporal_summary_0 <- read_all_temporal_summary_
      ## ---- brms_4b_pred_2
      newdata <- newdata_b
      true_sum <- benthos_reefs_temporal_summary_0
      brms_4b_pred <- pred_brms(mod_brms_4b,
        type = 2,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_4"
      )
      ## ----end
      brms_4b_pred
    }),
    tar_target(mse_2_mod_brms_4b_, {
      mse_brms <- mse_brms_
      mod_brms_4b <- mod_brms_4b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      ## ---- brms_4b_mse_2
      brms_4b_mse <- mse_brms(mod_brms_4b,
        newdata = newdata_b, type = 2, model_type = "covariates"
      )
      ## ----end
      brms_4b_mse
    }),
    
    ## stan -----------------------------------------------------------
    tar_target(mod_stan_4_, {
      source("model_functions.R")
      benthos_fixed_locs_obs_4 <- missing_years_data_prep_4_
      data_path <- missing_years_global_parameters_$data_path
      site_extra_functions_
      ## ---- stan_pre_4
      benthos_fixed_locs_obs_4 <-
        benthos_fixed_locs_obs_4 |>
        mutate(
          cYear = fYear,
          grid_id = factor(Reef),
          cSite = factor(interaction(Reef, Site)),
          cReplicate = ifelse(is.na(Transect),
            interaction(Site, Year),
            interaction(cSite, Transect)),
          Cover = cover,
          area = 1,
          sum = 1
        )
      saveRDS(benthos_fixed_locs_obs_4,
        file = paste0(data_path, "synthetic/saveRDS(benthos_fixed_locs_obs_4_forstan.rds")
      ) 
      stan_data <- prepare_data_for_stan(benthos_fixed_locs_obs_4, yrs = 1:12)
      model_stan <- cmdstanr::cmdstan_model(stan_file = "model2.stan")
      ## ----end
      ## ---- stan_4
      mod_stan_4 <- model_stan$sample(
        data = stan_data,
        seed = 123,
        iter_sampling = 5000,
        iter_warmup = 1000,
        thin = 5,
        chains = 3,
        parallel_chains = 3,
        adapt_delta = 0.99,
        output_dir = paste0(data_path, "synthetic/"),
      )
      saveRDS(mod_stan_4,
        file = paste0(data_path, "synthetic/mod_stan_4.rds")
      ) 
      ## ----end
      mod_stan_4
    }),
    tar_target(stan_trace_4_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_stan_4 <- mod_stan_4_
      ## ---- stan_trace_4
      mod_stan_4 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_4.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_4$draws(variables = c("beta", "phi", "sd_1", "sd_2", "sd_3")) |>
        mcmc_trace() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_trace_4.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ac_4_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_stan_4 <- mod_stan_4_
      ## ---- stan_ac_4
      mod_stan_4 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_4.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_4$draws(variables = c("beta", "phi", "sd_1", "sd_2", "sd_3")) |>
        mcmc_acf() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ac_4.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_rhat_4_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_stan_4 <- mod_stan_4_
      ## ---- stan_rhat_4
      mod_stan_4 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_4.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_4 |> bayesplot::rhat() |> 
        mcmc_rhat_hist() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_rhat_4.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ess_4_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_stan_4 <- mod_stan_4_
      ## ---- stan_ess_4
      mod_stan_4 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_4.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_4 |> bayesplot::neff_ratio() |> 
        mcmc_neff_hist() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ess_4.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ppc_4_, {
      benthos_fixed_locs_obs_4 <- missing_years_data_prep_4_
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_stan_4 <- mod_stan_4_
      ## ---- stan_ppc_4
      mod_stan_4 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_4.rds")
      )
      g <- 
        bayesplot::pp_check(
          benthos_fixed_locs_obs_4$cover,
          mod_stan_4$draws("ypred", format = "matrix")[1:100, ],
          ppc_dens_overlay
        ) +
        theme_classic()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ppc_4.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(pred_1_mod_stan_4_, {
      pred_stan <- pred_stan_
      mod_stan_4 <- mod_stan_4_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_4 <- missing_years_newdata_4_
      mod_simple_4 <- mod_simple_4_ #sampled_simple_raw_means_
      ## ---- stan_4_pred_1
      newdata <- newdata_4
      true_sum <- mod_simple_4
      stan_4_pred <- pred_stan(mod_stan_4,
        type = 1,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_stan_4"
      )
      ## ----end
      stan_4_pred
    }),
    tar_target(mse_1_mod_stan_4_, {
      mse_stan <- mse_stan_
      mod_stan_4 <- mod_stan_4_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_4b <- missing_years_newdata_4b_
      ## ---- stan_4_mse_1
      stan_4_mse <- mse_stan(mod_stan_4,
        newdata = newdata_4b, type = 1, model_type = ""
      )
      ## ----end
      stan_4_mse
    }),
    tar_target(pred_2_mod_stan_4_, {
      pred_stan <- pred_stan_
      mod_stan_4 <- mod_stan_4_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata <- site_replacements_newdata_
      newdata_4 <- missing_years_newdata_4_
      benthos_reefs_temporal_summary_0 <- read_all_temporal_summary_
      ## ---- stan_4_pred_2
      newdata <- newdata
      true_sum <- benthos_reefs_temporal_summary_0
      stan_4_pred <- pred_stan(mod_stan_4,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_stan_4",
        newdata_sampled = newdata_4
      )
      ## ----end
      stan_4_pred
    }),
    tar_target(mse_2_mod_stan_4_, {
      mse_stan <- mse_stan_
      mod_stan_4 <- mod_stan_4_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      newdata_4b <- missing_years_newdata_4b_
      ## ---- stan_4_mse_2
      newdata_4b <- newdata_b |>
        right_join(newdata_4b)
      stan_4_mse <- mse_stan(mod_stan_4,
        newdata = newdata_4b, type = 2, model_type = ""
      )
      stan_4b_mse_2 <- mse_stan(mod_stan_4,
        newdata = newdata_b, type = 2, model_type = "",
        sample_yrs_only = FALSE
      ) |> mutate(sample_yrs_only = FALSE)
      stan_4_mse <- bind_rows(stan_4_mse, stan_4b_mse_2)
      ## ----end
      stan_4_mse
    }),
    
    ## gbm ------------------------------------------------------------
    tar_target(mod_gbm_4_, {
      benthos_fixed_locs_obs_4 <- missing_years_data_prep_4_
      data_path <- missing_years_global_parameters_$data_path
      ## ---- gbm_4
      mod_gbm_4 <- gbm(cover ~ fYear,
        data =  benthos_fixed_locs_obs_4,
        distribution = "gaussian",
        n.trees = 10000,
        interaction.depth = 5,
        shrinkage = 0.001,
        bag.fraction = 0.5,
        cv.folds = 5, n.cores = 1,
        verbose = TRUE
      )
      saveRDS(mod_gbm_4,
        file = paste0(data_path, "synthetic/mod_gbm_4.rds")
      ) 
      ## ----end
      ## ---- gbm_post_4
      n.trees <- gbm.perf(mod_gbm_4, method = "cv")
      ## ----end
      list(mod_gbm_4 = mod_gbm_4, n.trees = n.trees)
    }),
    tar_target(pdp_gbm_4_, {
      mod_gbm_4 <- mod_gbm_4_$mod_gbm_4
      n.trees <- mod_gbm_4_$n.trees
      newdata_4 <- missing_years_newdata_4_
      data_path <- missing_years_global_parameters_$data_path
      ## ---- gbm_pdp_4
      gbm_4_sum <- newdata_4 |>
        mutate(median = predict(mod_gbm_4, newdata_4, n.trees = n.trees, type = "response")) |> 
        mutate(Year = as.numeric(as.character(fYear))) |>
        mutate(type = "gbm")
      saveRDS(gbm_4_sum,
        file = paste0(data_path, "synthetic/gbm_4_sum.rds")
      ) 
      ## ----end
      gbm_4_sum
    }),
    tar_target(pred_1_mod_gbm_4_, {
      pred_gbm <- pred_gbm_
      mod_gbm_4 <- mod_gbm_4_$mod_gbm_4
      n.trees <- mod_gbm_4_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_4 <- missing_years_newdata_4_
      mod_simple_4 <- mod_simple_4_ #sampled_simple_raw_means_
      ## ---- gbm_4_pred_1
      newdata <- newdata_4
      true_sum <- mod_simple_4
      gbm_4_pred <- pred_gbm(mod_gbm_4,
        n.trees = n.trees,
        type = 1,
        model_type = "",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_4"
      )
      ## ----end
      gbm_4_pred
    }),
    tar_target(mse_1_mod_gbm_4_, {
      mse_gbm <- mse_gbm_
      mod_gbm_4 <- mod_gbm_4_$mod_gbm_4
      n.trees <- mod_gbm_4_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_4b <- missing_years_newdata_4b_
      ## ---- gbm_4_mse_1
      gbm_4_mse <- mse_gbm(mod_gbm_4,
        n.trees = n.trees,
        newdata = newdata_4b, type = 1, model_type = ""
      )
      ## ----end
      gbm_4_mse
    }),
    tar_target(pred_2_mod_gbm_4_, {
      pred_gbm <- pred_gbm_
      mod_gbm_4 <- mod_gbm_4_$mod_gbm_4
      n.trees <- mod_gbm_4_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata <- site_replacements_newdata_
      newdata_4 <- missing_years_newdata_4_
      benthos_reefs_temporal_summary_0 <- read_all_temporal_summary_
      ## ---- gbm_4_pred_2
      newdata <- newdata
      true_sum <- benthos_reefs_temporal_summary_0
      gbm_4_pred <- pred_gbm(mod_gbm_4,
        n.trees =  n.trees,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_4",
        newdata_sampled = newdata_4
      )
      ## ----end
      gbm_4_pred
    }),
    tar_target(mse_2_mod_gbm_4_, {
      mse_gbm <- mse_gbm_
      mod_gbm_4 <- mod_gbm_4_$mod_gbm_4
      n.trees <- mod_gbm_4_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      ## ---- gbm_4_mse_2
      gbm_4_mse <- mse_gbm(mod_gbm_4,
        n.trees = n.trees,
        newdata = newdata_b, type = 2, model_type = ""
      )
      ## ----end
      gbm_4_mse
    }),
    
    ## gbm + covartiates ----------------------------------------------
    tar_target(mod_gbm_4b_, {
      benthos_fixed_locs_obs_4 <- missing_years_data_prep_4_
      data_path <- missing_years_global_parameters_$data_path
      ## ---- gbm_4b
      mod_gbm_4b <- gbm(cover ~ fYear + Latitude + Longitude + CYC + DHW + OTHER,
        data =  benthos_fixed_locs_obs_4,
        distribution = "gaussian",
        var.monotone = c(0, 0, 0, -1, -1, -1),
        n.trees = 10000,
        interaction.depth = 5,
        shrinkage = 0.001,
        bag.fraction = 0.5,
        cv.folds = 5, n.cores = 1,
        verbose = TRUE
      )
      saveRDS(mod_gbm_4b,
        file = paste0(data_path, "synthetic/mod_gbm_4b.rds")
      )
      ## ----end
      ## ---- gbm_post_4b
      n.trees <- gbm.perf(mod_gbm_4b, method = "cv")
      ## ----end
      list(mod_gbm_4b = mod_gbm_4b, n.trees = n.trees)
    }),
    tar_target(pred_1_mod_gbm_4b_, {
      pred_gbm <- pred_gbm_
      mod_gbm_4b <- mod_gbm_4b_$mod_gbm_4b
      n.trees <- mod_gbm_4b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_4b <- missing_years_newdata_4b_
      mod_simple_4 <- mod_simple_4_ #sampled_simple_raw_means_
      ## ---- gbm_4b_pred_1
      newdata <- newdata_4b
      true_sum <- mod_simple_4
      gbm_4b_pred <- pred_gbm(mod_gbm_4b,
        n.trees = n.trees,
        type = 1,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_4"
      )
      ## ----end
      gbm_4b_pred
    }),
    tar_target(mse_1_mod_gbm_4b_, {
      mse_gbm <- mse_gbm_
      mod_gbm_4b <- mod_gbm_4b_$mod_gbm_4b
      n.trees <- mod_gbm_4b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_4b <- missing_years_newdata_4b_
      ## ---- gbm_4b_mse_1
      gbm_4b_mse <- mse_gbm(mod_gbm_4b,
        n.trees = n.trees,
        newdata = newdata_4b, type = 1, model_type = "covariates"
      )
      ## ----end
      gbm_4b_mse
    }),
    tar_target(pred_2_mod_gbm_4b_, {
      pred_gbm <- pred_gbm_
      mod_gbm_4b <- mod_gbm_4b_$mod_gbm_4b
      n.trees <- mod_gbm_4b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      newdata_4 <- missing_years_newdata_4_
      benthos_reefs_temporal_summary_0 <- read_all_temporal_summary_
      ## ---- gbm_4b_pred_2
      newdata <- newdata_b
      true_sum <- benthos_reefs_temporal_summary_0
      gbm_4b_pred <- pred_gbm(mod_gbm_4b,
        n.trees =  n.trees,
        type = 2,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_4",
        newdata_sampled = newdata_4
      )
      ## ----end
      gbm_4b_pred
    }),
    tar_target(mse_2_mod_gbm_4b_, {
      mse_gbm <- mse_gbm_
      mod_gbm_4b <- mod_gbm_4b_$mod_gbm_4b
      n.trees <- mod_gbm_4b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      ## ---- gbm_4b_mse_2
      gbm_4b_mse <- mse_gbm(mod_gbm_4b,
        n.trees = n.trees,
        newdata = newdata_b, type = 2, model_type = "covariates"
      )
      ## ----end
      gbm_4b_mse
    }),
    
    ## dbarts --------------------------------------------------------
    tar_target(mod_dbarts_4_, {
      benthos_fixed_locs_obs_4 <- missing_years_data_prep_4_
      data_path <- missing_years_global_parameters_$data_path
      newdata_4 <- missing_years_newdata_4_
      newdata_4b <- missing_years_newdata_4b_
      newdata_4_a <- missing_years_data_prep_4_
      newdata <- site_replacements_newdata_
      newdata_a <- site_replacements_newdata_b_
      ## ---- dbarts_4
      print(head(benthos_fixed_locs_obs_4))
      mod_dbarts_4 <- bart2(log(cover) ~ fYear,
        data =   benthos_fixed_locs_obs_4,
        keepTrees = TRUE
      )
      saveRDS(mod_dbarts_4,
        file = paste0(data_path, "synthetic/mod_dbarts_4.rds")
      ) 
      ## ----end
      ## Unfortunately, the next part must be in the same tar_target
      ## due to the way dbarts stores pointers - they cannot be stored
      ## ---- dbarts_pred_4
      preds <- predict(mod_dbarts_4, newdata_4, type = "ev") |>
        exp() 
      saveRDS(preds,
        file = paste0(data_path, "synthetic/mod_dbarts_4_preds.rds")
      ) 
      preds_sum <- preds |>
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_4_preds_sum.rds")
      ) 

      newdata.1 <- newdata_4_a |>
        group_by(Year, Reef) |>
        summarise(
          HCC = mean(HCC, na.rm = TRUE),
          Latitude = mean(Latitude),
          Longitude = mean(Longitude),
          CYC = mean(CYC, na.rm = TRUE),
          DHW = mean(DHW, na.rm = TRUE),
          OTHER = mean(OTHER, na.rm = TRUE)
        ) |>
        mutate(fYear = factor(Year)) |> 
        ungroup()
      preds_a <- predict(mod_dbarts_4, newdata.1, type = "ev") |>
        exp() 
      saveRDS(preds_a,
        file = paste0(data_path, "synthetic/mod_dbarts_4_a_preds.rds")
      ) 

      newdata_4 <- newdata_4 |>
        mutate(Year = as.numeric(as.character(fYear))) 

      newdata <- newdata |>
        right_join(newdata_4, by = c("fYear")) |>
        droplevels()
      preds_0 <- predict(mod_dbarts_4, newdata, type = "ev") |>
        exp()
      saveRDS(preds_0,
        file = paste0(data_path, "synthetic/mod_dbarts_0_preds.rds")
      ) 
      preds_0_sum <- preds_0 |>
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_0_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_0_preds_sum.rds")
      ) 
      newdata_0_a <- newdata_a |>
        ungroup() |>
        mutate(fYear = factor(Year)) |> 
        right_join(newdata_4, by = c("fYear")) |> 
      droplevels()
      preds_0_a <- predict(mod_dbarts_4, newdata_0_a, type = "ev") |>
        exp()
      saveRDS(preds_0_a,
        file = paste0(data_path, "synthetic/mod_dbarts_0_a_preds.rds")
      ) 
      ## ----end
      list(
        mod_dbarts_4 = mod_dbarts_4,
        preds = preds,
        preds_sum = preds_sum,
        preds_a = preds_a,
        preds_0 = preds_0,
        preds_0_sum = preds_0_sum,
        preds_0_a = preds_0_a
      )
    }),
    tar_target(pred_1_mod_dbarts_4_, {
      pred_dbarts <- pred_dbarts_
      preds <- mod_dbarts_4_$preds_sum
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_4 <- missing_years_newdata_4_
      mod_simple_4 <- mod_simple_4_ #sampled_simple_raw_means_
      ## ---- dbarts_4_pred_1
      newdata <- newdata_4
      true_sum <- mod_simple_4
      dbarts_4_pred <- pred_dbarts(
        preds,
        type = 1,
        model_type = "",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_4"
      )
      ## ----end
      dbarts_4_pred
    }),
    tar_target(mse_1_mod_dbarts_4_, {
      mse_dbarts <- mse_dbarts_
      preds <- mod_dbarts_4_$preds_a
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      newdata_4b <- missing_years_newdata_4b_
      ## ---- dbarts_4b_mse_1
      dbarts_4_mse <- mse_dbarts(preds,
        newdata = newdata_4b, type = 1, model_type = ""
      )
      ## ----end
      dbarts_4_mse
    }),
    tar_target(pred_2_mod_dbarts_4_, {
      pred_dbarts <- pred_dbarts_
      preds <- mod_dbarts_4_$preds_0_sum
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## newdata <- site_replacements_newdata_
      newdata_4 <- missing_years_newdata_4_
      benthos_reefs_temporal_summary_0 <- read_all_temporal_summary_
      ## ---- dbarts_4_pred_2
      newdata <- newdata_4
      true_sum <- benthos_reefs_temporal_summary_0
      dbarts_4_pred <- pred_dbarts(
        preds,
        type = 2,
        model_type = "",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_4"
      )
      ## ----end
      dbarts_4_pred
    }),
    tar_target(mse_2_mod_dbarts_4_, {
      mse_dbarts <- mse_dbarts_
      preds <- mod_dbarts_4_$preds_a
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_ 
      ## ---- dbarts_4_mse_2
      dbarts_4_mse <- mse_dbarts(preds,
        newdata = newdata_b, type = 2, model_type = ""
      )
      ## ----end
      dbarts_4_mse
    }),
    
    ## dbarts + covariates --------------------------------------------
    tar_target(mod_dbarts_4b_, {
      benthos_fixed_locs_obs_4 <- missing_years_data_prep_4_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_4 <- missing_years_newdata_4_
      newdata_4b <- missing_years_newdata_4b_
      newdata_b <- site_replacements_newdata_b_
      ## ---- dbarts_4b
      print(head(benthos_fixed_locs_obs_4))
      mod_dbarts_4b <- bart2(log(cover) ~ fYear + Latitude + Longitude + CYC + DHW + OTHER,
        data =   benthos_fixed_locs_obs_4,
        keepTrees = TRUE
      )
      saveRDS(mod_dbarts_4b,
        file = paste0(data_path, "synthetic/mod_dbarts_4b.rds")
      ) 
      ## ----end
      ## Unfortunately, the next part must be in the same tar_target
      ## due to the way dbarts stores pointers - they cannot be stored
      ## ---- dbarts_pred_4b
      newdata_4b <- newdata_4b |>
        mutate(fYear = factor(Year))
      preds_4b <- predict(mod_dbarts_4b, newdata_4b, type = "ev") |>
        exp() 
      saveRDS(preds_4b,
        file = paste0(data_path, "synthetic/mod_dbarts_4b_preds_4b.rds")
      ) 
      preds_4b_sum <- preds_4b |> summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_4b_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_4b_preds_4b_sum.rds")
      )
      
      newdata_b <- newdata_b |>
        mutate(fYear = factor(Year)) |>
        right_join(newdata_4b |>
                     dplyr::select(fYear) |>
                     distinct(),
          by = c("fYear", "Year")) |>
        droplevels()
      preds_b <- predict(mod_dbarts_4b, newdata_b, type = "ev") |>
        exp() 
      saveRDS(preds_b,
        file = paste0(data_path, "synthetic/mod_dbarts_4b_preds_b.rds")
      ) 
      preds_b_sum <- preds_b |> 
        summarise_draws(median, HDInterval::hdi) |>
        mutate(fYear = newdata_b$fYear,
               Year = newdata_b$Year)
      saveRDS(preds_b_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_4b_preds_b_sum.rds")
      ) 

      ## ----end
      list(
        mod_dbarts_4 = mod_dbarts_4b,
        preds_4b = preds_4b,
        preds_4b_sum = preds_4b_sum,
        preds_b = preds_b,
        preds_b_sum = preds_b_sum
      )
    }),
    tar_target(pred_1_mod_dbarts_4b_, {
      pred_dbarts <- pred_dbarts_
      preds_4b <- mod_dbarts_4b_$preds_4b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_4b <- missing_years_newdata_4b_
      mod_simple_4 <- mod_simple_4_ #sampled_simple_raw_means_
      ## ---- dbarts_4b_pred_1
      newdata <- newdata_4b
      true_sum <- mod_simple_4
      dbarts_4b_pred <- pred_dbarts(
        preds_4b,
        type = 1,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_4"
      )
      ## ----end
      dbarts_4b_pred
    }),
    tar_target(mse_1_mod_dbarts_4b_, {
      mse_dbarts <- mse_dbarts_
      preds_4b <- mod_dbarts_4b_$preds_4b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_4b <- missing_years_newdata_4b_
      ## ---- dbarts_4b_mse_1
      dbarts_4_mse <- mse_dbarts(preds_4b,
        newdata = newdata_4b, type = 1, model_type = "covariates"
      )
      ## ----end
      dbarts_4_mse
    }),
    tar_target(pred_2_mod_dbarts_4b_, {
      pred_dbarts <- pred_dbarts_
      preds_b <- mod_dbarts_4b_$preds_b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      newdata_4b <- missing_years_newdata_4b_
      benthos_reefs_temporal_summary_0 <- read_all_temporal_summary_
      ## ---- dbarts_4b_pred_2
      newdata <- newdata_b
      true_sum <- benthos_reefs_temporal_summary_0
      newdata_b <- newdata_b |>
        mutate(fYear = factor(Year)) |>
        right_join(newdata_4b |>
                     dplyr::select(Year) |>
                     distinct(),
          by = c("Year")) |>
        droplevels()
      dbarts_4_pred <- pred_dbarts(
        preds_b, 
        type = 2, 
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_4",
        newdata_sampled = newdata_b
      )
      ## ----end
      dbarts_4_pred
    }),
    tar_target(mse_2_mod_dbarts_4b_, {
      mse_dbarts <- mse_dbarts_
      preds_b <- mod_dbarts_4b_$preds_b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      ## ---- dbarts_4b_mse_2
      dbarts_4_mse <- mse_dbarts(preds_b,
        newdata = newdata_b, type = 2, model_type = "covariates"
      )
      ## ----end
      dbarts_4_mse
    }),
    
    ## xgboost + covariates -------------------------------------------
    tar_target(mod_xgboost_4b_prep_, {
      benthos_fixed_locs_obs_4 <- missing_years_data_prep_4_
      data_path <- missing_years_global_parameters_$data_path
      ## ---- xgboost_4_prep
      data_train_4b <- benthos_fixed_locs_obs_4 |>
        dplyr::select(cover, Year, Latitude, Longitude, CYC, DHW, OTHER) 
      ## ----end
      data_train_4b
    }),
    tar_target(mod_xgboost_4b_tune_, {
      data_train_4b <- mod_xgboost_4b_prep_
      data_path <- missing_years_global_parameters_$data_path
      ## ---- xgboost_4_tune
      ## Define the recipe
      tune_recipe <- recipe(cover ~ ., data = data_train_4b) |> 
        step_dummy(all_nominal_predictors())
      ##Define the model
      tune_model <- boost_tree(learn_rate = tune(),
        trees = tune(), 
        min_n = tune(), 
        tree_depth = tune()) |>    # Model type
        set_engine("xgboost") |>   # Model engine
        set_mode("regression")     # Model mode
      ## Define the workflow
      tune_workflow <- workflow() |> 
        add_recipe(tune_recipe) |> 
        add_model(tune_model)
      tune_grid_values <-
        grid_space_filling(learn_rate(),
          trees(),
          tree_depth(),
          min_n(),
          size = 20,
          type = "max_entropy")
      ##Run the hyper parameters tuning
      tuned_results <- tune_grid(tune_workflow,
        resamples = vfold_cv(data_train_4b, v = 5),
        grid = tune_grid_values)
      ## Get best set of parameters
      model_hyperparams <-
        select_best(tuned_results, metric = "rmse") |> 
        select(-".config") |> 
        as_tibble() |>
        mutate(nb_training = nrow(data_train_4b),
          grid_size = nrow(tune_grid_values)) 
      ## ----end
     list(
        tune_recipe = tune_recipe,
        tune_model = tune_model,
        tune_workflow = tune_workflow,
        tune_grid_values = tune_grid_values,
        tuned_results = tuned_results,
        model_hyperparams = model_hyperparams
     ) 
    }),
    tar_target(mod_xgboost_4b_fit_, {
      data_train_4b <- mod_xgboost_4b_prep_
      data_path <- missing_years_global_parameters_$data_path
      tune_recipe = mod_xgboost_4b_tune_$tune_recipe
      tune_model = mod_xgboost_4b_tune_$tune_model
      tune_workflow = mod_xgboost_4b_tune_$tune_workflow
      tune_grid_values = mod_xgboost_4b_tune_$tune_grid_values
      tuned_results = mod_xgboost_4b_tune_$tuned_results
      model_hyperparams = mod_xgboost_4b_tune_$model_hyperparams
      ## ---- xgboost_4_fit
      ## Redefine the model (with hyper parameters)
      tune_model <-
        boost_tree(learn_rate = model_hyperparams$learn_rate,
          trees = model_hyperparams$trees, 
          min_n = model_hyperparams$min_n, 
          tree_depth = model_hyperparams$tree_depth) |> # Model type
        set_engine("xgboost") |> # Model engine
        set_mode("regression") # Model mode
      ## Redefine the workflow
      tune_workflow <- workflow() |>
        add_recipe(tune_recipe) |> 
        add_model(tune_model)
      ## Fit the final model
      final_fitted_4b <- tune_workflow |>
        fit(data_train_4b)
      ## ----end
      final_fitted_4b
    }),
    tar_target(pred_1_mod_xgboost_4b_, {
      pred_xgboost <- pred_xgboost_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_4b <- missing_years_newdata_4b_
      final_fitted_4b <- mod_xgboost_4b_fit_
      mod_simple_4 <- mod_simple_4_ #sampled_simple_raw_means_
      ## ---- xgboost_4b_pred_1
      true_sum <- mod_simple_4
      newdata <- newdata_4b
      xgboost_4b_pred <- pred_xgboost(
        final_fitted_4b,
        type = 1,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_xgboost_4"
      )
      ## ----end
      xgboost_4b_pred
    }),
    tar_target(mse_1_mod_xgboost_4b_, {
      mse_xgboost <- mse_xgboost_
      final_fitted_4b <- mod_xgboost_4b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_4b <- missing_years_newdata_4b_
      ## ---- xgboost_4b_mse_1
      xgboost_4_mse <- mse_xgboost(final_fitted_4b,
        newdata = newdata_4b, type = 1, model_type = "covariates"
      )
      ## ----end
      xgboost_4_mse
    }),
    tar_target(pred_2_mod_xgboost_4b_, {
      pred_xgboost <- pred_xgboost_
      final_fitted_4b <- mod_xgboost_4b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      newdata_4 <- missing_years_newdata_4_
      benthos_reefs_temporal_summary_0 <- read_all_temporal_summary_
      ## ---- xgboost_4b_pred_2
      newdata <- newdata_b
      newdata_4 <- newdata_4 |>
        mutate(Year = as.numeric(as.character(fYear)))
      true_sum <- benthos_reefs_temporal_summary_0
      xgboost_4_pred <- pred_xgboost(final_fitted_4b,
        type = 2,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_xgboost_4",
        newdata_sampled = newdata_4
      )
      ## ----end
      xgboost_4_pred
    }),
    tar_target(mse_2_mod_xgboost_4b_, {
      mse_xgboost <- mse_xgboost_
      final_fitted_4b <- mod_xgboost_4b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_b <- site_replacements_newdata_b_
      ## ---- xgboost_4b_mse_2
      xgboost_4_mse <- mse_xgboost(final_fitted_4b,
        newdata = newdata_b, type = 2, model_type = "covariates"
      )
      ## ----end
      xgboost_4_mse
    }),
    
    ## Comparisons ----------------------------------------------------
    tar_target(mse_1_mod_pymc_barts_4b_file_, {
      fl <- paste0(
        missing_years_global_parameters_$data_path,
        "modelled/pymc_bart_4b_mse_1.csv"
      )
      if (!file.exists(fl)) return(NULL)
      fl
    }, format = "file"),
    tar_target(mse_1_mod_pymc_barts_4b_, {
      if(is.null(mse_1_mod_pymc_barts_4b_file_)) return(NULL)
      read_csv(file = mse_1_mod_pymc_barts_4b_file_) 
    }, cue = tar_cue(mode = "always")),
    tar_target(mse_2_mod_pymc_barts_4b_file_, {
      fl <- paste0(
        missing_years_global_parameters_$data_path,
        "modelled/pymc_bart_4b_mse_2.csv"
      )
      if (!file.exists(fl)) return(NULL)
      fl
    },
      format = "file"),
    tar_target(mse_2_mod_pymc_barts_4b_, {
      if(is.null(mse_2_mod_pymc_barts_4b_file_)) return(NULL)
      read_csv(file = mse_2_mod_pymc_barts_4b_file_) 
    }, cue = tar_cue(mode = "always")),
    tar_target(comparisons_4_, {
      data_path <- missing_years_global_parameters_$data_path
      mse_1_mod_glmmTMB_4 <- mse_1_mod_glmmTMB_4_
      mse_2_mod_glmmTMB_4 <- mse_2_mod_glmmTMB_4_
      mse_1_mod_glmmTMB_4b <- mse_1_mod_glmmTMB_4b_
      mse_2_mod_glmmTMB_4b <- mse_2_mod_glmmTMB_4b_

      mse_1_mod_brms_4 <- mse_1_mod_brms_4_
      mse_2_mod_brms_4 <- mse_2_mod_brms_4_
      mse_1_mod_brms_4b <- mse_1_mod_brms_4b_
      mse_2_mod_brms_4b <- mse_2_mod_brms_4b_

      mse_1_mod_stan_4 <- mse_1_mod_stan_4_
      mse_2_mod_stan_4 <- mse_2_mod_stan_4_

      mse_1_mod_gbm_4 <- mse_1_mod_gbm_4_
      mse_2_mod_gbm_4 <- mse_2_mod_gbm_4_
      mse_1_mod_gbm_4b <- mse_1_mod_gbm_4b_
      mse_2_mod_gbm_4b <- mse_2_mod_gbm_4b_

      mse_1_mod_dbarts_4 <- mse_1_mod_dbarts_4_
      mse_2_mod_dbarts_4 <- mse_2_mod_dbarts_4_
      mse_1_mod_dbarts_4b <- mse_1_mod_dbarts_4b_
      mse_2_mod_dbarts_4b <- mse_2_mod_dbarts_4b_

      mse_1_mod_xgboost_4b <- mse_1_mod_xgboost_4b_
      mse_2_mod_xgboost_4b <- mse_2_mod_xgboost_4b_

      mse_1_mod_pymc_barts_4b <- mse_1_mod_pymc_barts_4b_
      mse_2_mod_pymc_barts_4b <- mse_2_mod_pymc_barts_4b_
      ## ---- comparisons_4
      mse_1_mod_pymc_barts_4b <- read_csv(
        file = paste0(data_path, "modelled/pymc_bart_4b_mse_1.csv")
      )   
      mse_2_mod_pymc_barts_4b <- read_csv(
        file = paste0(data_path, "modelled/pymc_bart_4b_mse_2.csv")
      )   
      comparisons_4 <- bind_rows(
        mse_1_mod_glmmTMB_4,
        mse_2_mod_glmmTMB_4,
        mse_1_mod_glmmTMB_4b,
        mse_2_mod_glmmTMB_4b,

        mse_1_mod_brms_4,
        mse_2_mod_brms_4,
        mse_1_mod_brms_4b,
        mse_2_mod_brms_4b,
        
        mse_1_mod_stan_4,
        mse_2_mod_stan_4,
        
        mse_1_mod_gbm_4,
        mse_2_mod_gbm_4,
        mse_1_mod_gbm_4b,
        mse_2_mod_gbm_4b,

        mse_1_mod_dbarts_4,
        mse_2_mod_dbarts_4,
        mse_1_mod_dbarts_4b,
        mse_2_mod_dbarts_4b,

        mse_1_mod_xgboost_4b,
        mse_2_mod_xgboost_4b,

        mse_1_mod_pymc_barts_4b,
        mse_2_mod_pymc_barts_4b,
        ## mse_4_mod_pymc_barts_4b
        )
      saveRDS(comparisons_4,
        file = paste0(data_path, "synthetic/comparisons_4.rds")
      ) 
       
      comps_4 <- 
        comparisons_4 |>
        dplyr::select(-lower, -upper) |> 
        pivot_longer(
          cols = c("mse_mean", "mse_median", "acc_mean", "acc_median"),
          names_to = "metric",
          values_to = "value"
        ) |>
        mutate(model_type = ifelse(model_type == "", "Without covariates", "With covariates")) |>
        separate(metric, into = c("metric", "stat"), sep = "_") |>
        mutate(type = case_when(
          type == 1 ~ "Predicting the sample data",
          type == 2 ~ "Predicting All reefs"
        )) 
      saveRDS(comps_4,
        file = paste0(data_path, "synthetic/comps_4.rds")
      )  
      ## ----end
      comps_4
    }),
    tar_target(comparisons_4_plots_, {
      comps_4 <- comparisons_4_
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- comparisons_4_tab
      g <-
        comps_4 |> 
        filter(is.na(sample_yrs_only)) |>
        filter(stat == "mean", metric == "mse") |>
        ggplot(aes(x = value, y = model)) +
        ## geom_segment(aes(xend = 0, yend = model, color = sample_yrs_only), position = position_dodge(width = 0.4)) +
        ## geom_point(aes(color = sample_yrs_only), position = position_dodge(width = 0.4)) +
        geom_segment(aes(xend = 0, yend = model), color = "black") +
        geom_point() +
        facet_grid(model_type~type, scales = "free_y") +
        scale_y_discrete("") +
        scale_x_continuous("MSE") +
        theme_bw() +
        ggtitle("Mean Square Error of models built on Northern reefs data (25 reefs)")
 
      ggsave(
        filename = paste0(
          fig_path, "mse_4_1.png"
        ),
        g,
        width = 10, height = 6, dpi = 100
      )
      g <- 
        comps_4 |> 
        filter(is.na(sample_yrs_only)) |>
        filter(stat == "mean", metric == "acc") |>
        ggplot(aes(x = value, y = model)) +
        geom_segment(aes(xend = 0, yend = model), color = "black") +
        geom_point() +
        facet_grid(model_type~type, scales = "free_y") +
        scale_y_discrete("") +
        scale_x_continuous("Mean inaccuracy (%)", labels = function(x) sprintf("%0.1f%%", x*100)) +
        theme_bw() +
        ggtitle("Mean accuracy of models built on Northern reefs data (25 reefs)")

      ggsave(
        filename = paste0(
          fig_path, "mse_4_2.png"
        ),
        g,
        width = 10, height = 6, dpi = 100
      )
      ## ----end
      comps_4
    })
    
    
    
    
    
   
    
    

  )
  return(targets)
}
