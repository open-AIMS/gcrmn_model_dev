source("model_functions.R")

incomplete_spatial <- function() {
  targets <- list(
    tar_target(incomplete_spatial_link,
      lnk <- synthetic_covariates_),

    # Target: Load raw data
    tar_target(
      incomplete_spatial_libraries_,
      {
        ## ---- site replacement libraries
        library(tidyverse) # for data manipulation and visualisation
        library(sf) # for spatial data handling and visualisation
        library(lwgeom)  # for spatial data handling
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
      incomplete_spatial_extra_functions_,
      {
        ## ---- site replacement functions
        source("model_functions.R")
        ## source("helper_functions.R")
        ## ----end
      }
    ),
    tar_target(
      incomplete_spatial_global_parameters_,
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
    
    ## Northern subdomain (All reefs) 
    tar_target(
      read_sampled_reefs_data_10_,
      {
        benthos_fixed_locs_obs_10 <- synthetic_northern_reefs_benthos_reefs_locs_
        data_path <- incomplete_spatial_global_parameters_$data_path
        ## ---- read sampled reefs data 10
        benthos_fixed_locs_obs_10 <- readRDS(
          file = paste0(
            data_path,
            "synthetic/benthos_reefs_pts_northern.rds"
          )
        ) |>
          mutate(
            HCC = plogis(HCC)
          )
        ## ----end
        benthos_fixed_locs_obs_10 
      }
    ),
    tar_target(
      read_sampled_reefs_data_10_plot_1_,
      {
        benthos_fixed_locs_obs_10 <- read_sampled_reefs_data_10_
        data_path <- incomplete_spatial_global_parameters_$data_path
        fig_path <- incomplete_spatial_global_parameters_$fig_path
        ## ---- read sampled reefs data 10 plot
        g <- benthos_fixed_locs_obs_10 |>
          ggplot() +
          geom_sf(aes(fill = HCC, color = HCC)) +
          facet_wrap(~Year) +
          scale_fill_viridis_c() +
          scale_color_viridis_c() +
          theme_bw()
        ggsave(
          filename = paste0(
            fig_path, "R_sampled_reefs_10_plot.png"
          ),
          g,
          width = 8, height = 6, dpi = 72
        )
        ## ----end
      }
    ),
    tar_target(
      read_sampled_reefs_data_10_plot_2_,
      {
        benthos_fixed_locs_obs_10 <- read_sampled_reefs_data_10_
        data_path <- incomplete_spatial_global_parameters_$data_path
        fig_path <- incomplete_spatial_global_parameters_$fig_path
        ## ---- read sampled reefs data 10 plot 2
        benthos_fixed_locs_obs_10 <- benthos_fixed_locs_obs_10 |>
          mutate(
            fYear = as.factor(Year),
            Reef = as.factor(Reef),
            cover = HCC#/100
          )

        benthos_reefs_temporal_summary_10 <- benthos_fixed_locs_obs_10 |>
          st_drop_geometry() |>
          group_by(Year) |>
          summarise(Mean = mean(HCC),
            Median = median(HCC),
            SD = sd(HCC),
            Lower = quantile(HCC, 0.025),
            Upper = quantile(HCC, 0.975))
        saveRDS(benthos_reefs_temporal_summary_10,
          file = paste0(data_path, "synthetic/benthos_reefs_temporal_summary_10.rds")
        )
        g <-
          benthos_reefs_temporal_summary_10 |>
          ggplot() +
          geom_ribbon(aes(x = Year, ymin = Lower, ymax = Upper), alpha = 0.2) +
          geom_line(aes(x = Year, y = Mean, colour = "mean")) +
          geom_line(aes(x = Year, y = Median, colour = "median")) +
          theme_bw()
        ggsave(
          filename = paste0(
            fig_path, "R_all_temporal_summary_plot_10.png"
          ),
          g,
          width = 8, height = 6, dpi = 72
        )
        ## ----end
        benthos_reefs_temporal_summary_10
      }
    ),

    ## Northern subdomain (Sampled reefs) 
    tar_target(
      read_sampled_reefs_data_12_,
      {
      benthos_fixed_locs_obs_12 <- synthetic_northern_sampled_reefs_benthos_reefs_locs_obs_disturb_
        data_path <- incomplete_spatial_global_parameters_$data_path
        ## ---- read sampled reefs data 12
        benthos_fixed_locs_obs_12 <- readRDS(
          file = paste0(
            data_path,
            "synthetic/benthos_fixed_locs_northern_obs_disturb.rds"
          )) |>
          mutate(
            HCC = HCC/100
          ) 
        ## ----end
        benthos_fixed_locs_obs_12 
      }
    ),
    tar_target(
      read_sampled_reefs_data_12_plot_1_,
      {
        benthos_fixed_locs_obs_12 <- read_sampled_reefs_data_12_
        data_path <- incomplete_spatial_global_parameters_$data_path
        fig_path <- incomplete_spatial_global_parameters_$fig_path
        ## ---- read sampled reefs data 12 plot
        g <- benthos_fixed_locs_obs_12 |>
          st_as_sf(coords = c("Longitude", "Latitude")) |>
          ggplot() +
          geom_sf(aes(fill = HCC, color = HCC)) +
          facet_wrap(~Year) +
          scale_fill_viridis_c() +
          scale_color_viridis_c() +
          theme_bw()
        ggsave(
          filename = paste0(
            fig_path, "R_sampled_reefs_12_plot.png"
          ),
          g,
          width = 8, height = 6, dpi = 72
        )
        ## ----end
      }
    ),
    tar_target(
      read_sampled_reefs_data_12_plot_2_,
      {
        benthos_fixed_locs_obs_12 <- read_sampled_reefs_data_12_
        data_path <- incomplete_spatial_global_parameters_$data_path
        fig_path <- incomplete_spatial_global_parameters_$fig_path
        ## ---- read sampled reefs data 12 plot 2
        ## benthos_fixed_locs_obs_12 <- benthos_fixed_locs_obs_12 |>
        ##   mutate(
        ##     fYear = as.factor(Year),
        ##     Reef = as.factor(Reef),
        ##     cover = HCC/100
        ##   )

        benthos_reefs_temporal_summary_12 <- benthos_fixed_locs_obs_12 |>
          st_drop_geometry() |>
          group_by(Year) |>
          summarise(Mean = mean(HCC),
            Median = median(HCC),
            SD = sd(HCC),
            Lower = quantile(HCC, 0.025),
            Upper = quantile(HCC, 0.975))
        saveRDS(benthos_reefs_temporal_summary_12,
          file = paste0(data_path, "synthetic/benthos_reefs_temporal_summary_12.rds")
        )
        g <-
          benthos_reefs_temporal_summary_12 |>
          ggplot() +
          geom_ribbon(aes(x = Year, ymin = Lower, ymax = Upper), alpha = 0.2) +
          geom_line(aes(x = Year, y = Mean, colour = "mean")) +
          geom_line(aes(x = Year, y = Median, colour = "median")) +
          theme_bw()
        ggsave(
          filename = paste0(
            fig_path, "R_all_temporal_summary_plot_12.png"
          ),
          g,
          width = 8, height = 6, dpi = 72
        )
        ## ----end
        benthos_reefs_temporal_summary_12
      }
    ),

    ## Southern subdomain (All reefs) 
    tar_target(
      read_sampled_reefs_data_11_,
      {
        benthos_fixed_locs_obs_11 <- synthetic_southern_reefs_benthos_reefs_locs_
        data_path <- incomplete_spatial_global_parameters_$data_path
        ## ---- read sampled reefs data 11
        benthos_fixed_locs_obs_11 <- readRDS(
          file = paste0(
            data_path,
            "synthetic/benthos_reefs_pts_southern.rds"
          )
        ) |>
          mutate(
            HCC = plogis(HCC)
          )
        ## ----end
        benthos_fixed_locs_obs_11 
      }
    ),
    tar_target(
      read_sampled_reefs_data_11_plot_1_,
      {
        benthos_fixed_locs_obs_11 <- read_sampled_reefs_data_11_
        data_path <- incomplete_spatial_global_parameters_$data_path
        fig_path <- incomplete_spatial_global_parameters_$fig_path
        ## ---- read sampled reefs data 11 plot
        g <- benthos_fixed_locs_obs_11 |>
          ggplot() +
          geom_sf(aes(fill = HCC, color = HCC)) +
          facet_wrap(~Year) +
          scale_fill_viridis_c() +
          scale_color_viridis_c() +
          theme_bw()
        ggsave(
          filename = paste0(
            fig_path, "R_sampled_reefs_11_plot.png"
          ),
          g,
          width = 8, height = 6, dpi = 72
        )
        ## ----end
      }
    ),
    tar_target(
      read_sampled_reefs_data_11_plot_2_,
      {
        benthos_fixed_locs_obs_11 <- read_sampled_reefs_data_11_
        data_path <- incomplete_spatial_global_parameters_$data_path
        fig_path <- incomplete_spatial_global_parameters_$fig_path
        ## ---- read sampled reefs data 11 plot 2
        benthos_fixed_locs_obs_11 <- benthos_fixed_locs_obs_11 |>
          mutate(
            fYear = as.factor(Year),
            Reef = as.factor(Reef),
            cover = HCC
          )

        benthos_reefs_temporal_summary_11 <- benthos_fixed_locs_obs_11 |>
          st_drop_geometry() |>
          group_by(Year) |>
          summarise(Mean = mean(HCC),
            Median = median(HCC),
            SD = sd(HCC),
            Lower = quantile(HCC, 0.025),
            Upper = quantile(HCC, 0.975))
        saveRDS(benthos_reefs_temporal_summary_11,
          file = paste0(data_path, "synthetic/benthos_reefs_temporal_summary_11.rds")
        )
        g <-
          benthos_reefs_temporal_summary_11 |>
          ggplot() +
          geom_ribbon(aes(x = Year, ymin = Lower, ymax = Upper), alpha = 0.2) +
          geom_line(aes(x = Year, y = Mean, colour = "mean")) +
          geom_line(aes(x = Year, y = Median, colour = "median")) +
          theme_bw()
        ggsave(
          filename = paste0(
            fig_path, "R_all_temporal_summary_plot_11.png"
          ),
          g,
          width = 8, height = 6, dpi = 72
        )
        ## ----end
        benthos_reefs_temporal_summary_11
      }
    ),

    ## Southern subdomain (Sampled reefs) 
    tar_target(
      read_sampled_reefs_data_13_,
      {
        data_path <- incomplete_spatial_global_parameters_$data_path
      benthos_fixed_locs_obs_13 <- synthetic_southern_sampled_reefs_benthos_reefs_locs_obs_disturb_
        ## ---- read sampled reefs data 13
        benthos_fixed_locs_obs_13 <- readRDS(
          file = paste0(
            data_path,
            "synthetic/benthos_fixed_locs_southern_obs_disturb.rds"

          )
        ) |>
          mutate(
            HCC = HCC/100
          )
        ## ----end
        benthos_fixed_locs_obs_13 
      }
    ),
    tar_target(
      read_sampled_reefs_data_13_plot_1_,
      {
        benthos_fixed_locs_obs_13 <- read_sampled_reefs_data_13_
        data_path <- incomplete_spatial_global_parameters_$data_path
        fig_path <- incomplete_spatial_global_parameters_$fig_path
        ## ---- read sampled reefs data 13 plot
        g <- benthos_fixed_locs_obs_13 |>
          st_as_sf(coords = c("Longitude", "Latitude")) |>
          ggplot() +
          geom_sf(aes(fill = HCC, color = HCC)) +
          facet_wrap(~Year) +
          scale_fill_viridis_c() +
          scale_color_viridis_c() +
          theme_bw()
        ggsave(
          filename = paste0(
            fig_path, "R_sampled_reefs_13_plot.png"
          ),
          g,
          width = 8, height = 6, dpi = 72
        )
        ## ----end
      }
    ),
    tar_target(
      read_sampled_reefs_data_13_plot_2_,
      {
        benthos_fixed_locs_obs_13 <- read_sampled_reefs_data_13_
        data_path <- incomplete_spatial_global_parameters_$data_path
        fig_path <- incomplete_spatial_global_parameters_$fig_path
        ## ---- read sampled reefs data 13 plot 2
        benthos_reefs_temporal_summary_13 <- benthos_fixed_locs_obs_13 |>
          st_drop_geometry() |>
          group_by(Year) |>
          summarise(Mean = mean(HCC),
            Median = median(HCC),
            SD = sd(HCC),
            Lower = quantile(HCC, 0.025),
            Upper = quantile(HCC, 0.975))
        saveRDS(benthos_reefs_temporal_summary_13,
          file = paste0(data_path, "synthetic/benthos_reefs_temporal_summary_13.rds")
        )

        g <-
          benthos_reefs_temporal_summary_13 |>
          ggplot() +
          geom_ribbon(aes(x = Year, ymin = Lower, ymax = Upper), alpha = 0.2) +
          geom_line(aes(x = Year, y = Mean, colour = "mean")) +
          geom_line(aes(x = Year, y = Median, colour = "median")) +
          theme_bw()
        ggsave(
          filename = paste0(
            fig_path, "R_all_temporal_summary_plot_13.png"
          ),
          g,
          width = 8, height = 6, dpi = 72
        )
        ## ----end
        benthos_reefs_temporal_summary_13
      }
    ),

    ## Western subdomain (All reefs) 
    tar_target(
      read_sampled_reefs_data_14_,
      {
        benthos_fixed_locs_obs_14 <- synthetic_western_reefs_benthos_reefs_locs_
        data_path <- incomplete_spatial_global_parameters_$data_path
        ## ---- read sampled reefs data 14
        benthos_fixed_locs_obs_14 <- readRDS(
          file = paste0(
            data_path,
            "synthetic/benthos_reefs_pts_western.rds"
          )
        ) |>
          mutate(
            HCC = plogis(HCC)
          )
        ## ----end
        benthos_fixed_locs_obs_14 
      }
    ),
    tar_target(
      read_sampled_reefs_data_14_plot_1_,
      {
        benthos_fixed_locs_obs_14 <- read_sampled_reefs_data_14_
        data_path <- incomplete_spatial_global_parameters_$data_path
        fig_path <- incomplete_spatial_global_parameters_$fig_path
        ## ---- read sampled reefs data 14 plot
        g <- benthos_fixed_locs_obs_14 |>
          ggplot() +
          geom_sf(aes(fill = HCC, color = HCC)) +
          facet_wrap(~Year) +
          scale_fill_viridis_c() +
          scale_color_viridis_c() +
          theme_bw()
        ggsave(
          filename = paste0(
            fig_path, "R_sampled_reefs_14_plot.png"
          ),
          g,
          width = 8, height = 6, dpi = 72
        )
        ## ----end
      }
    ),
    tar_target(
      read_sampled_reefs_data_14_plot_2_,
      {
        benthos_fixed_locs_obs_14 <- read_sampled_reefs_data_14_
        data_path <- incomplete_spatial_global_parameters_$data_path
        fig_path <- incomplete_spatial_global_parameters_$fig_path
        ## ---- read sampled reefs data 14 plot 2
        benthos_fixed_locs_obs_14 <- benthos_fixed_locs_obs_14 |>
          mutate(
            fYear = as.factor(Year),
            Reef = as.factor(Reef),
            cover = HCC#/100
          )

        benthos_reefs_temporal_summary_14 <- benthos_fixed_locs_obs_14 |>
          st_drop_geometry() |>
          group_by(Year) |>
          summarise(Mean = mean(HCC),
            Median = median(HCC),
            SD = sd(HCC),
            Lower = quantile(HCC, 0.025),
            Upper = quantile(HCC, 0.975))
        saveRDS(benthos_reefs_temporal_summary_14,
          file = paste0(data_path, "synthetic/benthos_reefs_temporal_summary_14.rds")
        )
        g <-
          benthos_reefs_temporal_summary_14 |>
          ggplot() +
          geom_ribbon(aes(x = Year, ymin = Lower, ymax = Upper), alpha = 0.2) +
          geom_line(aes(x = Year, y = Mean, colour = "mean")) +
          geom_line(aes(x = Year, y = Median, colour = "median")) +
          theme_bw()
        ggsave(
          filename = paste0(
            fig_path, "R_all_temporal_summary_plot_14.png"
          ),
          g,
          width = 8, height = 6, dpi = 72
        )
        ## ----end
        benthos_reefs_temporal_summary_14
      }
    ),

    ## Western subdomain (Sampled reefs) 
    tar_target(
      read_sampled_reefs_data_16_,
      {
      benthos_fixed_locs_obs_16 <- synthetic_western_sampled_reefs_benthos_reefs_locs_obs_disturb_
        data_path <- incomplete_spatial_global_parameters_$data_path
        ## ---- read sampled reefs data 16
        benthos_fixed_locs_obs_16 <- readRDS(
          file = paste0(
            data_path,
            "synthetic/benthos_fixed_locs_western_obs_disturb.rds"
          )) |>
          mutate(
            HCC = HCC/100
          ) 
        ## ----end
        benthos_fixed_locs_obs_16 
      }
    ),
    tar_target(
      read_sampled_reefs_data_16_plot_1_,
      {
        benthos_fixed_locs_obs_16 <- read_sampled_reefs_data_16_
        data_path <- incomplete_spatial_global_parameters_$data_path
        fig_path <- incomplete_spatial_global_parameters_$fig_path
        ## ---- read sampled reefs data 16 plot
        g <- benthos_fixed_locs_obs_16 |>
          st_as_sf(coords = c("Longitude", "Latitude")) |>
          ggplot() +
          geom_sf(aes(fill = HCC, color = HCC)) +
          facet_wrap(~Year) +
          scale_fill_viridis_c() +
          scale_color_viridis_c() +
          theme_bw()
        ggsave(
          filename = paste0(
            fig_path, "R_sampled_reefs_16_plot.png"
          ),
          g,
          width = 8, height = 6, dpi = 72
        )
        ## ----end
      }
    ),
    tar_target(
      read_sampled_reefs_data_16_plot_2_,
      {
        benthos_fixed_locs_obs_16 <- read_sampled_reefs_data_16_
        data_path <- incomplete_spatial_global_parameters_$data_path
        fig_path <- incomplete_spatial_global_parameters_$fig_path
        ## ---- read sampled reefs data 16 plot 2
        ## benthos_fixed_locs_obs_16 <- benthos_fixed_locs_obs_16 |>
        ##   mutate(
        ##     fYear = as.factor(Year),
        ##     Reef = as.factor(Reef),
        ##     cover = HCC/100
        ##   )

        benthos_reefs_temporal_summary_16 <- benthos_fixed_locs_obs_16 |>
          st_drop_geometry() |>
          group_by(Year) |>
          summarise(Mean = mean(HCC),
            Median = median(HCC),
            SD = sd(HCC),
            Lower = quantile(HCC, 0.025),
            Upper = quantile(HCC, 0.975))
        saveRDS(benthos_reefs_temporal_summary_16,
          file = paste0(data_path, "synthetic/benthos_reefs_temporal_summary_16.rds")
        )
        g <-
          benthos_reefs_temporal_summary_16 |>
          ggplot() +
          geom_ribbon(aes(x = Year, ymin = Lower, ymax = Upper), alpha = 0.2) +
          geom_line(aes(x = Year, y = Mean, colour = "mean")) +
          geom_line(aes(x = Year, y = Median, colour = "median")) +
          theme_bw()
        ggsave(
          filename = paste0(
            fig_path, "R_all_temporal_summary_plot_16.png"
          ),
          g,
          width = 8, height = 6, dpi = 72
        )
        ## ----end
        benthos_reefs_temporal_summary_16
      }
    ),
    
    ## Eastern subdomain (All reefs) 
    tar_target(
      read_sampled_reefs_data_15_,
      {
        benthos_fixed_locs_obs_15 <- synthetic_eastern_reefs_benthos_reefs_locs_
        data_path <- incomplete_spatial_global_parameters_$data_path
        ## ---- read sampled reefs data 15
        benthos_fixed_locs_obs_15 <- readRDS(
          file = paste0(
            data_path,
            "synthetic/benthos_reefs_pts_eastern.rds"
          )
        ) |>
          mutate(
            HCC = plogis(HCC)
          )
        ## ----end
        benthos_fixed_locs_obs_15 
      }
    ),
    tar_target(
      read_sampled_reefs_data_15_plot_1_,
      {
        benthos_fixed_locs_obs_15 <- read_sampled_reefs_data_15_
        data_path <- incomplete_spatial_global_parameters_$data_path
        fig_path <- incomplete_spatial_global_parameters_$fig_path
        ## ---- read sampled reefs data 15 plot
        g <- benthos_fixed_locs_obs_15 |>
          ggplot() +
          geom_sf(aes(fill = HCC, color = HCC)) +
          facet_wrap(~Year) +
          scale_fill_viridis_c() +
          scale_color_viridis_c() +
          theme_bw()
        ggsave(
          filename = paste0(
            fig_path, "R_sampled_reefs_15_plot.png"
          ),
          g,
          width = 8, height = 6, dpi = 72
        )
        ## ----end
      }
    ),
    tar_target(
      read_sampled_reefs_data_15_plot_2_,
      {
        benthos_fixed_locs_obs_15 <- read_sampled_reefs_data_15_
        data_path <- incomplete_spatial_global_parameters_$data_path
        fig_path <- incomplete_spatial_global_parameters_$fig_path
        ## ---- read sampled reefs data 15 plot 2
        benthos_fixed_locs_obs_15 <- benthos_fixed_locs_obs_15 |>
          mutate(
            fYear = as.factor(Year),
            Reef = as.factor(Reef),
            cover = HCC#/100
          )

        benthos_reefs_temporal_summary_15 <- benthos_fixed_locs_obs_15 |>
          st_drop_geometry() |>
          group_by(Year) |>
          summarise(Mean = mean(HCC),
            Median = median(HCC),
            SD = sd(HCC),
            Lower = quantile(HCC, 0.025),
            Upper = quantile(HCC, 0.975))
        saveRDS(benthos_reefs_temporal_summary_15,
          file = paste0(data_path, "synthetic/benthos_reefs_temporal_summary_15.rds")
        )
        g <-
          benthos_reefs_temporal_summary_15 |>
          ggplot() +
          geom_ribbon(aes(x = Year, ymin = Lower, ymax = Upper), alpha = 0.2) +
          geom_line(aes(x = Year, y = Mean, colour = "mean")) +
          geom_line(aes(x = Year, y = Median, colour = "median")) +
          theme_bw()
        ggsave(
          filename = paste0(
            fig_path, "R_all_temporal_summary_plot_15.png"
          ),
          g,
          width = 8, height = 6, dpi = 72
        )
        ## ----end
        benthos_reefs_temporal_summary_15
      }
    ),

    ## Eastern subdomain (Sampled reefs) 
    tar_target(
      read_sampled_reefs_data_17_,
      {
      benthos_fixed_locs_obs_17 <- synthetic_eastern_sampled_reefs_benthos_reefs_locs_obs_disturb_
        data_path <- incomplete_spatial_global_parameters_$data_path
        ## ---- read sampled reefs data 17
        benthos_fixed_locs_obs_17 <- readRDS(
          file = paste0(
            data_path,
            "synthetic/benthos_fixed_locs_eastern_obs_disturb.rds"
          )) |>
          mutate(
            HCC = HCC/100
          ) 
        ## ----end
        benthos_fixed_locs_obs_17 
      }
    ),
    tar_target(
      read_sampled_reefs_data_17_plot_1_,
      {
        benthos_fixed_locs_obs_17 <- read_sampled_reefs_data_17_
        data_path <- incomplete_spatial_global_parameters_$data_path
        fig_path <- incomplete_spatial_global_parameters_$fig_path
        ## ---- read sampled reefs data 17 plot
        g <- benthos_fixed_locs_obs_17 |>
          st_as_sf(coords = c("Longitude", "Latitude")) |>
          ggplot() +
          geom_sf(aes(fill = HCC, color = HCC)) +
          facet_wrap(~Year) +
          scale_fill_viridis_c() +
          scale_color_viridis_c() +
          theme_bw()
        ggsave(
          filename = paste0(
            fig_path, "R_sampled_reefs_17_plot.png"
          ),
          g,
          width = 8, height = 6, dpi = 72
        )
        ## ----end
      }
    ),
    tar_target(
      read_sampled_reefs_data_17_plot_2_,
      {
        benthos_fixed_locs_obs_17 <- read_sampled_reefs_data_17_
        data_path <- incomplete_spatial_global_parameters_$data_path
        fig_path <- incomplete_spatial_global_parameters_$fig_path
        ## ---- read sampled reefs data 17 plot 2
        ## benthos_fixed_locs_obs_17 <- benthos_fixed_locs_obs_17 |>
        ##   mutate(
        ##     fYear = as.factor(Year),
        ##     Reef = as.factor(Reef),
        ##     cover = HCC/100
        ##   )

        benthos_reefs_temporal_summary_17 <- benthos_fixed_locs_obs_17 |>
          st_drop_geometry() |>
          group_by(Year) |>
          summarise(Mean = mean(HCC),
            Median = median(HCC),
            SD = sd(HCC),
            Lower = quantile(HCC, 0.025),
            Upper = quantile(HCC, 0.975))
        saveRDS(benthos_reefs_temporal_summary_17,
          file = paste0(data_path, "synthetic/benthos_reefs_temporal_summary_17.rds")
        )
        g <-
          benthos_reefs_temporal_summary_17 |>
          ggplot() +
          geom_ribbon(aes(x = Year, ymin = Lower, ymax = Upper), alpha = 0.2) +
          geom_line(aes(x = Year, y = Mean, colour = "mean")) +
          geom_line(aes(x = Year, y = Median, colour = "median")) +
          theme_bw()
        ggsave(
          filename = paste0(
            fig_path, "R_all_temporal_summary_plot_17.png"
          ),
          g,
          width = 8, height = 6, dpi = 72
        )
        ## ----end
        benthos_reefs_temporal_summary_17
      }
    ),
    
    ## Data preparations ================================================  

    ## Northern subdomain (All reefs) 
    tar_target(incomplete_spatial_data_prep_10_, {
      benthos_fixed_locs_obs_10 <- read_sampled_reefs_data_10_
      ## ---- sampled data prep 10
      benthos_fixed_locs_obs_10 <- benthos_fixed_locs_obs_10 |>
        mutate(
          fYear = as.factor(Year),
          Reef = as.factor(Reef),
          cover = plogis(HCC)
        )
      benthos_fixed_locs_obs_10 
      ## ----end
      benthos_fixed_locs_obs_10 
    }),
    tar_target(incomplete_spatial_newdata_10_, {
      benthos_fixed_locs_obs_10 <- incomplete_spatial_data_prep_10_
      ## ---- newdata 10
      newdata_10 <-
        benthos_fixed_locs_obs_10 |>
        tidyr::expand(fYear)
      newdata_10
      ## ----end
      newdata_10
    }),
    tar_target(incomplete_spatial_newdata_10b_, {
      benthos_fixed_locs_obs_10 <- incomplete_spatial_data_prep_10_
      ## ---- newdata 10b
      newdata_10b <-
        benthos_fixed_locs_obs_10 |>
        mutate(Longitude = st_coordinates(geometry)[, 1],
          Latitude = st_coordinates(geometry)[, 2]) |>
        st_drop_geometry() |>
        group_by(Year, Reef) |>
        summarise(
          HCC = mean(HCC),
          CYC = mean(CYC),
          DHW = mean(DHW),
          OTHER = mean(OTHER),
          Longitude = mean(Longitude),
          Latitude = mean(Latitude)
          )
      newdata_10b
      ## ----end
      newdata_10b
    }),
    
    ## Northern subdomain (Sampled reefs) 
    tar_target(incomplete_spatial_data_prep_12_, {
      benthos_fixed_locs_obs_12 <- read_sampled_reefs_data_12_
      ## ---- sampled data prep 12
      benthos_fixed_locs_obs_12 <- benthos_fixed_locs_obs_12 |>
        mutate(
          fYear = as.factor(Year),
          Reef = as.factor(Reef),
          Site = interaction(Reef, Site),
          Transect = interaction(Site, Transect),
          cover = HCC
        )
      benthos_fixed_locs_obs_12 
      ## ----end
      benthos_fixed_locs_obs_12 
    }),
    tar_target(incomplete_spatial_newdata_12_, {
      benthos_fixed_locs_obs_12 <- incomplete_spatial_data_prep_12_
      ## ---- newdata 12
      newdata_12 <-
        benthos_fixed_locs_obs_12 |>
        tidyr::expand(fYear)
      newdata_12
      ## ----end
      newdata_12
    }),
    tar_target(incomplete_spatial_newdata_12b_, {
      benthos_fixed_locs_obs_12 <- incomplete_spatial_data_prep_12_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- newdata 12b
      newdata_12b <-
        benthos_fixed_locs_obs_12 |> 
        mutate(Longitude = st_coordinates(geometry)[, 1],
          Latitude = st_coordinates(geometry)[, 2]) |>
        st_drop_geometry() |>
        group_by(Year, Reef) |>
        summarise(
          HCC = mean(HCC),
          CYC = mean(CYC),
          DHW = mean(DHW),
          OTHER = mean(OTHER),
          Longitude = mean(Longitude),
          Latitude = mean(Latitude)
          )
      newdata_12b
      saveRDS(newdata_12b,
        file = paste0(
          data_path,
          "synthetic/newdata_12b.rds"
        )
      )
      ## ----end
      newdata_12b
    }),

    ## Southern subdomain (All reefs) 
    tar_target(incomplete_spatial_data_prep_11_, {
      benthos_fixed_locs_obs_11 <- read_sampled_reefs_data_11_
      ## ---- sampled data prep 11
      benthos_fixed_locs_obs_11 <- benthos_fixed_locs_obs_11 |>
        mutate(
          fYear = as.factor(Year),
          Reef = as.factor(Reef),
          cover = plogis(HCC)
        )
      benthos_fixed_locs_obs_11 
      ## ----end
      benthos_fixed_locs_obs_11 
    }),
    tar_target(incomplete_spatial_newdata_11_, {
      benthos_fixed_locs_obs_11 <- incomplete_spatial_data_prep_11_
      ## ---- newdata 11
      newdata_11 <-
        benthos_fixed_locs_obs_11 |>
        tidyr::expand(fYear)
      newdata_11
      ## ----end
      newdata_11
    }),
    tar_target(incomplete_spatial_newdata_11b_, {
      benthos_fixed_locs_obs_11 <- incomplete_spatial_data_prep_11_
      ## ---- newdata 11b
      newdata_11b <-
        benthos_fixed_locs_obs_11 |>
        mutate(Longitude = st_coordinates(geometry)[, 1],
          Latitude = st_coordinates(geometry)[, 2]) |>
        st_drop_geometry() |>
        group_by(Year, Reef) |>
        summarise(
          HCC = mean(HCC),
          CYC = mean(CYC),
          DHW = mean(DHW),
          OTHER = mean(OTHER),
          Longitude = mean(Longitude),
          Latitude = mean(Latitude)
          )
      newdata_11b
      ## ----end
      newdata_11b
    }),
    
    ## Southern subdomain (Sampled reefs) 
    tar_target(incomplete_spatial_data_prep_13_, {
      benthos_fixed_locs_obs_13 <- read_sampled_reefs_data_13_
      ## ---- sampled data prep 13
      benthos_fixed_locs_obs_13 <- benthos_fixed_locs_obs_13 |>
        mutate(
          fYear = as.factor(Year),
          Reef = as.factor(Reef),
          Site = interaction(Reef, Site),
          Transect = interaction(Site, Transect),
          cover = HCC
        )
      benthos_fixed_locs_obs_13 
      ## ----end
      benthos_fixed_locs_obs_13 
    }),
    tar_target(incomplete_spatial_newdata_13_, {
      benthos_fixed_locs_obs_13 <- incomplete_spatial_data_prep_13_
      ## ---- newdata 13
      newdata_13 <-
        benthos_fixed_locs_obs_13 |>
        tidyr::expand(fYear)
      newdata_13
      ## ----end
      newdata_13
    }),
    tar_target(incomplete_spatial_newdata_13b_, {
      benthos_fixed_locs_obs_13 <- incomplete_spatial_data_prep_13_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- newdata 13b
      newdata_13b <-
        benthos_fixed_locs_obs_13 |>
        mutate(Longitude = st_coordinates(geometry)[, 1],
          Latitude = st_coordinates(geometry)[, 2]) |>
        st_drop_geometry() |>
        group_by(Year, Reef) |>
        summarise(
          HCC = mean(HCC),
          CYC = mean(CYC),
          DHW = mean(DHW),
          OTHER = mean(OTHER),
          Longitude = mean(Longitude),
          Latitude = mean(Latitude)
          )
      saveRDS(newdata_13b,
        file = paste0(
          data_path,
          "synthetic/newdata_13b.rds"
        )
      )
      newdata_13b
      ## ----end
      newdata_13b
    }),

    ## Western subdomain (All reefs) 
    tar_target(incomplete_spatial_data_prep_14_, {
      benthos_fixed_locs_obs_14 <- read_sampled_reefs_data_14_
      ## ---- sampled data prep 14
      benthos_fixed_locs_obs_14 <- benthos_fixed_locs_obs_14 |>
        mutate(
          fYear = as.factor(Year),
          Reef = as.factor(Reef),
          cover = plogis(HCC)
        )
      benthos_fixed_locs_obs_14 
      ## ----end
      benthos_fixed_locs_obs_14 
    }),
    tar_target(incomplete_spatial_newdata_14_, {
      benthos_fixed_locs_obs_14 <- incomplete_spatial_data_prep_14_
      ## ---- newdata 14
      newdata_14 <-
        benthos_fixed_locs_obs_14 |>
        tidyr::expand(fYear)
      newdata_14
      ## ----end
      newdata_14
    }),
    tar_target(incomplete_spatial_newdata_14b_, {
      benthos_fixed_locs_obs_14 <- incomplete_spatial_data_prep_14_
      ## ---- newdata 14b
      newdata_14b <-
        benthos_fixed_locs_obs_14 |>
        mutate(Longitude = st_coordinates(geometry)[, 1],
          Latitude = st_coordinates(geometry)[, 2]) |>
        st_drop_geometry() |>
        group_by(Year, Reef) |>
        summarise(
          HCC = mean(HCC),
          CYC = mean(CYC),
          DHW = mean(DHW),
          OTHER = mean(OTHER),
          Longitude = mean(Longitude),
          Latitude = mean(Latitude)
          )
      newdata_14b
      ## ----end
      newdata_14b
    }),

    ## Western subdomain (Sampled reefs) 
    tar_target(incomplete_spatial_data_prep_16_, {
      benthos_fixed_locs_obs_16 <- read_sampled_reefs_data_16_
      ## ---- sampled data prep 16
      benthos_fixed_locs_obs_16 <- benthos_fixed_locs_obs_16 |>
        mutate(
          fYear = as.factor(Year),
          Reef = as.factor(Reef),
          Site = interaction(Reef, Site),
          Transect = interaction(Site, Transect),
          cover = HCC
        )
      benthos_fixed_locs_obs_16 
      ## ----end
      benthos_fixed_locs_obs_16 
    }),
    tar_target(incomplete_spatial_newdata_16_, {
      benthos_fixed_locs_obs_16 <- incomplete_spatial_data_prep_16_
      ## ---- newdata 16
      newdata_16 <-
        benthos_fixed_locs_obs_16 |>
        tidyr::expand(fYear)
      newdata_16
      ## ----end
      newdata_16
    }),
    tar_target(incomplete_spatial_newdata_16b_, {
      benthos_fixed_locs_obs_16 <- incomplete_spatial_data_prep_16_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- newdata 16b
      newdata_16b <-
        benthos_fixed_locs_obs_16 |> 
        mutate(Longitude = st_coordinates(geometry)[, 1],
          Latitude = st_coordinates(geometry)[, 2]) |>
        st_drop_geometry() |>
        group_by(Year, Reef) |>
        summarise(
          HCC = mean(HCC),
          CYC = mean(CYC),
          DHW = mean(DHW),
          OTHER = mean(OTHER),
          Longitude = mean(Longitude),
          Latitude = mean(Latitude)
          )
      newdata_16b
      saveRDS(newdata_16b,
        file = paste0(
          data_path,
          "synthetic/newdata_16b.rds"
        )
      )
      ## ----end
      newdata_16b
    }),

    ## Eastern subdomain (All reefs) 
    tar_target(incomplete_spatial_data_prep_15_, {
      benthos_fixed_locs_obs_15 <- read_sampled_reefs_data_15_
      ## ---- sampled data prep 15
      benthos_fixed_locs_obs_15 <- benthos_fixed_locs_obs_15 |>
        mutate(
          fYear = as.factor(Year),
          Reef = as.factor(Reef),
          cover = plogis(HCC)
        )
      benthos_fixed_locs_obs_15 
      ## ----end
      benthos_fixed_locs_obs_15 
    }),
    tar_target(incomplete_spatial_newdata_15_, {
      benthos_fixed_locs_obs_15 <- incomplete_spatial_data_prep_15_
      ## ---- newdata 15
      newdata_15 <-
        benthos_fixed_locs_obs_15 |>
        tidyr::expand(fYear)
      newdata_15
      ## ----end
      newdata_15
    }),
    tar_target(incomplete_spatial_newdata_15b_, {
      benthos_fixed_locs_obs_15 <- incomplete_spatial_data_prep_15_
      ## ---- newdata 15b
      newdata_15b <-
        benthos_fixed_locs_obs_15 |>
        mutate(Longitude = st_coordinates(geometry)[, 1],
          Latitude = st_coordinates(geometry)[, 2]) |>
        st_drop_geometry() |>
        group_by(Year, Reef) |>
        summarise(
          HCC = mean(HCC),
          CYC = mean(CYC),
          DHW = mean(DHW),
          OTHER = mean(OTHER),
          Longitude = mean(Longitude),
          Latitude = mean(Latitude)
          )
      newdata_15b
      ## ----end
      newdata_15b
    }),

    ## Eastern subdomain (Sampled reefs) 
    tar_target(incomplete_spatial_data_prep_17_, {
      benthos_fixed_locs_obs_17 <- read_sampled_reefs_data_17_
      ## ---- sampled data prep 17
      benthos_fixed_locs_obs_17 <- benthos_fixed_locs_obs_17 |>
        mutate(
          fYear = as.factor(Year),
          Reef = as.factor(Reef),
          Site = interaction(Reef, Site),
          Transect = interaction(Site, Transect),
          cover = HCC
        )
      benthos_fixed_locs_obs_17 
      ## ----end
      benthos_fixed_locs_obs_17 
    }),
    tar_target(incomplete_spatial_newdata_17_, {
      benthos_fixed_locs_obs_17 <- incomplete_spatial_data_prep_17_
      ## ---- newdata 17
      newdata_17 <-
        benthos_fixed_locs_obs_17 |>
        tidyr::expand(fYear)
      newdata_17
      ## ----end
      newdata_17
    }),
    tar_target(incomplete_spatial_newdata_17b_, {
      benthos_fixed_locs_obs_17 <- incomplete_spatial_data_prep_17_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- newdata 17b
      newdata_17b <-
        benthos_fixed_locs_obs_17 |> 
        mutate(Longitude = st_coordinates(geometry)[, 1],
          Latitude = st_coordinates(geometry)[, 2]) |>
        st_drop_geometry() |>
        group_by(Year, Reef) |>
        summarise(
          HCC = mean(HCC),
          CYC = mean(CYC),
          DHW = mean(DHW),
          OTHER = mean(OTHER),
          Longitude = mean(Longitude),
          Latitude = mean(Latitude)
          )
      newdata_17b
      saveRDS(newdata_17b,
        file = paste0(
          data_path,
          "synthetic/newdata_17b.rds"
        )
      )
      ## ----end
      newdata_17b
    }),

    ## Models =========================================================

    ## Northern subdomain (Sampled reefs) 
    ## simple ---------------------------------------------------------
    tar_target(mod_simple_12_, {
      benthos_fixed_locs_obs_12 <- incomplete_spatial_data_prep_12_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- simple_12
      mod_simple_12 <- benthos_fixed_locs_obs_12 |>
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
      mod_simple_12 <-
        benthos_fixed_locs_obs_12 |>
        group_by(Year) |>
        summarise(
          Mean = mean(cover),
          Median = median(cover)
        ) 
      saveRDS(mod_simple_12,
        file = paste0(data_path, "synthetic/mod_simple_12.rds")
      ) 
      ## ----end
      mod_simple_12
    }),

    ## glmmTMB --------------------------------------------------------
    tar_target(mod_glmmTMB_12_, {
      benthos_fixed_locs_obs_12 <- incomplete_spatial_data_prep_12_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_12
      mod_glmmTMB_12 <- glmmTMB(cover ~ fYear + (1 | Site) + (1 | Transect),
        data = benthos_fixed_locs_obs_12,
        family = "beta_family"
      )
      saveRDS(mod_glmmTMB_12,
        file = paste0(data_path, "synthetic/mod_glmmTMB_12.rds")
      ) 
      ## ----end
      mod_glmmTMB_12
    }),
    tar_target(mod_glmmTMB_12_sample_data_, {
      benthos_fixed_locs_obs_12 <- incomplete_spatial_data_prep_12_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_12 sample data
      benthos_fixed_locs_obs_12 
      benthos_fixed_locs_obs_12 |>
        group_by(Year) |>
        summarise(
          Mean = mean(cover),
          Median = median(cover)
          )
      ## ----end
      benthos_fixed_locs_obs_12
    }),
    tar_target(mod_glmmTMB_12_sample_data_summary_10_, {
      benthos_reefs_temporal_summary_10 <- read_sampled_reefs_data_10_plot_2_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_12 sample data summary 10
      benthos_reefs_temporal_summary_10 <- readRDS(
        file = paste0(data_path, "synthetic/benthos_reefs_temporal_summary_10.rds")
      )
      benthos_reefs_temporal_summary_10
      ## ----end
      benthos_reefs_temporal_summary_10
    }),
    tar_target(mod_glmmTMB_12_sample_data_summary_11_, {
      benthos_reefs_temporal_summary_11 <- read_sampled_reefs_data_11_plot_2_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_12 sample data summary 11
      benthos_reefs_temporal_summary_11 <- readRDS(
        file = paste0(data_path, "synthetic/benthos_reefs_temporal_summary_11.rds")
      )
      benthos_reefs_temporal_summary_11
      ## ----end
      benthos_reefs_temporal_summary_11
    }),
    tar_target(mod_glmmTMB_12_newdata_12_, {
      newdata_12 <- incomplete_spatial_newdata_12_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_12 newdata 12
      newdata_12 
      ## ----end
      newdata_12
    }),
    tar_target(mod_glmmTMB_12_newdata_12b_, {
      newdata_12b <- incomplete_spatial_newdata_12b_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_12 newdata 12b
      newdata_12b 
      ## ----end
      newdata_12b
    }),
    tar_target(mod_glmmTMB_12_newdata_13_, {
      newdata_13 <- incomplete_spatial_newdata_13_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_12 newdata 13
      newdata_13 
      ## ----end
      newdata_13
    }),
    tar_target(mod_glmmTMB_12_newdata_13b_, {
      newdata_13b <- incomplete_spatial_newdata_13_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_12 newdata 13b
      newdata_13b 
      ## ----end
      newdata_13b
    }),
    tar_target(mod_glmmTMB_12_newdata_10_, {
      newdata_10 <- incomplete_spatial_newdata_10_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_12 newdata 10
      newdata_10 
      ## ----end
      newdata_10
    }),
    tar_target(mod_glmmTMB_12_newdata_10b_, {
      newdata_10b <- incomplete_spatial_newdata_10b_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_12 newdata 10b
      newdata_10b 
      ## ----end
      newdata_10b
    }),
    tar_target(mod_glmmTMB_12_newdata_11_, {
      newdata_11 <- incomplete_spatial_newdata_11_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_12 newdata 11
      newdata_11 
      ## ----end
      newdata_11
    }),
    tar_target(mod_glmmTMB_12_newdata_11b_, {
      newdata_11b <- incomplete_spatial_newdata_11b_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_12 newdata 11b
      newdata_11b 
      ## ----end
      newdata_11b
    }),
    tar_target(dharma_mod_glmmTMB_12_, {
      DHARMa_glmmTMB <- DHARMa_glmmTMB_
      mod_glmmTMB_12 <- mod_glmmTMB_12_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- glmmTMB_12_dharma
      glmmTMB_12_dharma <- DHARMa_glmmTMB(mod_glmmTMB_12,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_12")
      ## ----end
      glmmTMB_12_dharma
    }),
    tar_target(pred_1_mod_glmmTMB_12_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_12 <- mod_glmmTMB_12_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_12 <- incomplete_spatial_newdata_12_
      mod_simple_12 <- mod_simple_12_ #sampled_simple_raw_means_
      ## ---- glmmTMB_12_pred_1
      newdata <- newdata_12
      true_sum <- mod_simple_12
      glmmTMB_12_pred <- pred_glmmTMB(mod_glmmTMB_12,
        type = 1,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_12"
      )
      ## ----end
      glmmTMB_12_pred
    }),
    tar_target(pred_2_mod_glmmTMB_12_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_12 <- mod_glmmTMB_12_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10 <- incomplete_spatial_newdata_10_
      benthos_reefs_temporal_summary_10 <- read_sampled_reefs_data_10_plot_2_
      ## ---- glmmTMB_12_pred_2
      newdata <- newdata_10
      true_sum <- benthos_reefs_temporal_summary_10
      glmmTMB_12_pred <- pred_glmmTMB(mod_glmmTMB_12,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_12"
      )
      ## ----end
      glmmTMB_12_pred
    }),
    tar_target(pred_3_mod_glmmTMB_12_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_12 <- mod_glmmTMB_12_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11 <- incomplete_spatial_newdata_11_
      benthos_reefs_temporal_summary_11 <- read_sampled_reefs_data_11_plot_2_
      ## ---- glmmTMB_12_pred_3
      newdata <- newdata_11
      true_sum <- benthos_reefs_temporal_summary_11
      glmmTMB_12_pred <- pred_glmmTMB(mod_glmmTMB_12,
        type = 3,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_12"
      )
      ## ----end
      glmmTMB_12_pred
    }),
    tar_target(mse_1_mod_glmmTMB_12_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_12 <- mod_glmmTMB_12_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_12b <- incomplete_spatial_newdata_12b_
      ## ---- glmmTMB_12_mse_1
      glmmTMB_12_mse <- mse_glmmTMB(mod_glmmTMB_12,
        newdata = newdata_12b, type = 1, model_type = ""
      )
      ## ----end
      glmmTMB_12_mse
    }),
    tar_target(mse_2_mod_glmmTMB_12_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_12 <- mod_glmmTMB_12_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      ## ---- glmmTMB_12_mse_2
      glmmTMB_12_mse <- mse_glmmTMB(mod_glmmTMB_12,
        newdata = newdata_10b, type = 2, model_type = ""
      )
      ## ----end
      glmmTMB_12_mse
    }),
    tar_target(mse_3_mod_glmmTMB_12_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_12 <- mod_glmmTMB_12_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      ## ---- glmmTMB_12_mse_3
      glmmTMB_12_mse <- mse_glmmTMB(mod_glmmTMB_12,
        newdata = newdata_11b, type = 3, model_type = ""
      )
      ## ----end
      glmmTMB_12_mse
    }),

    ## glmmTMB + covariates -------------------------------------------
    tar_target(mod_glmmTMB_12b_, {
      benthos_fixed_locs_obs_12 <- incomplete_spatial_data_prep_12_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_12b
      mod_glmmTMB_12b <- glmmTMB(cover ~ fYear + CYC + DHW + OTHER + Latitude + Longitude +
                                   (1 | Site) + (1 | Transect),
        data = benthos_fixed_locs_obs_12,
        family = "beta_family"
      )
      saveRDS(mod_glmmTMB_12b,
        file = paste0(data_path, "synthetic/mod_glmmTMB_12b.rds")
      ) 
      ## ----end
      mod_glmmTMB_12b
    }),
    tar_target(dharma_mod_glmmTMB_12b_, {
      DHARMa_glmmTMB <- DHARMa_glmmTMB_
      mod_glmmTMB_12b <- mod_glmmTMB_12b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- glmmTMB_12b_dharma
      glmmTMB_12b_dharma <- DHARMa_glmmTMB(mod_glmmTMB_12b,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_12b")
      ## ----end
      glmmTMB_12b_dharma
    }),
    tar_target(pred_1_mod_glmmTMB_12b_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_12b <- mod_glmmTMB_12b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_12b <- incomplete_spatial_newdata_12b_
      mod_simple_12 <- mod_simple_12_ #sampled_simple_raw_means_
      ## ---- glmmTMB_12b_pred_1
      true_sum <- mod_simple_12
      newdata <- newdata_12b
      glmmTMB_12b_pred <- pred_glmmTMB(mod_glmmTMB_12b,
        type = 1,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_12"
      )
      ## ----end
      glmmTMB_12b_pred
    }),
    tar_target(pred_2_mod_glmmTMB_12b_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_12 <- mod_glmmTMB_12b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      benthos_reefs_temporal_summary_10 <- read_sampled_reefs_data_10_plot_2_
      ## ---- glmmTMB_12b_pred_2
      newdata <- newdata_10b
      true_sum <- benthos_reefs_temporal_summary_10
      glmmTMB_12b_pred <- pred_glmmTMB(mod_glmmTMB_12,
        type = 2,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_12"
      ) 
      ## ----end
      glmmTMB_12b_pred
    }),
    tar_target(pred_3_mod_glmmTMB_12b_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_12 <- mod_glmmTMB_12b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      benthos_reefs_temporal_summary_11 <- read_sampled_reefs_data_11_plot_2_
      ## ---- glmmTMB_12b_pred_3
      newdata <- newdata_11b
      true_sum <- benthos_reefs_temporal_summary_11
      glmmTMB_12b_pred <- pred_glmmTMB(mod_glmmTMB_12,
        type = 3,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_12"
      )
      ## ----end
      glmmTMB_12b_pred
    }),
    tar_target(mse_1_mod_glmmTMB_12b_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_12 <- mod_glmmTMB_12b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_12b <- incomplete_spatial_newdata_12b_
      ## ---- glmmTMB_12b_mse_1
      glmmTMB_12_mse <- mse_glmmTMB(mod_glmmTMB_12,
        newdata = newdata_12b, type = 1, model_type = "covariates"
      )
      ## ----end
      glmmTMB_12_mse
    }),
    tar_target(mse_2_mod_glmmTMB_12b_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_12 <- mod_glmmTMB_12b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      ## ---- glmmTMB_12b_mse_2
      glmmTMB_12_mse <- mse_glmmTMB(mod_glmmTMB_12,
        newdata = newdata_10b, type = 2, model_type = "covariates"
      )
      ## ----end
      glmmTMB_12_mse
    }),
    tar_target(mse_3_mod_glmmTMB_12b_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_12 <- mod_glmmTMB_12b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      ## ---- glmmTMB_12b_mse_3
      glmmTMB_12_mse <- mse_glmmTMB(mod_glmmTMB_12,
        newdata = newdata_11b, type = 3, model_type = "covariates"
      )
      ## ----end
      glmmTMB_12_mse
    }),
    
    ## brms -----------------------------------------------------------
    tar_target(mod_brms_12_, {
      benthos_fixed_locs_obs_12 <- incomplete_spatial_data_prep_12_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- brms_pre_12
      benthos_fixed_locs_obs_12 |>
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
      ## ---- brms_12
      mod_brms_12 <- brm(mod_form,
        data = benthos_fixed_locs_obs_12,
        iter = 5000,
        warmup = 1000,
        chains = 3,
        cores = 3,
        prior = priors,
        thin =  5,
        control = list(adapt_delta = 0.99),
        backend = "cmdstanr"
      )
      saveRDS(mod_brms_12,
        file = paste0(data_path, "synthetic/mod_brms_12.rds")
      ) 
      ## ----end
      mod_brms_12
    }),
    tar_target(pred_1_mod_brms_12_, {
      pred_brms <- pred_brms_
      mod_brms_12 <- mod_brms_12_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_12 <- incomplete_spatial_newdata_12_
      mod_simple_12 <- mod_simple_12_ #sampled_simple_raw_means_
      ## ---- brms_12_pred_1
      newdata <- newdata_12
      true_sum <- mod_simple_12
      brms_12_pred <- pred_brms(mod_brms_12,
        type = 1,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_12"
      )
      ## ----end
      brms_12_pred
    }),
    tar_target(pred_2_mod_brms_12_, {
      pred_brms <- pred_brms_
      mod_brms_12 <- mod_brms_12_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10 <- incomplete_spatial_newdata_10_
      benthos_reefs_temporal_summary_10 <- read_sampled_reefs_data_10_plot_2_
      ## ---- brms_12_pred_2
      newdata <- newdata_10
      true_sum <- benthos_reefs_temporal_summary_10
      brms_12_pred <- pred_brms(mod_brms_12,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_12"
      )
      ## ----end
      brms_12_pred
    }),
    tar_target(pred_3_mod_brms_12_, {
      pred_brms <- pred_brms_
      mod_brms_12 <- mod_brms_12_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11 <- incomplete_spatial_newdata_11_
      benthos_reefs_temporal_summary_11 <- read_sampled_reefs_data_11_plot_2_
      ## ---- brms_12_pred_3
      newdata <- newdata_11
      true_sum <- benthos_reefs_temporal_summary_11
      brms_12_pred <- pred_brms(mod_brms_12,
        type = 3,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_12"
      )
      ## ----end
      brms_12_pred
    }),
    tar_target(brms_trace_12_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_12 <- mod_brms_12_
      ## ---- brms_trace_12
      mod_brms_12 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_12.rds")
      )
      vars <- mod_brms_12 |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_12$fit |> stan_trace(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_trace_12.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ac_12_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_12 <- mod_brms_12_
      ## ---- brms_ac_12
      mod_brms_12 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_12.rds")
      )
      vars <- mod_brms_12 |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_12$fit |> stan_ac(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ac_12.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_rhat_12_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_12 <- mod_brms_12_
      ## ---- brms_rhat_12
      mod_brms_12 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_12.rds")
      )
      g <- mod_brms_12$fit |> stan_rhat()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_rhat_12.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ess_12_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_12 <- mod_brms_12_
      ## ---- brms_ess_12
      mod_brms_12 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_12.rds")
      )
      g <- mod_brms_12$fit |> stan_ess()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ess_12.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ppc_12_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_12 <- mod_brms_12_
      ## ---- brms_ppc_12
      mod_brms_12 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_12.rds")
      )
      g <- mod_brms_12 |> pp_check( type='dens_overlay', ndraws=100)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ppc_12.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(mse_1_mod_brms_12_, {
      mse_brms <- mse_brms_
      mod_brms_12 <- mod_brms_12_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_12b <- incomplete_spatial_newdata_12b_
      ## ---- brms_12_mse_1
      brms_12_mse <- mse_brms(mod_brms_12,
        newdata = newdata_12b, type = 1, model_type = ""
      )
      ## ----end
      brms_12_mse
    }),
    tar_target(mse_2_mod_brms_12_, {
      mse_brms <- mse_brms_
      mod_brms_12 <- mod_brms_12_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      ## ---- brms_12_mse_2
      brms_12_mse <- mse_brms(mod_brms_12,
        newdata = newdata_10b, type = 2, model_type = ""
      )
      ## ----end
      brms_12_mse
    }),
    tar_target(mse_3_mod_brms_12_, {
      mse_brms <- mse_brms_
      mod_brms_12 <- mod_brms_12_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      ## ---- brms_12_mse_3
      brms_12_mse <- mse_brms(mod_brms_12,
        newdata = newdata_11b, type = 3, model_type = ""
      )
      ## ----end
      brms_12_mse
    }),

    ## brms + covariates ----------------------------------------------
    tar_target(mod_brms_12b_, {
      benthos_fixed_locs_obs_12 <- incomplete_spatial_data_prep_12_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- brms_pre_12b
      benthos_fixed_locs_obs_12 |>
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
      mod_form <- bf(cover ~ fYear + scale(Longitude) + scale(Latitude) +
                       scale(CYC) + scale(DHW) + scale(OTHER) +
                       (1 | Site) + (1 | Transect),
        family = "Beta"
      )
      ## ----end
      ## ---- brms_12
      mod_brms_12b <- brm(mod_form,
        data = benthos_fixed_locs_obs_12,
        iter = 5000,
        warmup = 1000,
        chains = 3,
        cores = 3,
        prior = priors,
        thin =  5,
        control = list(adapt_delta = 0.99),
        backend = "cmdstanr"
      )
      saveRDS(mod_brms_12b,
        file = paste0(data_path, "synthetic/mod_brms_12b.rds")
      ) 
      ## ----end
      mod_brms_12b
    }),
    tar_target(pred_1_mod_brms_12b_, {
      pred_brms <- pred_brms_
      mod_brms_12b <- mod_brms_12b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_12b <- incomplete_spatial_newdata_12b_
      mod_simple_12 <- mod_simple_12_ #sampled_simple_raw_means_
      ## ---- brms_12b_pred_1
      true_sum <- mod_simple_12
      newdata <- newdata_12b
      brms_12b_pred <- pred_brms(mod_brms_12b,
        type = 1,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_12"
      )
      ## ----end
      brms_12b_pred
    }),
    tar_target(pred_2_mod_brms_12b_, {
      pred_brms <- pred_brms_
      mod_brms_12 <- mod_brms_12b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      benthos_reefs_temporal_summary_10 <- read_sampled_reefs_data_10_plot_2_
      ## ---- brms_12b_pred_2
      newdata <- newdata_10b
      true_sum <- benthos_reefs_temporal_summary_10
      brms_12b_pred <- pred_brms(mod_brms_12,
        type = 2,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_12"
      )
      ## ----end
      brms_12b_pred
    }),
    tar_target(pred_3_mod_brms_12b_, {
      pred_brms <- pred_brms_
      mod_brms_12 <- mod_brms_12b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      benthos_reefs_temporal_summary_11 <- read_sampled_reefs_data_11_plot_2_
      ## ---- brms_12b_pred_3
      newdata <- newdata_11b
      true_sum <- benthos_reefs_temporal_summary_11
      brms_12b_pred <- pred_brms(mod_brms_12,
        type = 3,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_12"
      )
      ## ----end
      brms_12b_pred
    }),
    tar_target(brms_trace_12b_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_12b <- mod_brms_12b_
      ## ---- brms_trace_12
      mod_brms_12b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_12b.rds")
      )
      vars <- mod_brms_12b |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_12b$fit |> stan_trace(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_trace_12b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ac_12b_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_12b <- mod_brms_12b_
      ## ---- brms_ac_12b
      mod_brms_12b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_12b.rds")
      )
      vars <- mod_brms_12b |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_12b$fit |> stan_ac(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ac_12b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_rhat_12b_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_12b <- mod_brms_12b_
      ## ---- brms_rhat_12b
      mod_brms_12b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_12b.rds")
      )
      g <- mod_brms_12b$fit |> stan_rhat()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_rhat_12b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ess_12b_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_12b <- mod_brms_12b_
      ## ---- brms_ess_12b
      mod_brms_12b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_12b.rds")
      )
      g <- mod_brms_12b$fit |> stan_ess()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ess_12b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ppc_12b_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_12b <- mod_brms_12b_
      ## ---- brms_ppc_12b
      mod_brms_12b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_12b.rds")
      )
      g <- mod_brms_12b |> pp_check( type='dens_overlay', ndraws=100)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ppc_12b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(mse_1_mod_brms_12b_, {
      mse_brms <- mse_brms_
      mod_brms_12 <- mod_brms_12b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_12b <- incomplete_spatial_newdata_12b_
      ## ---- brms_12b_mse_1
      brms_12_mse <- mse_brms(mod_brms_12,
        newdata = newdata_12b, type = 1, model_type = "covariates"
      )
      ## ----end
      brms_12_mse
    }),
    tar_target(mse_2_mod_brms_12b_, {
      mse_brms <- mse_brms_
      mod_brms_12 <- mod_brms_12b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      ## ---- brms_12b_mse_2
      brms_12_mse <- mse_brms(mod_brms_12,
        newdata = newdata_10b, type = 2, model_type = "covariates"
      )
      ## ----end
      brms_12_mse
    }),
    tar_target(mse_3_mod_brms_12b_, {
      mse_brms <- mse_brms_
      mod_brms_12 <- mod_brms_12b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      ## ---- brms_12b_mse_3
      brms_12_mse <- mse_brms(mod_brms_12,
        newdata = newdata_11b, type = 3, model_type = "covariates"
      )
      ## ----end
      brms_12_mse
    }),

    ## stan -----------------------------------------------------------
    tar_target(mod_stan_12_, {
      source("model_functions.R")
      benthos_fixed_locs_obs_12 <- incomplete_spatial_data_prep_12_
      data_path <- incomplete_spatial_global_parameters_$data_path
      site_extra_functions_
      ## ---- stan_pre_12
      benthos_fixed_locs_obs_12 <-
        benthos_fixed_locs_obs_12 |>
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
      saveRDS(benthos_fixed_locs_obs_12,
        file = paste0(data_path, "synthetic/saveRDS(benthos_fixed_locs_obs_12_forstan.rds")
      ) 
      stan_data <- prepare_data_for_stan(benthos_fixed_locs_obs_12, yrs = NULL)
      model_stan <- cmdstanr::cmdstan_model(stan_file = "model1.stan")
      ## ----end
      ## ---- stan_12
      mod_stan_12 <- model_stan$sample(
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
      saveRDS(mod_stan_12,
        file = paste0(data_path, "synthetic/mod_stan_12.rds")
      ) 
      ## ----end
      mod_stan_12
    }),
    tar_target(pred_1_mod_stan_12_, {
      pred_stan <- pred_stan_
      mod_stan_12 <- mod_stan_12_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_12 <- incomplete_spatial_newdata_12_
      mod_simple_12 <- mod_simple_12_ #sampled_simple_raw_means_
      ## ---- stan_12_pred_1
      newdata <- newdata_12
      true_sum <- mod_simple_12
      stan_12_pred <- pred_stan(mod_stan_12,
        type = 1,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_stan_12"
      )
      ## ----end
      stan_12_pred
    }),
    tar_target(pred_2_mod_stan_12_, {
      pred_stan <- pred_stan_
      mod_stan_12 <- mod_stan_12_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10 <- incomplete_spatial_newdata_10_
      benthos_reefs_temporal_summary_10 <- read_sampled_reefs_data_10_plot_2_
      ## ---- stan_12_pred_2
      newdata <- newdata_10
      true_sum <- benthos_reefs_temporal_summary_10
      stan_12_pred <- pred_stan(mod_stan_12,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_stan_12"
      )
      ## ----end
      stan_12_pred
    }),
    tar_target(pred_3_mod_stan_12_, {
      pred_stan <- pred_stan_
      mod_stan_12 <- mod_stan_12_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11 <- incomplete_spatial_newdata_11_
      benthos_reefs_temporal_summary_11 <- read_sampled_reefs_data_11_plot_2_
      ## ---- stan_12_pred_3
      newdata <- newdata_11
      true_sum <- benthos_reefs_temporal_summary_11
      stan_12_pred <- pred_stan(mod_stan_12,
        type = 3,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_stan_12"
      )
      ## ----end
      stan_12_pred
    }),
    tar_target(stan_trace_12_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_stan_12 <- mod_stan_12_
      ## ---- stan_trace_12
      mod_stan_12 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_12.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_12$draws(variables = c("beta", "phi", "sd_1", "sd_2", "sd_3")) |>
        mcmc_trace() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_trace_12.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ac_12_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_stan_12 <- mod_stan_12_
      ## ---- stan_ac_12
      mod_stan_12 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_12.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_12$draws(variables = c("beta", "phi", "sd_1", "sd_2", "sd_3")) |>
        mcmc_acf() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ac_12.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_rhat_12_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_stan_12 <- mod_stan_12_
      ## ---- stan_rhat_12
      mod_stan_12 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_12.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_12 |> bayesplot::rhat() |> 
        mcmc_rhat_hist() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_rhat_12.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ess_12_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_stan_12 <- mod_stan_12_
      ## ---- stan_ess_12
      mod_stan_12 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_12.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_12 |> bayesplot::neff_ratio() |> 
        mcmc_neff_hist() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ess_12.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ppc_12_, {
      benthos_fixed_locs_obs_12 <- incomplete_spatial_data_prep_12_
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_stan_12 <- mod_stan_12_
      ## ---- stan_ppc_12
      mod_stan_12 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_12.rds")
      )
      g <- 
        bayesplot::pp_check(
          benthos_fixed_locs_obs_12$cover,
          mod_stan_12$draws("ypred", format = "matrix")[1:100, ],
          ppc_dens_overlay
        ) +
        theme_classic()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ppc_12.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(mse_1_mod_stan_12_, {
      mse_stan <- mse_stan_
      mod_stan_12 <- mod_stan_12_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_12b <- incomplete_spatial_newdata_12b_
      ## ---- stan_12_mse_1
      stan_12_mse <- mse_stan(mod_stan_12,
        newdata = newdata_12b, type = 1, model_type = ""
      )
      ## ----end
      stan_12_mse
    }),
    tar_target(mse_2_mod_stan_12_, {
      mse_stan <- mse_stan_
      mod_stan_12 <- mod_stan_12_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      ## ---- stan_12_mse_2
      stan_12_mse <- mse_stan(mod_stan_12,
        newdata = newdata_10b, type = 2, model_type = ""
      )
      ## ----end
      stan_12_mse
    }),
    tar_target(mse_3_mod_stan_12_, {
      mse_stan <- mse_stan_
      mod_stan_12 <- mod_stan_12_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      ## ---- stan_12_mse_3
      stan_12_mse <- mse_stan(mod_stan_12,
        newdata = newdata_11b, type = 3, model_type = ""
      )
      ## ----end
      stan_12_mse
    }),
    
    ## gbm ------------------------------------------------------------
    tar_target(mod_gbm_12_, {
      benthos_fixed_locs_obs_12 <- incomplete_spatial_data_prep_12_
      data_path <- incomplete_spatial_global_parameters_$data_path
      ## ---- gbm_12
      mod_gbm_12 <- gbm(cover ~ fYear,
        data =  benthos_fixed_locs_obs_12,
        distribution = "gaussian",
        n.trees = 10000,
        interaction.depth = 5,
        shrinkage = 0.001,
        bag.fraction = 0.5,
        cv.folds = 5,
        verbose = TRUE
      )
      saveRDS(mod_gbm_12,
        file = paste0(data_path, "synthetic/mod_gbm_12.rds")
      ) 
      ## ----end
      ## ---- gbm_post_12
      n.trees <- gbm.perf(mod_gbm_12, method = "cv")
      ## ----end
      list(mod_gbm_12 = mod_gbm_12, n.trees = n.trees)
    }),
    tar_target(pred_1_mod_gbm_12_, {
      pred_gbm <- pred_gbm_
      mod_gbm_12 <- mod_gbm_12_$mod_gbm_12
      n.trees <- mod_gbm_12_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_12 <- incomplete_spatial_newdata_12_
      mod_simple_12 <- mod_simple_12_ #sampled_simple_raw_means_
      ## ---- gbm_12_pred_1
      newdata <- newdata_12
      true_sum <- mod_simple_12
      gbm_12_pred <- pred_gbm(mod_gbm_12,
        n.trees = n.trees,
        type = 1,
        model_type = "",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_12"
      )
      ## ----end
      gbm_12_pred
    }),
    tar_target(pred_2_mod_gbm_12_, {
      pred_gbm <- pred_gbm_
      mod_gbm_12 <- mod_gbm_12_$mod_gbm_12
      n.trees <- mod_gbm_12_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10 <- incomplete_spatial_newdata_10_
      benthos_reefs_temporal_summary_10 <- read_sampled_reefs_data_10_plot_2_
      ## ---- gbm_12_pred_2
      newdata <- newdata_10
      true_sum <- benthos_reefs_temporal_summary_10
      gbm_12_pred <- pred_gbm(mod_gbm_12,
        n.trees = n.trees,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_12"
      )
      ## ----end
      gbm_12_pred
    }),
    tar_target(pred_3_mod_gbm_12_, {
      pred_gbm <- pred_gbm_
      mod_gbm_12 <- mod_gbm_12_$mod_gbm_12
      n.trees <- mod_gbm_12_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11 <- incomplete_spatial_newdata_11_
      benthos_reefs_temporal_summary_11 <- read_sampled_reefs_data_11_plot_2_
      ## ---- gbm_12_pred_3
      newdata <- newdata_11
      true_sum <- benthos_reefs_temporal_summary_11
      gbm_12_pred <- pred_gbm(mod_gbm_12,
        n.trees = n.trees,
        type = 3,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_12"
      )
      ## ----end
      gbm_12_pred
    }),
    tar_target(mse_1_mod_gbm_12_, {
      mse_gbm <- mse_gbm_
      mod_gbm_12 <- mod_gbm_12_$mod_gbm_12
      n.trees <- mod_gbm_12_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_12b <- incomplete_spatial_newdata_12b_
      ## ---- gbm_12_mse_1
      gbm_12_mse <- mse_gbm(mod_gbm_12,
        n.trees = n.trees,
        newdata = newdata_12b, type = 1, model_type = ""
      )
      ## ----end
      gbm_12_mse
    }),
    tar_target(mse_2_mod_gbm_12_, {
      mse_gbm <- mse_gbm_
      mod_gbm_12 <- mod_gbm_12_$mod_gbm_12
      n.trees <- mod_gbm_12_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      ## ---- gbm_12_mse_2
      gbm_12_mse <- mse_gbm(mod_gbm_12,
        n.trees = n.trees,
        newdata = newdata_10b, type = 2, model_type = ""
      )
      ## ----end
      gbm_12_mse
    }),
    tar_target(mse_3_mod_gbm_12_, {
      mse_gbm <- mse_gbm_
      mod_gbm_12 <- mod_gbm_12_$mod_gbm_12
      n.trees <- mod_gbm_12_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      ## ---- gbm_12_mse_3
      gbm_12_mse <- mse_gbm(mod_gbm_12,
        n.trees = n.trees,
        newdata = newdata_11b, type = 3, model_type = ""
      )
      ## ----end
      gbm_12_mse
    }),

    ## gbm + covartiates ----------------------------------------------
    tar_target(mod_gbm_12b_, {
      benthos_fixed_locs_obs_12 <- incomplete_spatial_data_prep_12_
      data_path <- incomplete_spatial_global_parameters_$data_path
      ## ---- gbm_12b
      mod_gbm_12b <- gbm(cover ~ fYear + Latitude + Longitude + CYC + DHW + OTHER,
        data =  benthos_fixed_locs_obs_12,
        distribution = "gaussian",
        var.monotone = c(0, 0, 0, -1, -1, -1),
        n.trees = 10000,
        interaction.depth = 5,
        shrinkage = 0.001,
        bag.fraction = 0.5,
        cv.folds = 5,
        verbose = TRUE
      )
      saveRDS(mod_gbm_12b,
        file = paste0(data_path, "synthetic/mod_gbm_12b.rds")
      )
      ## ----end
      ## ---- gbm_post_12b
      n.trees <- gbm.perf(mod_gbm_12b, method = "cv")
      ## ----end
      list(mod_gbm_12b = mod_gbm_12b, n.trees = n.trees)
    }),
    tar_target(pred_1_mod_gbm_12b_, {
      pred_gbm <- pred_gbm_
      mod_gbm_12b <- mod_gbm_12b_$mod_gbm_12
      n.trees <- mod_gbm_12b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_12b <- incomplete_spatial_newdata_12b_
      mod_simple_12 <- mod_simple_12_ #sampled_simple_raw_means_
      ## ---- gbm_12b_pred_1
      true_sum <- mod_simple_12
      newdata <- newdata_12b
      gbm_12_pred <- pred_gbm(mod_gbm_12b,
        n.trees = n.trees,
        type = 1,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_12"
      )
      ## ----end
      gbm_12_pred
    }),
    tar_target(pred_2_mod_gbm_12b_, {
      pred_gbm <- pred_gbm_
      mod_gbm_12b <- mod_gbm_12b_$mod_gbm_12
      n.trees <- mod_gbm_12b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      benthos_reefs_temporal_summary_10 <- read_sampled_reefs_data_10_plot_2_
      ## ---- gbm_12b_pred_2
      newdata <- newdata_10b
      true_sum <- benthos_reefs_temporal_summary_10
      gbm_12_pred <- pred_gbm(mod_gbm_12b,
        n.trees = n.trees,
        type = 2,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_12"
      )
      ## ----end
      gbm_12_pred
    }),
    tar_target(pred_3_mod_gbm_12b_, {
      pred_gbm <- pred_gbm_
      mod_gbm_12b <- mod_gbm_12b_$mod_gbm_12
      n.trees <- mod_gbm_12b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      benthos_reefs_temporal_summary_11 <- read_sampled_reefs_data_11_plot_2_
      ## ---- gbm_12b_pred_3
      newdata <- newdata_11b
      true_sum <- benthos_reefs_temporal_summary_11
      gbm_12_pred <- pred_gbm(mod_gbm_12b,
        n.trees = n.trees,
        type = 3,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_12"
      )
      ## ----end
      gbm_12_pred
    }),
    tar_target(infl_gbm_12b_, {
      mod_gbm_12b <- mod_gbm_12b_$mod_gbm_12b
      n.trees <- mod_gbm_12b_$n.trees
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      ## ---- gbm_infl_12b
      ## gbm_3b_infl <- gbm::relative.influence(mod_gbm_3b, n.trees = n.trees, scale = TRUE, sort = TRUE)
      infl <- summary(mod_gbm_12b, n.trees =  n.trees, plot = FALSE)
      g <-
        infl |>
        as.data.frame() |>
        arrange(rel.inf) |>
        mutate(var = factor(var, levels = unique(var))) |>
        ggplot(aes(y = var, x = rel.inf)) +
        geom_bar(stat = "identity") +
        scale_y_discrete("") +
        scale_x_continuous("Relative influence") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_infl_mod_gbm_12b.png"
        ),
        g,
        width = 6, height = 4, dpi = 72
      )
      ## ----end
    }),
    tar_target(mse_1_mod_gbm_12b_, {
      mse_gbm <- mse_gbm_
      mod_gbm_12 <- mod_gbm_12b_$mod_gbm_12b
      n.trees <- mod_gbm_12b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_12b <- incomplete_spatial_newdata_12b_
      ## ---- gbm_12b_mse_1
      gbm_12_mse <- mse_gbm(mod_gbm_12,
        n.trees = n.trees,
        newdata = newdata_12b, type = 1, model_type = "covariates"
      )
      ## ----end
      gbm_12_mse
    }),
    tar_target(mse_2_mod_gbm_12b_, {
      mse_gbm <- mse_gbm_
      mod_gbm_12 <- mod_gbm_12b_$mod_gbm_12b
      n.trees <- mod_gbm_12b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      ## ---- gbm_12b_mse_2
      gbm_12_mse <- mse_gbm(mod_gbm_12,
        n.trees = n.trees,
        newdata = newdata_10b, type = 2, model_type = "covariates"
      )
      ## ----end
      gbm_12_mse
    }),
    tar_target(mse_3_mod_gbm_12b_, {
      mse_gbm <- mse_gbm_
      mod_gbm_12 <- mod_gbm_12b_$mod_gbm_12b
      n.trees <- mod_gbm_12b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      ## ---- gbm_12b_mse_3
      gbm_12_mse <- mse_gbm(mod_gbm_12,
        n.trees = n.trees,
        newdata = newdata_11b, type = 3, model_type = "covariates"
      )
      ## ----end
      gbm_12_mse
    }),

    ## dbarts --------------------------------------------------------
    tar_target(mod_dbarts_12_, {
      benthos_fixed_locs_obs_12 <- incomplete_spatial_data_prep_12_
      benthos_fixed_locs_obs_10 <- incomplete_spatial_data_prep_10_
      benthos_fixed_locs_obs_11 <- incomplete_spatial_data_prep_11_
      data_path <- incomplete_spatial_global_parameters_$data_path
      newdata_12 <- incomplete_spatial_newdata_12_
      newdata_12b <- incomplete_spatial_newdata_12b_
      newdata_12_a <- incomplete_spatial_data_prep_12_#incomplete_spatial_newdata_12b_
      newdata_10 <- incomplete_spatial_newdata_10_
      newdata_10_a <- incomplete_spatial_newdata_10b_
      newdata_11 <- incomplete_spatial_newdata_11_
      newdata_11_a <- incomplete_spatial_newdata_11b_
      ## ---- dbarts_12
      print(head(benthos_fixed_locs_obs_12))
      mod_dbarts_12 <- bart2(log(cover) ~ fYear,
        data =   benthos_fixed_locs_obs_12,
        keepTrees = TRUE
      )
      saveRDS(mod_dbarts_12,
        file = paste0(data_path, "synthetic/mod_dbarts_12.rds")
      ) 
      ## ----end
      ## Unfortunately, the next part must be in the same tar_target
      ## due to the way dbarts stores pointers - they cannot be stored
      ## ---- dbarts_pred_12
      preds <- predict(mod_dbarts_12, newdata_12, type = "ev") |>
        exp()
      saveRDS(preds,
        file = paste0(data_path, "synthetic/mod_dbarts_12_preds.rds")
      ) 
      preds_sum <- preds |>
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_12_preds_sum.rds")
      ) 

      newdata <- newdata_12_a |>
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
      preds_a <- predict(mod_dbarts_12, newdata, type = "ev") |>
        exp() 
      saveRDS(preds_a,
        file = paste0(data_path, "synthetic/mod_dbarts_12_a_preds.rds")
      ) 

      newdata_12b <- newdata_12b |>
        mutate(fYear = factor(Year))
      preds_12b <- predict(mod_dbarts_12, newdata_12b, type = "ev") |>
        exp() 
      saveRDS(preds_12b,
        file = paste0(data_path, "synthetic/mod_dbarts_12b_preds.rds")
      ) 


      
      preds_10 <- predict(mod_dbarts_12, newdata_10, type = "ev") |>
        exp()
      saveRDS(preds_10,
        file = paste0(data_path, "synthetic/mod_dbarts_10_preds.rds")
      ) 
      preds_10_sum <- preds_10 |>
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_10_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_10_preds_sum.rds")
      ) 
      newdata_10_a <- newdata_10_a |>
        ungroup() |>
        mutate(fYear = factor(Year))
      preds_10_a <- predict(mod_dbarts_12, newdata_10_a, type = "ev") |>
        exp()
      saveRDS(preds_10_a,
        file = paste0(data_path, "synthetic/mod_dbarts_10_a_preds.rds")
      ) 

      preds_11 <- predict(mod_dbarts_12, newdata_11, type = "ev") |>
        exp()
      saveRDS(preds_11,
        file = paste0(data_path, "synthetic/mod_dbarts_11_preds.rds")
      ) 
      preds_11_sum <- preds_11 |> 
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_11_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_11_preds_sum.rds")
      ) 
      newdata_11_a <- newdata_11_a |>
        ungroup() |>
        mutate(fYear = factor(Year))
      preds_11_a <- predict(mod_dbarts_12, newdata_11_a, type = "ev") |>
        exp()
      saveRDS(preds_11_a,
        file = paste0(data_path, "synthetic/mod_dbarts_11_a_preds.rds")
      ) 
      ## ----end
      list(
        mod_dbarts_12 = mod_dbarts_12,
        preds = preds,
        preds_sum = preds_sum,
        preds_a = preds_a,
        preds_10 = preds_10,
        preds_10_sum = preds_10_sum,
        preds_10_a = preds_10_a,
        preds_11 = preds_11,
        preds_11_sum = preds_11_sum,
        preds_11_a = preds_11_a
      )
    }),
    tar_target(pred_1_mod_dbarts_12_, {
      pred_dbarts <- pred_dbarts_
      preds <- mod_dbarts_12_$preds_sum
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_12 <- incomplete_spatial_newdata_12_
      mod_simple_12 <- mod_simple_12_ #sampled_simple_raw_means_
      ## ---- dbarts_12_pred_1
      newdata <- newdata_12
      true_sum <- mod_simple_12
      dbarts_12_pred <- pred_dbarts(
        preds,
        type = 1,
        model_type = "",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_12"
      )
      ## ----end
      dbarts_12_pred
    }),
    tar_target(pred_2_mod_dbarts_12_, {
      pred_dbarts <- pred_dbarts_
      preds <- mod_dbarts_12_$preds_10_sum
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10 <- incomplete_spatial_newdata_10_
      benthos_reefs_temporal_summary_10 <- read_sampled_reefs_data_10_plot_2_
      ## ---- dbarts_12_pred_2
      newdata <- newdata_10
      true_sum <- benthos_reefs_temporal_summary_10
      dbarts_12_pred <- pred_dbarts(
        preds,
        type = 2,
        model_type = "",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_12"
      )
      ## ----end
      dbarts_12_pred
    }),
    tar_target(pred_3_mod_dbarts_12_, {
      pred_dbarts <- pred_dbarts_
      preds <- mod_dbarts_12_$preds_11_sum
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11 <- incomplete_spatial_newdata_11_
      benthos_reefs_temporal_summary_11 <- read_sampled_reefs_data_11_plot_2_
      ## ---- dbarts_12_pred_3
      newdata <- newdata_11
      true_sum <- benthos_reefs_temporal_summary_11
      dbarts_12_pred <- pred_dbarts(
        preds,
        type = 3,
        model_type = "",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_12"
      )
      ## ----end
      dbarts_12_pred
    }),
    tar_target(mse_1_mod_dbarts_12_, {
      mse_dbarts <- mse_dbarts_
      preds <- mod_dbarts_12_$preds_a
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_12b <- incomplete_spatial_newdata_12b_
      ## ---- dbarts_12_mse_1
      newdata_12b <- newdata_12b |>
        mutate(fYear = factor(Year)) |> 
        ungroup()
      dbarts_12_mse <- mse_dbarts(preds,
        newdata = newdata_12b, type = 1, model_type = ""
      )
      ## ----end
      dbarts_12_mse
    }),
    tar_target(mse_2_mod_dbarts_12_, {
      mse_dbarts <- mse_dbarts_
      preds <- mod_dbarts_12_$preds_10_a
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_ #incomplete_spatial_newdata_12b_
      ## ---- dbarts_12_mse_2
      dbarts_12_mse <- mse_dbarts(preds,
        newdata = newdata_10b, type = 2, model_type = ""
      )
      ## ----end
      dbarts_12_mse
    }),
    tar_target(mse_3_mod_dbarts_12_, {
      mse_dbarts <- mse_dbarts_
      preds <- mod_dbarts_12_$preds_11_a
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_ #incomplete_spatial_newdata_12b_
      ## ---- dbarts_12_mse_3
      ## newdata <- newdata |>
      ##   ungroup() |>
      ##   mutate(fYear = factor(Year))
      dbarts_12_mse <- mse_dbarts(preds,
        newdata = newdata_11b, type = 3, model_type = ""
      )
      ## ----end
      dbarts_12_mse
    }),

    ## dbarts + covariates --------------------------------------------
    tar_target(mod_dbarts_12b_, {
      benthos_fixed_locs_obs_12 <- incomplete_spatial_data_prep_12_
      benthos_fixed_locs_obs_10 <- incomplete_spatial_data_prep_10_
      benthos_fixed_locs_obs_11 <- incomplete_spatial_data_prep_11_
      data_path <- incomplete_spatial_global_parameters_$data_path
      newdata_12b <- incomplete_spatial_newdata_12b_
      ## newdata_10 <- incomplete_spatial_newdata_10_
      newdata_10b <- incomplete_spatial_newdata_10b_
      newdata_11b <- incomplete_spatial_newdata_11b_
      ## ---- dbarts_12b
      print(head(benthos_fixed_locs_obs_12))
      mod_dbarts_12b <- bart2(log(cover) ~ fYear + Latitude + Longitude + CYC + DHW + OTHER,
        data =   benthos_fixed_locs_obs_12,
        keepTrees = TRUE
      )
      saveRDS(mod_dbarts_12b,
        file = paste0(data_path, "synthetic/mod_dbarts_12b.rds")
      ) 
      ## ----end
      ## Unfortunately, the next part must be in the same tar_target
      ## due to the way dbarts stores pointers - they cannot be stored
      ## ---- dbarts_pred_12b
      newdata_12b <- newdata_12b |>
        mutate(fYear = factor(Year))
      preds_12b <- predict(mod_dbarts_12b, newdata_12b, type = "ev") |>
        exp() 
      saveRDS(preds_12b,
        file = paste0(data_path, "synthetic/mod_dbarts_12b_preds_12b.rds")
      ) 
      preds_12b_sum <- preds_12b |> summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_12b_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_12b_preds_12b_sum.rds")
      )
      
      ## newdata_10 <- newdata_10 |>
      ##   mutate(fYear = factor(Year))
      ## preds_10 <- predict(mod_dbarts_12b, newdata_10, type = "ev") |>
      ##   exp() 
      ## saveRDS(preds_10,
      ##   file = paste0(data_path, "synthetic/mod_dbarts_10b_preds.rds")
      ## ) 
      ## preds_10_sum <- preds_10 |> 
      ##   summarise_draws(median, HDInterval::hdi)
      ## saveRDS(preds_10_sum,
      ##   file = paste0(data_path, "synthetic/mod_dbarts_10b_preds_sum.rds")
      ## ) 

      newdata_10b <- newdata_10b |>
        mutate(fYear = factor(Year))
      preds_10b <- predict(mod_dbarts_12b, newdata_10b, type = "ev") |>
        exp() 
      saveRDS(preds_10b,
        file = paste0(data_path, "synthetic/mod_dbarts_12b_preds_10b.rds")
      ) 
      preds_10b_sum <- preds_10b |> 
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_10b_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_12b_preds_10b_sum.rds")
      ) 

      newdata_11b <- newdata_11b |>
        mutate(fYear = factor(Year))
      preds_11b <- predict(mod_dbarts_12b, newdata_11b, type = "ev") |>
        exp()
      saveRDS(preds_11b,
        file = paste0(data_path, "synthetic/mod_dbarts_12b_preds_11b.rds")
      ) 
      preds_11b_sum <- preds_11b |>
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_11b_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_12b_preds_11b_sum.rds")
      ) 
      ## ----end
      list(
        mod_dbarts_12 = mod_dbarts_12b,
        preds_12b = preds_12b,
        preds_12b_sum = preds_12b_sum,
        preds_10b = preds_10b,
        preds_10b_sum = preds_10b_sum,
        preds_11b = preds_11b,
        preds_11b_sum = preds_11b_sum
      )
    }),
    tar_target(pred_1_mod_dbarts_12b_, {
      pred_dbarts <- pred_dbarts_
      preds_12b <- mod_dbarts_12b_$preds_12b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_12b <- incomplete_spatial_newdata_12b_
      mod_simple_12 <- mod_simple_12_ #sampled_simple_raw_means_
      ## ---- dbarts_12b_pred_1
      true_sum <- mod_simple_12
      newdata <- newdata_12b
      dbarts_12_pred <- pred_dbarts(
        preds_12b,
        type = 1,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_12"
      )
      ## ----end
      dbarts_12_pred
    }),
    tar_target(pred_2_mod_dbarts_12b_, {
      pred_dbarts <- pred_dbarts_
      preds_10b <- mod_dbarts_12b_$preds_10b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      benthos_reefs_temporal_summary_10 <- read_sampled_reefs_data_10_plot_2_
      ## ---- dbarts_12b_pred_2
      newdata <- newdata_10b
      true_sum <- benthos_reefs_temporal_summary_10
      dbarts_12_pred <- pred_dbarts(
        preds_10b, 
        type = 2, 
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_12"
      )
      ## ----end
      dbarts_12_pred
    }),
    tar_target(pred_3_mod_dbarts_12b_, {
      pred_dbarts <- pred_dbarts_
      preds_11b <- mod_dbarts_12b_$preds_11b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      benthos_reefs_temporal_summary_11 <- read_sampled_reefs_data_11_plot_2_
      ## ---- dbarts_12b_pred_3
      newdata <- newdata_11b
      true_sum <- benthos_reefs_temporal_summary_11
      dbarts_12_pred <- pred_dbarts(
        preds_11b,
        type = 3,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_12"
      )
      ## ----end
      dbarts_12_pred
    }),
    tar_target(mse_1_mod_dbarts_12b_, {
      mse_dbarts <- mse_dbarts_
      preds_12b <- mod_dbarts_12b_$preds_12b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata <- incomplete_spatial_data_prep_12_ #incomplete_spatial_newdata_12b_ #incomplete_spatial_newdata_12b_
      newdata_12b <- incomplete_spatial_newdata_12b_
      ## ---- dbarts_12b_mse_1
      dbarts_12_mse <- mse_dbarts(preds_12b,
        newdata = newdata_12b, type = 1, model_type = "covariates"
      )
      ## ----end
      dbarts_12_mse
    }),
    tar_target(mse_2_mod_dbarts_12b_, {
      mse_dbarts <- mse_dbarts_
      preds_10b <- mod_dbarts_12b_$preds_10b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      ## ---- dbarts_12b_mse_2
      dbarts_12_mse <- mse_dbarts(preds_10b,
        newdata = newdata_10b, type = 2, model_type = "covariates"
      )
      ## ----end
      dbarts_12_mse
    }),
    tar_target(mse_3_mod_dbarts_12b_, {
      mse_dbarts <- mse_dbarts_
      preds_11b <- mod_dbarts_12b_$preds_11b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      ## ---- dbarts_12b_mse_3
      dbarts_12_mse <- mse_dbarts(preds_11b,
        newdata = newdata_11b, type = 3, model_type = "covariates"
      )
      ## ----end
      dbarts_12_mse
    }),
    
    ## xgboost + covariates -------------------------------------------
    tar_target(mod_xgboost_12b_prep_, {
      benthos_fixed_locs_obs_12 <- incomplete_spatial_data_prep_12_
      data_path <- incomplete_spatial_global_parameters_$data_path
      ## ---- xgboost_12_prep
      data_train_12b <- benthos_fixed_locs_obs_12 |>
        dplyr::select(cover, Year, Latitude, Longitude, CYC, DHW, OTHER) 
      ## ----end
      data_train_12b
    }),
    tar_target(mod_xgboost_12b_tune_, {
      data_train_12b <- mod_xgboost_12b_prep_
      data_path <- incomplete_spatial_global_parameters_$data_path
      ## ---- xgboost_12_tune
      ## Define the recipe
      tune_recipe <- recipe(cover ~ ., data = data_train_12b) |> 
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
        resamples = vfold_cv(data_train_12b, v = 5),
        grid = tune_grid_values)
      ## Get best set of parameters
      model_hyperparams <-
        select_best(tuned_results, metric = "rmse") |> 
        select(-".config") |> 
        as_tibble() |>
        mutate(nb_training = nrow(data_train_12b),
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
    tar_target(mod_xgboost_12b_fit_, {
      data_train_12b <- mod_xgboost_12b_prep_
      data_path <- incomplete_spatial_global_parameters_$data_path
      tune_recipe = mod_xgboost_12b_tune_$tune_recipe
      tune_model = mod_xgboost_12b_tune_$tune_model
      tune_workflow = mod_xgboost_12b_tune_$tune_workflow
      tune_grid_values = mod_xgboost_12b_tune_$tune_grid_values
      tuned_results = mod_xgboost_12b_tune_$tuned_results
      model_hyperparams = mod_xgboost_12b_tune_$model_hyperparams
      ## ---- xgboost_12_fit
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
      final_fitted_12b <- tune_workflow |>
        fit(data_train_12b)
      ## ----end
      final_fitted_12b
    }),
    tar_target(pred_1_mod_xgboost_12b_, {
      pred_xgboost <- pred_xgboost_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_12b <- incomplete_spatial_newdata_12b_
      final_fitted_12b <- mod_xgboost_12b_fit_
      mod_simple_12 <- mod_simple_12_ #sampled_simple_raw_means_
      ## ---- xgboost_12b_pred_1
      true_sum <- mod_simple_12
      newdata <- newdata_12b
      xgboost_12b_pred <- pred_xgboost(
        final_fitted_12b,
        type = 1,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_xgboost_12"
      )
      ## ----end
      xgboost_12b_pred
    }),
    tar_target(mse_1_mod_xgboost_12b_, {
      mse_xgboost <- mse_xgboost_
      final_fitted_12b <- mod_xgboost_12b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_12b <- incomplete_spatial_newdata_12b_
      ## ---- xgboost_12b_mse_1
      xgboost_12_mse <- mse_xgboost(final_fitted_12b,
        newdata = newdata_12b, type = 1, model_type = "covariates"
      )
      ## ----end
      xgboost_12_mse
    }),
    tar_target(pred_2_mod_xgboost_12b_, {
      pred_xgboost <- pred_xgboost_
      final_fitted_12b <- mod_xgboost_12b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      benthos_reefs_temporal_summary_10 <- read_sampled_reefs_data_10_plot_2_
      ## ---- xgboost_12b_pred_2
      newdata <- newdata_10b
      true_sum <- benthos_reefs_temporal_summary_10
      xgboost_12_pred <- pred_xgboost(final_fitted_12b,
        type = 2,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_xgboost_12"
      )
      ## ----end
      xgboost_12_pred
    }),
    tar_target(mse_2_mod_xgboost_12b_, {
      mse_xgboost <- mse_xgboost_
      final_fitted_12b <- mod_xgboost_12b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      ## ---- xgboost_12b_mse_2
      xgboost_12_mse <- mse_xgboost(final_fitted_12b,
        newdata = newdata_10b, type = 2, model_type = "covariates"
      )
      ## ----end
      xgboost_12_mse
    }),
    tar_target(pred_3_mod_xgboost_12b_, {
      pred_xgboost <- pred_xgboost_
      final_fitted_12b <- mod_xgboost_12b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      benthos_reefs_temporal_summary_11 <- read_sampled_reefs_data_11_plot_2_
      ## ---- xgboost_12b_pred_3
      newdata <- newdata_11b
      true_sum <- benthos_reefs_temporal_summary_11
      xgboost_12_pred <- pred_xgboost(final_fitted_12b,
        type = 3,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_xgboost_12"
      )
      ## ----end
      xgboost_12_pred
    }),
    tar_target(mse_3_mod_xgboost_12b_, {
      mse_xgboost <- mse_xgboost_
      final_fitted_12b <- mod_xgboost_12b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      ## ---- xgboost_12b_mse_3
      xgboost_12_mse <- mse_xgboost(final_fitted_12b,
        newdata = newdata_11b, type = 3, model_type = "covariates"
      )
      ## ----end
      xgboost_12_mse
    }),

    ## Comparisons ----------------------------------------------------
    tar_target(mse_1_mod_pymc_barts_12b_file_,
      paste0(
        incomplete_spatial_global_parameters_$data_path,
        "modelled/pymc_bart_12b_mse_1.csv"
      ), format = "file"),
    tar_target(mse_1_mod_pymc_barts_12b_, {
      read_csv(file = mse_1_mod_pymc_barts_12b_file_) 
    }),
    tar_target(mse_2_mod_pymc_barts_12b_file_,
      paste0(
        incomplete_spatial_global_parameters_$data_path,
        "modelled/pymc_bart_12b_mse_2.csv"
      ), format = "file"),
    tar_target(mse_2_mod_pymc_barts_12b_, {
      read_csv(file = mse_2_mod_pymc_barts_12b_file_) 
    }),
    tar_target(mse_3_mod_pymc_barts_12b_file_,
      paste0(
        incomplete_spatial_global_parameters_$data_path,
        "modelled/pymc_bart_12b_mse_3.csv"
      ), format = "file"),
    tar_target(mse_3_mod_pymc_barts_12b_, {
      read_csv(file = mse_3_mod_pymc_barts_12b_file_) 
    }),
    tar_target(comparisons_12_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      mse_1_mod_glmmTMB_12 <- mse_1_mod_glmmTMB_12_
      mse_2_mod_glmmTMB_12 <- mse_2_mod_glmmTMB_12_
      mse_3_mod_glmmTMB_12 <- mse_3_mod_glmmTMB_12_
      mse_1_mod_glmmTMB_12b <- mse_1_mod_glmmTMB_12b_
      mse_2_mod_glmmTMB_12b <- mse_2_mod_glmmTMB_12b_
      mse_3_mod_glmmTMB_12b <- mse_3_mod_glmmTMB_12b_

      mse_1_mod_brms_12 <- mse_1_mod_brms_12_
      mse_2_mod_brms_12 <- mse_2_mod_brms_12_
      mse_3_mod_brms_12 <- mse_3_mod_brms_12_
      mse_1_mod_brms_12b <- mse_1_mod_brms_12b_
      mse_2_mod_brms_12b <- mse_2_mod_brms_12b_
      mse_3_mod_brms_12b <- mse_3_mod_brms_12b_

      mse_1_mod_stan_12 <- mse_1_mod_stan_12_
      mse_2_mod_stan_12 <- mse_2_mod_stan_12_
      mse_3_mod_stan_12 <- mse_3_mod_stan_12_

      mse_1_mod_gbm_12 <- mse_1_mod_gbm_12_
      mse_2_mod_gbm_12 <- mse_2_mod_gbm_12_
      mse_3_mod_gbm_12 <- mse_3_mod_gbm_12_
      mse_1_mod_gbm_12b <- mse_1_mod_gbm_12b_
      mse_2_mod_gbm_12b <- mse_2_mod_gbm_12b_
      mse_3_mod_gbm_12b <- mse_3_mod_gbm_12b_

      mse_1_mod_dbarts_12 <- mse_1_mod_dbarts_12_
      mse_2_mod_dbarts_12 <- mse_2_mod_dbarts_12_
      mse_3_mod_dbarts_12 <- mse_3_mod_dbarts_12_
      mse_1_mod_dbarts_12b <- mse_1_mod_dbarts_12b_
      mse_2_mod_dbarts_12b <- mse_2_mod_dbarts_12b_
      mse_3_mod_dbarts_12b <- mse_3_mod_dbarts_12b_

      mse_1_mod_xgboost_12b <- mse_1_mod_xgboost_12b_
      mse_2_mod_xgboost_12b <- mse_2_mod_xgboost_12b_
      mse_3_mod_xgboost_12b <- mse_3_mod_xgboost_12b_

      mse_1_mod_pymc_barts_12b <- mse_1_mod_pymc_barts_12b_
      mse_2_mod_pymc_barts_12b <- mse_2_mod_pymc_barts_12b_
      mse_3_mod_pymc_barts_12b <- mse_3_mod_pymc_barts_12b_
      ## ---- comparisons_12
      mse_1_mod_pymc_barts_12b <- read_csv(
        file = paste0(data_path, "modelled/pymc_bart_12b_mse_1.csv")
      )   
      mse_2_mod_pymc_barts_12b <- read_csv(
        file = paste0(data_path, "modelled/pymc_bart_12b_mse_2.csv")
      )   
      mse_3_mod_pymc_barts_12b <- read_csv(
        file = paste0(data_path, "modelled/pymc_bart_12b_mse_3.csv")
      )   
      comparisons_12 <- bind_rows(
        mse_1_mod_glmmTMB_12,
        mse_2_mod_glmmTMB_12,
        mse_3_mod_glmmTMB_12,
        mse_1_mod_glmmTMB_12b,
        mse_2_mod_glmmTMB_12b,
        mse_3_mod_glmmTMB_12b,

        mse_1_mod_brms_12,
        mse_2_mod_brms_12,
        mse_3_mod_brms_12,
        mse_1_mod_brms_12b,
        mse_2_mod_brms_12b,
        mse_3_mod_brms_12b,
        
        mse_1_mod_stan_12,
        mse_2_mod_stan_12,
        mse_3_mod_stan_12,
        
        mse_1_mod_gbm_12,
        mse_2_mod_gbm_12,
        mse_3_mod_gbm_12,
        mse_1_mod_gbm_12b,
        mse_2_mod_gbm_12b,
        mse_3_mod_gbm_12b,

        mse_1_mod_dbarts_12,
        mse_2_mod_dbarts_12,
        mse_3_mod_dbarts_12,
        mse_1_mod_dbarts_12b,
        mse_2_mod_dbarts_12b,
        mse_3_mod_dbarts_12b,

        mse_1_mod_xgboost_12b,
        mse_2_mod_xgboost_12b,
        mse_3_mod_xgboost_12b,

        mse_1_mod_pymc_barts_12b,
        mse_2_mod_pymc_barts_12b,
        mse_3_mod_pymc_barts_12b
        )
      saveRDS(comparisons_12,
        file = paste0(data_path, "synthetic/comparisons_12.rds")
      ) 
       
      comps_12 <- 
        comparisons_12 |>
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
          type == 2 ~ "Predicting Southern reefs",
          type == 3 ~ "Predicting Northern reefs"
        )) 
      saveRDS(comps_12,
        file = paste0(data_path, "synthetic/comps_12.rds")
      )  
      ## ----end
      comps_12
    }),
    tar_target(comparisons_12_plots_, {
      comps_12 <- comparisons_12_
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- comparisons_12_tab
      g <-
        comps_12 |> 
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
          fig_path, "mse_12_1.png"
        ),
        g,
        width = 10, height = 6, dpi = 100
      )
      g <- 
        comps_12 |> 
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
          fig_path, "mse_12_2.png"
        ),
        g,
        width = 10, height = 6, dpi = 100
      )
      ## ----end
      comps_12
    }),
    
    ## Southern subdomain (Sampled reefs) 
    ## simple ---------------------------------------------------------
    tar_target(mod_simple_13_, {
      benthos_fixed_locs_obs_13 <- incomplete_spatial_data_prep_13_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- simple_12
      mod_simple_13 <- benthos_fixed_locs_obs_13 |>
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
      saveRDS(mod_simple_13,
        file = paste0(data_path, "synthetic/mod_simple_13.rds")
      ) 
      ## ----end
      mod_simple_13
    }),

    ## glmmTMB --------------------------------------------------------
    tar_target(mod_glmmTMB_13_, {
      benthos_fixed_locs_obs_13 <- incomplete_spatial_data_prep_13_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_13
      mod_glmmTMB_13 <- glmmTMB(cover ~ fYear + (1 | Site) + (1 | Transect),
        data = benthos_fixed_locs_obs_13,
        family = "beta_family"
      )
      saveRDS(mod_glmmTMB_13,
        file = paste0(data_path, "synthetic/mod_glmmTMB_13.rds")
      ) 
      ## ----end
      mod_glmmTMB_13
    }),
    tar_target(mod_glmmTMB_13_sample_data_, {
      benthos_fixed_locs_obs_13 <- incomplete_spatial_data_prep_13_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_13 sample data
      benthos_fixed_locs_obs_13 
      benthos_fixed_locs_obs_13 |>
        group_by(Year) |>
        summarise(
          Mean = mean(cover),
          Median = median(cover)
          )
      ## ----end
      benthos_fixed_locs_obs_13
    }),
    tar_target(mod_glmmTMB_13_sample_data_summary_10_, {
      benthos_reefs_temporal_summary_10 <- read_sampled_reefs_data_10_plot_2_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_13 sample data summary 10
      benthos_reefs_temporal_summary_10 <- readRDS(
        file = paste0(data_path, "synthetic/benthos_reefs_temporal_summary_10.rds")
      )
      benthos_reefs_temporal_summary_10
      ## ----end
      benthos_reefs_temporal_summary_10
    }),
    tar_target(mod_glmmTMB_13_sample_data_summary_11_, {
      benthos_reefs_temporal_summary_11 <- read_sampled_reefs_data_11_plot_2_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_13 sample data summary 11
      benthos_reefs_temporal_summary_11 <- readRDS(
        file = paste0(data_path, "synthetic/benthos_reefs_temporal_summary_11.rds")
      )
      benthos_reefs_temporal_summary_11
      ## ----end
      benthos_reefs_temporal_summary_11
    }),
    tar_target(mod_glmmTMB_13_newdata_12_, {
      newdata_12 <- incomplete_spatial_newdata_12_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_13 newdata 12
      newdata_12 
      ## ----end
      newdata_12
    }),
    tar_target(mod_glmmTMB_13_newdata_12b_, {
      newdata_12b <- incomplete_spatial_newdata_12b_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_13 newdata 12b
      newdata_12b 
      ## ----end
      newdata_12b
    }),
    tar_target(mod_glmmTMB_13_newdata_13_, {
      newdata_13 <- incomplete_spatial_newdata_13_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_13 newdata 13
      newdata_13 
      ## ----end
      newdata_13
    }),
    tar_target(mod_glmmTMB_13_newdata_13b_, {
      newdata_13b <- incomplete_spatial_newdata_13_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_13 newdata 13b
      newdata_13b 
      ## ----end
      newdata_13b
    }),
    tar_target(mod_glmmTMB_13_newdata_10_, {
      newdata_10 <- incomplete_spatial_newdata_10_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_13 newdata 10
      newdata_10 
      ## ----end
      newdata_10
    }),
    tar_target(mod_glmmTMB_13_newdata_10b_, {
      newdata_10b <- incomplete_spatial_newdata_10b_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_13 newdata 10b
      newdata_10b 
      ## ----end
      newdata_10b
    }),
    tar_target(mod_glmmTMB_13_newdata_11_, {
      newdata_11 <- incomplete_spatial_newdata_11_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_13 newdata 11
      newdata_11 
      ## ----end
      newdata_11
    }),
    tar_target(mod_glmmTMB_13_newdata_11b_, {
      newdata_11b <- incomplete_spatial_newdata_11b_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_13 newdata 11b
      newdata_11b 
      ## ----end
      newdata_11b
    }),

    tar_target(dharma_mod_glmmTMB_13_, {
      DHARMa_glmmTMB <- DHARMa_glmmTMB_
      mod_glmmTMB_13 <- mod_glmmTMB_13_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- glmmTMB_13_dharma
      glmmTMB_13_dharma <- DHARMa_glmmTMB(mod_glmmTMB_13,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_13")
      ## ----end
      glmmTMB_13_dharma
    }),
    tar_target(pred_1_mod_glmmTMB_13_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_13 <- mod_glmmTMB_13_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_13 <- incomplete_spatial_newdata_13_
      mod_simple_13 <- mod_simple_13_ #sampled_simple_raw_means_
      ## ---- glmmTMB_13_pred_1
      newdata <- newdata_13
      true_sum <- mod_simple_13
      glmmTMB_13_pred <- pred_glmmTMB(mod_glmmTMB_13,
        type = 1,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_13"
      )
      ## ----end
      glmmTMB_13_pred
    }),
    tar_target(pred_2_mod_glmmTMB_13_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_13 <- mod_glmmTMB_13_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11 <- incomplete_spatial_newdata_11_
      benthos_reefs_temporal_summary_11 <- read_sampled_reefs_data_11_plot_2_
      ## ---- glmmTMB_13_pred_2
      newdata <- newdata_11
      true_sum <- benthos_reefs_temporal_summary_11
      glmmTMB_13_pred <- pred_glmmTMB(mod_glmmTMB_13,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_13"
      )
      ## ----end
      glmmTMB_13_pred
    }),
    tar_target(pred_3_mod_glmmTMB_13_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_13 <- mod_glmmTMB_13_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10 <- incomplete_spatial_newdata_10_
      benthos_reefs_temporal_summary_10 <- read_sampled_reefs_data_10_plot_2_
      ## ---- glmmTMB_13_pred_3
      newdata <- newdata_10
      true_sum <- benthos_reefs_temporal_summary_10
      glmmTMB_13_pred <- pred_glmmTMB(mod_glmmTMB_13,
        type = 3,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_13"
      )
      ## ----end
      glmmTMB_13_pred
    }),
    tar_target(mse_1_mod_glmmTMB_13_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_13 <- mod_glmmTMB_13_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_13b <- incomplete_spatial_newdata_13b_
      ## ---- glmmTMB_13_mse_1
      glmmTMB_13_mse <- mse_glmmTMB(mod_glmmTMB_13,
        newdata = newdata_13b, type = 1, model_type = ""
      )
      ## ----end
      glmmTMB_13_mse
    }),
    tar_target(mse_2_mod_glmmTMB_13_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_13 <- mod_glmmTMB_13_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      ## ---- glmmTMB_13_mse_2
      glmmTMB_13_mse <- mse_glmmTMB(mod_glmmTMB_13,
        newdata = newdata_11b, type = 2, model_type = ""
      )
      ## ----end
      glmmTMB_13_mse
    }),
    tar_target(mse_3_mod_glmmTMB_13_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_13 <- mod_glmmTMB_13_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      ## ---- glmmTMB_13_mse_3
      glmmTMB_13_mse <- mse_glmmTMB(mod_glmmTMB_13,
        newdata = newdata_10b, type = 3, model_type = ""
      )
      ## ----end
      glmmTMB_13_mse
    }),

    ## glmmTMB + covariates -------------------------------------------
    tar_target(mod_glmmTMB_13b_, {
      benthos_fixed_locs_obs_13 <- incomplete_spatial_data_prep_13_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_13b
      mod_glmmTMB_13b <- glmmTMB(cover ~ fYear + CYC + DHW + OTHER + Latitude + Longitude +
                                   (1 | Site) + (1 | Transect),
        data = benthos_fixed_locs_obs_13,
        family = "beta_family"
      )
      saveRDS(mod_glmmTMB_13b,
        file = paste0(data_path, "synthetic/mod_glmmTMB_13b.rds")
      ) 
      ## ----end
      mod_glmmTMB_13b
    }),
    tar_target(dharma_mod_glmmTMB_13b_, {
      DHARMa_glmmTMB <- DHARMa_glmmTMB_
      mod_glmmTMB_13b <- mod_glmmTMB_13b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- glmmTMB_13b_dharma
      glmmTMB_13b_dharma <- DHARMa_glmmTMB(mod_glmmTMB_13b,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_13b")
      ## ----end
      glmmTMB_13b_dharma
    }),
    tar_target(pred_1_mod_glmmTMB_13b_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_13b <- mod_glmmTMB_13b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_13b <- incomplete_spatial_newdata_13b_
      mod_simple_13 <- mod_simple_13_ #sampled_simple_raw_means_
      ## ---- glmmTMB_13b_pred_1
      true_sum <- mod_simple_13
      newdata <- newdata_13b
      glmmTMB_13b_pred <- pred_glmmTMB(mod_glmmTMB_13b,
        type = 1,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_13"
      )
      ## ----end
      glmmTMB_13b_pred
    }),
    tar_target(pred_2_mod_glmmTMB_13b_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_13 <- mod_glmmTMB_13b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      benthos_reefs_temporal_summary_11 <- read_sampled_reefs_data_11_plot_2_
      ## ---- glmmTMB_13b_pred_2
      newdata <- newdata_11b
      true_sum <- benthos_reefs_temporal_summary_11
      glmmTMB_13b_pred <- pred_glmmTMB(mod_glmmTMB_13,
        type = 2,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_13"
      ) 
      ## ----end
      glmmTMB_13b_pred
    }),
    tar_target(pred_3_mod_glmmTMB_13b_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_13 <- mod_glmmTMB_13b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      benthos_reefs_temporal_summary_10 <- read_sampled_reefs_data_10_plot_2_
      ## ---- glmmTMB_13b_pred_3
      newdata <- newdata_10b
      true_sum <- benthos_reefs_temporal_summary_10
      glmmTMB_13b_pred <- pred_glmmTMB(mod_glmmTMB_13,
        type = 3,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_13"
      )
      ## ----end
      glmmTMB_13b_pred
    }),
    tar_target(mse_1_mod_glmmTMB_13b_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_13 <- mod_glmmTMB_13b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_13b <- incomplete_spatial_newdata_13b_
      ## ---- glmmTMB_13b_mse_1
      glmmTMB_13_mse <- mse_glmmTMB(mod_glmmTMB_13,
        newdata = newdata_13b, type = 1, model_type = "covariates"
      )
      ## ----end
      glmmTMB_13_mse
    }),
    tar_target(mse_2_mod_glmmTMB_13b_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_13 <- mod_glmmTMB_13b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      ## ---- glmmTMB_13b_mse_2
      glmmTMB_13_mse <- mse_glmmTMB(mod_glmmTMB_13,
        newdata = newdata_11b, type = 2, model_type = "covariates"
      )
      ## ----end
      glmmTMB_13_mse
    }),
    tar_target(mse_3_mod_glmmTMB_13b_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_13 <- mod_glmmTMB_13b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      ## ---- glmmTMB_13b_mse_3
      glmmTMB_13_mse <- mse_glmmTMB(mod_glmmTMB_13,
        newdata = newdata_10b, type = 3, model_type = "covariates"
      )
      ## ----end
      glmmTMB_13_mse
    }),
    
    ## brms -----------------------------------------------------------
    tar_target(mod_brms_13_, {
      benthos_fixed_locs_obs_13 <- incomplete_spatial_data_prep_13_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- brms_pre_13
      benthos_fixed_locs_obs_13 |>
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
      ## ---- brms_13
      mod_brms_13 <- brm(mod_form,
        data = benthos_fixed_locs_obs_13,
        iter = 5000,
        warmup = 1000,
        chains = 3,
        cores = 3,
        prior = priors,
        thin =  5,
        control = list(adapt_delta = 0.99),
        backend = "cmdstanr"
      )
      saveRDS(mod_brms_13,
        file = paste0(data_path, "synthetic/mod_brms_13.rds")
      ) 
      ## ----end
      mod_brms_13
    }),
    tar_target(pred_1_mod_brms_13_, {
      pred_brms <- pred_brms_
      mod_brms_13 <- mod_brms_13_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_13 <- incomplete_spatial_newdata_13_
      mod_simple_13 <- mod_simple_13_ #sampled_simple_raw_means_
      ## ---- brms_13_pred_1
      newdata <- newdata_13
      true_sum <- mod_simple_13
      brms_13_pred <- pred_brms(mod_brms_13,
        type = 1,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_13"
      )
      ## ----end
      brms_13_pred
    }),
    tar_target(pred_2_mod_brms_13_, {
      pred_brms <- pred_brms_
      mod_brms_13 <- mod_brms_13_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11 <- incomplete_spatial_newdata_11_
      benthos_reefs_temporal_summary_11 <- read_sampled_reefs_data_11_plot_2_
      ## ---- brms_13_pred_2
      newdata <- newdata_11
      true_sum <- benthos_reefs_temporal_summary_11
      brms_13_pred <- pred_brms(mod_brms_13,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_13"
      )
      ## ----end
      brms_13_pred
    }),
    tar_target(pred_3_mod_brms_13_, {
      pred_brms <- pred_brms_
      mod_brms_13 <- mod_brms_13_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10 <- incomplete_spatial_newdata_10_
      benthos_reefs_temporal_summary_10 <- read_sampled_reefs_data_10_plot_2_
      ## ---- brms_13_pred_3
      newdata <- newdata_10
      true_sum <- benthos_reefs_temporal_summary_10
      brms_13_pred <- pred_brms(mod_brms_13,
        type = 3,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_13"
      )
      ## ----end
      brms_13_pred
    }),
    tar_target(brms_trace_13_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_13 <- mod_brms_13_
      ## ---- brms_trace_13
      mod_brms_13 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_13.rds")
      )
      vars <- mod_brms_13 |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_13$fit |> stan_trace(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_trace_13.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ac_13_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_13 <- mod_brms_13_
      ## ---- brms_ac_13
      mod_brms_13 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_13.rds")
      )
      vars <- mod_brms_13 |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_13$fit |> stan_ac(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ac_13.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_rhat_13_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_13 <- mod_brms_13_
      ## ---- brms_rhat_13
      mod_brms_13 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_13.rds")
      )
      g <- mod_brms_13$fit |> stan_rhat()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_rhat_13.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ess_13_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_13 <- mod_brms_13_
      ## ---- brms_ess_13
      mod_brms_13 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_13.rds")
      )
      g <- mod_brms_13$fit |> stan_ess()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ess_13.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ppc_13_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_13 <- mod_brms_13_
      ## ---- brms_ppc_13
      mod_brms_13 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_13.rds")
      )
      g <- mod_brms_13 |> pp_check( type='dens_overlay', ndraws=100)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ppc_13.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(mse_1_mod_brms_13_, {
      mse_brms <- mse_brms_
      mod_brms_13 <- mod_brms_13_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_13b <- incomplete_spatial_newdata_13b_
      ## ---- brms_13_mse_1
      brms_13_mse <- mse_brms(mod_brms_13,
        newdata = newdata_13b, type = 1, model_type = ""
      )
      ## ----end
      brms_13_mse
    }),
    tar_target(mse_2_mod_brms_13_, {
      mse_brms <- mse_brms_
      mod_brms_13 <- mod_brms_13_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      ## ---- brms_13_mse_2
      brms_13_mse <- mse_brms(mod_brms_13,
        newdata = newdata_11b, type = 2, model_type = ""
      )
      ## ----end
      brms_13_mse
    }),
    tar_target(mse_3_mod_brms_13_, {
      mse_brms <- mse_brms_
      mod_brms_13 <- mod_brms_13_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      ## ---- brms_13_mse_3
      brms_13_mse <- mse_brms(mod_brms_13,
        newdata = newdata_10b, type = 3, model_type = ""
      )
      ## ----end
      brms_13_mse
    }),
    
    ## brms + covariates ----------------------------------------------
    tar_target(mod_brms_13b_, {
      benthos_fixed_locs_obs_13 <- incomplete_spatial_data_prep_13_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- brms_pre_13b
      benthos_fixed_locs_obs_13 |>
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
      mod_form <- bf(cover ~ fYear + scale(Longitude) + scale(Latitude) +
                       scale(CYC) + scale(DHW) + scale(OTHER) +
                       (1 | Site) + (1 | Transect),
        family = "Beta"
      )
      ## ----end
      ## ---- brms_13
      mod_brms_13b <- brm(mod_form,
        data = benthos_fixed_locs_obs_13,
        iter = 5000,
        warmup = 1000,
        chains = 3,
        cores = 3,
        prior = priors,
        thin =  5,
        control = list(adapt_delta = 0.99),
        backend = "cmdstanr"
      )
      saveRDS(mod_brms_13b,
        file = paste0(data_path, "synthetic/mod_brms_13b.rds")
      ) 
      ## ----end
      mod_brms_13b
    }),
    tar_target(pred_1_mod_brms_13b_, {
      pred_brms <- pred_brms_
      mod_brms_13b <- mod_brms_13b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_13b <- incomplete_spatial_newdata_13b_
      mod_simple_13 <- mod_simple_13_ #sampled_simple_raw_means_
      ## ---- brms_13b_pred_1
      true_sum <- mod_simple_13
      newdata <- newdata_13b
      brms_13b_pred <- pred_brms(mod_brms_13b,
        type = 1,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_13"
      )
      ## ----end
      brms_13b_pred
    }),
    tar_target(pred_2_mod_brms_13b_, {
      pred_brms <- pred_brms_
      mod_brms_13 <- mod_brms_13b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      benthos_reefs_temporal_summary_11 <- read_sampled_reefs_data_11_plot_2_
      ## ---- brms_13b_pred_2
      newdata <- newdata_11b
      true_sum <- benthos_reefs_temporal_summary_11
      brms_13b_pred <- pred_brms(mod_brms_13,
        type = 2,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_13"
      )
      ## ----end
      brms_13b_pred
    }),
    tar_target(pred_3_mod_brms_13b_, {
      pred_brms <- pred_brms_
      mod_brms_13 <- mod_brms_13b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      benthos_reefs_temporal_summary_10 <- read_sampled_reefs_data_10_plot_2_
      ## ---- brms_13b_pred_3
      newdata <- newdata_10b
      true_sum <- benthos_reefs_temporal_summary_10
      brms_13b_pred <- pred_brms(mod_brms_13,
        type = 3,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_13"
      )
      ## ----end
      brms_13b_pred
    }),
    tar_target(brms_trace_13b_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_13b <- mod_brms_13b_
      ## ---- brms_trace_13
      mod_brms_13b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_13b.rds")
      )
      vars <- mod_brms_13b |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_13b$fit |> stan_trace(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_trace_13b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ac_13b_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_13b <- mod_brms_13b_
      ## ---- brms_ac_13b
      mod_brms_13b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_13b.rds")
      )
      vars <- mod_brms_13b |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_13b$fit |> stan_ac(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ac_13b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_rhat_13b_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_13b <- mod_brms_13b_
      ## ---- brms_rhat_13b
      mod_brms_13b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_13b.rds")
      )
      g <- mod_brms_13b$fit |> stan_rhat()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_rhat_13b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ess_13b_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_13b <- mod_brms_13b_
      ## ---- brms_ess_13b
      mod_brms_13b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_13b.rds")
      )
      g <- mod_brms_13b$fit |> stan_ess()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ess_13b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ppc_13b_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_13b <- mod_brms_13b_
      ## ---- brms_ppc_13b
      mod_brms_13b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_13b.rds")
      )
      g <- mod_brms_13b |> pp_check( type='dens_overlay', ndraws=100)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ppc_13b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(mse_1_mod_brms_13b_, {
      mse_brms <- mse_brms_
      mod_brms_13 <- mod_brms_13b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_13b <- incomplete_spatial_newdata_13b_
      ## ---- brms_13b_mse_1
      brms_13_mse <- mse_brms(mod_brms_13,
        newdata = newdata_13b, type = 1, model_type = "covariates"
      )
      ## ----end
      brms_13_mse
    }),
    tar_target(mse_2_mod_brms_13b_, {
      mse_brms <- mse_brms_
      mod_brms_13 <- mod_brms_13b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      ## ---- brms_13b_mse_2
      brms_13_mse <- mse_brms(mod_brms_13,
        newdata = newdata_11b, type = 2, model_type = "covariates"
      )
      ## ----end
      brms_13_mse
    }),
    tar_target(mse_3_mod_brms_13b_, {
      mse_brms <- mse_brms_
      mod_brms_13 <- mod_brms_13b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      ## ---- brms_13b_mse_3
      brms_13_mse <- mse_brms(mod_brms_13,
        newdata = newdata_10b, type = 3, model_type = "covariates"
      )
      ## ----end
      brms_13_mse
    }),

    ## stan -----------------------------------------------------------
    tar_target(mod_stan_13_, {
      source("model_functions.R")
      benthos_fixed_locs_obs_13 <- incomplete_spatial_data_prep_13_
      data_path <- incomplete_spatial_global_parameters_$data_path
      site_extra_functions_
      ## ---- stan_pre_13
      benthos_fixed_locs_obs_13 <-
        benthos_fixed_locs_obs_13 |>
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
      saveRDS(benthos_fixed_locs_obs_13,
        file = paste0(data_path, "synthetic/saveRDS(benthos_fixed_locs_obs_13_forstan.rds")
      ) 
      stan_data <- prepare_data_for_stan(benthos_fixed_locs_obs_13, yrs = NULL)
      model_stan <- cmdstanr::cmdstan_model(stan_file = "model1.stan")
      ## ----end
      ## ---- stan_13
      mod_stan_13 <- model_stan$sample(
        data = stan_data,
        seed = 133,
        iter_sampling = 5000,
        iter_warmup = 1000,
        thin = 5,
        chains = 3,
        parallel_chains = 3,
        adapt_delta = 0.99,
        output_dir = paste0(data_path, "synthetic/"),
      )
      saveRDS(mod_stan_13,
        file = paste0(data_path, "synthetic/mod_stan_13.rds")
      ) 
      ## ----end
      mod_stan_13
    }),
    tar_target(pred_1_mod_stan_13_, {
      pred_stan <- pred_stan_
      mod_stan_13 <- mod_stan_13_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_13 <- incomplete_spatial_newdata_13_
      mod_simple_13 <- mod_simple_13_ #sampled_simple_raw_means_
      ## ---- stan_13_pred_1
      newdata <- newdata_13
      true_sum <- mod_simple_13
      stan_13_pred <- pred_stan(mod_stan_13,
        type = 1,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_stan_13"
      )
      ## ----end
      stan_13_pred
    }),
    tar_target(pred_2_mod_stan_13_, {
      pred_stan <- pred_stan_
      mod_stan_13 <- mod_stan_13_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11 <- incomplete_spatial_newdata_11_
      benthos_reefs_temporal_summary_11 <- read_sampled_reefs_data_11_plot_2_
      ## ---- stan_13_pred_2
      newdata <- newdata_11
      true_sum <- benthos_reefs_temporal_summary_11
      stan_13_pred <- pred_stan(mod_stan_13,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_stan_13"
      )
      ## ----end
      stan_13_pred
    }),
    tar_target(pred_3_mod_stan_13_, {
      pred_stan <- pred_stan_
      mod_stan_13 <- mod_stan_13_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10 <- incomplete_spatial_newdata_10_
      benthos_reefs_temporal_summary_10 <- read_sampled_reefs_data_10_plot_2_
      ## ---- stan_13_pred_3
      newdata <- newdata_10
      true_sum <- benthos_reefs_temporal_summary_10
      stan_13_pred <- pred_stan(mod_stan_13,
        type = 3,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_stan_13"
      )
      ## ----end
      stan_13_pred
    }),
    tar_target(stan_trace_13_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_stan_13 <- mod_stan_13_
      ## ---- stan_trace_13
      mod_stan_13 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_13.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_13$draws(variables = c("beta", "phi", "sd_1", "sd_2", "sd_3")) |>
        mcmc_trace() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_trace_13.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ac_13_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_stan_13 <- mod_stan_13_
      ## ---- stan_ac_13
      mod_stan_13 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_13.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_13$draws(variables = c("beta", "phi", "sd_1", "sd_2", "sd_3")) |>
        mcmc_acf() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ac_13.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_rhat_13_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_stan_13 <- mod_stan_13_
      ## ---- stan_rhat_13
      mod_stan_13 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_13.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_13 |> bayesplot::rhat() |> 
        mcmc_rhat_hist() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_rhat_13.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ess_13_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_stan_13 <- mod_stan_13_
      ## ---- stan_ess_13
      mod_stan_13 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_13.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_13 |> bayesplot::neff_ratio() |> 
        mcmc_neff_hist() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ess_13.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ppc_13_, {
      benthos_fixed_locs_obs_13 <- incomplete_spatial_data_prep_13_
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_stan_13 <- mod_stan_13_
      ## ---- stan_ppc_13
      mod_stan_13 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_13.rds")
      )
      g <- 
        bayesplot::pp_check(
          benthos_fixed_locs_obs_13$cover,
          mod_stan_13$draws("ypred", format = "matrix")[1:100, ],
          ppc_dens_overlay
        ) +
        theme_classic()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ppc_13.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(mse_1_mod_stan_13_, {
      mse_stan <- mse_stan_
      mod_stan_13 <- mod_stan_13_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_13b <- incomplete_spatial_newdata_13b_
      ## ---- stan_13_mse_1
      stan_13_mse <- mse_stan(mod_stan_13,
        newdata = newdata_13b, type = 1, model_type = ""
      )
      ## ----end
      stan_13_mse
    }),
    tar_target(mse_2_mod_stan_13_, {
      mse_stan <- mse_stan_
      mod_stan_13 <- mod_stan_13_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      ## ---- stan_13_mse_2
      stan_13_mse <- mse_stan(mod_stan_13,
        newdata = newdata_11b, type = 2, model_type = ""
      )
      ## ----end
      stan_13_mse
    }),
    tar_target(mse_3_mod_stan_13_, {
      mse_stan <- mse_stan_
      mod_stan_13 <- mod_stan_13_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      ## ---- stan_13_mse_3
      stan_13_mse <- mse_stan(mod_stan_13,
        newdata = newdata_10b, type = 3, model_type = ""
      )
      ## ----end
      stan_13_mse
    }),
    
    ## gbm ------------------------------------------------------------
    tar_target(mod_gbm_13_, {
      benthos_fixed_locs_obs_13 <- incomplete_spatial_data_prep_13_
      data_path <- incomplete_spatial_global_parameters_$data_path
      ## ---- gbm_13
      mod_gbm_13 <- gbm(cover ~ fYear,
        data =  benthos_fixed_locs_obs_13,
        distribution = "gaussian",
        n.trees = 10000,
        interaction.depth = 5,
        shrinkage = 0.001,
        bag.fraction = 0.5,
        cv.folds = 5,
        verbose = TRUE
      )
      saveRDS(mod_gbm_13,
        file = paste0(data_path, "synthetic/mod_gbm_13.rds")
      ) 
      ## ----end
      ## ---- gbm_post_13
      n.trees <- gbm.perf(mod_gbm_13, method = "cv")
      ## ----end
      list(mod_gbm_13 = mod_gbm_13, n.trees = n.trees)
    }),
    tar_target(pred_1_mod_gbm_13_, {
      pred_gbm <- pred_gbm_
      mod_gbm_13 <- mod_gbm_13_$mod_gbm_13
      n.trees <- mod_gbm_13_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_13 <- incomplete_spatial_newdata_13_
      mod_simple_13 <- mod_simple_13_ #sampled_simple_raw_means_
      ## ---- gbm_13_pred_1
      newdata <- newdata_13
      true_sum <- mod_simple_13
      gbm_13_pred <- pred_gbm(mod_gbm_13,
        n.trees = n.trees,
        type = 1,
        model_type = "",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_13"
      )
      ## ----end
      gbm_13_pred
    }),
    tar_target(pred_2_mod_gbm_13_, {
      pred_gbm <- pred_gbm_
      mod_gbm_13 <- mod_gbm_13_$mod_gbm_13
      n.trees <- mod_gbm_13_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11 <- incomplete_spatial_newdata_11_
      benthos_reefs_temporal_summary_11 <- read_sampled_reefs_data_11_plot_2_
      ## ---- gbm_13_pred_2
      newdata <- newdata_11
      true_sum <- benthos_reefs_temporal_summary_11
      gbm_13_pred <- pred_gbm(mod_gbm_13,
        n.trees = n.trees,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_13"
      )
      ## ----end
      gbm_13_pred
    }),
    tar_target(pred_3_mod_gbm_13_, {
      pred_gbm <- pred_gbm_
      mod_gbm_13 <- mod_gbm_13_$mod_gbm_13
      n.trees <- mod_gbm_13_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10 <- incomplete_spatial_newdata_10_
      benthos_reefs_temporal_summary_10 <- read_sampled_reefs_data_10_plot_2_
      ## ---- gbm_13_pred_3
      newdata <- newdata_10
      true_sum <- benthos_reefs_temporal_summary_10
      gbm_13_pred <- pred_gbm(mod_gbm_13,
        n.trees = n.trees,
        type = 3,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_13"
      )
      ## ----end
      gbm_13_pred
    }),
    tar_target(mse_1_mod_gbm_13_, {
      mse_gbm <- mse_gbm_
      mod_gbm_13 <- mod_gbm_13_$mod_gbm_13
      n.trees <- mod_gbm_13_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_13b <- incomplete_spatial_newdata_13b_
      ## ---- gbm_13_mse_1
      ## )
      gbm_13_mse <- mse_gbm(mod_gbm_13,
        n.trees = n.trees,
        newdata = newdata_13b, type = 1, model_type = ""
      )
      ## ----end
      gbm_13_mse
    }),
    tar_target(mse_2_mod_gbm_13_, {
      mse_gbm <- mse_gbm_
      mod_gbm_13 <- mod_gbm_13_$mod_gbm_13
      n.trees <- mod_gbm_13_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      ## ---- gbm_13_mse_2
      gbm_13_mse <- mse_gbm(mod_gbm_13,
        n.trees = n.trees,
        newdata = newdata_11b, type = 2, model_type = ""
      )
      ## ----end
      gbm_13_mse
    }),
    tar_target(mse_3_mod_gbm_13_, {
      mse_gbm <- mse_gbm_
      mod_gbm_13 <- mod_gbm_13_$mod_gbm_13
      n.trees <- mod_gbm_13_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      ## ---- gbm_13_mse_3
      gbm_13_mse <- mse_gbm(mod_gbm_13,
        n.trees = n.trees,
        newdata = newdata_10b, type = 3, model_type = ""
      )
      ## ----end
      gbm_13_mse
    }),

    ## gbm + covartiates ----------------------------------------------
    tar_target(mod_gbm_13b_, {
      benthos_fixed_locs_obs_13 <- incomplete_spatial_data_prep_13_
      data_path <- incomplete_spatial_global_parameters_$data_path
      ## ---- gbm_13b
      mod_gbm_13b <- gbm(cover ~ fYear + Latitude + Longitude + CYC + DHW + OTHER,
        data =  benthos_fixed_locs_obs_13,
        distribution = "gaussian",
        var.monotone = c(0, 0, 0, -1, -1, -1),
        n.trees = 10000,
        interaction.depth = 5,
        shrinkage = 0.001,
        bag.fraction = 0.5,
        cv.folds = 5,
        verbose = TRUE
      )
      saveRDS(mod_gbm_13b,
        file = paste0(data_path, "synthetic/mod_gbm_13b.rds")
      )
      ## ----end
      ## ---- gbm_post_13b
      n.trees <- gbm.perf(mod_gbm_13b, method = "cv")
      ## ----end
      list(mod_gbm_13b = mod_gbm_13b, n.trees = n.trees)
    }),
    tar_target(pred_1_mod_gbm_13b_, {
      pred_gbm <- pred_gbm_
      mod_gbm_13b <- mod_gbm_13b_$mod_gbm_13b
      n.trees <- mod_gbm_13b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_13b <- incomplete_spatial_newdata_13b_
      mod_simple_13 <- mod_simple_13_ #sampled_simple_raw_means_
      ## ---- gbm_13b_pred_1
      true_sum <- mod_simple_13
      newdata <- newdata_13b
      gbm_13_pred <- pred_gbm(mod_gbm_13b,
        n.trees = n.trees,
        type = 1,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_13"
      )
      ## ----end
      gbm_13_pred
    }),
    tar_target(pred_2_mod_gbm_13b_, {
      pred_gbm <- pred_gbm_
      mod_gbm_13b <- mod_gbm_13b_$mod_gbm_13b
      n.trees <- mod_gbm_13b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      benthos_reefs_temporal_summary_11 <- read_sampled_reefs_data_11_plot_2_
      ## ---- gbm_13b_pred_2
      newdata <- newdata_11b
      true_sum <- benthos_reefs_temporal_summary_11
      gbm_13_pred <- pred_gbm(mod_gbm_13b,
        n.trees = n.trees,
        type = 2,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_13"
      )
      ## ----end
      gbm_13_pred
    }),
    tar_target(pred_3_mod_gbm_13b_, {
      pred_gbm <- pred_gbm_
      mod_gbm_13b <- mod_gbm_13b_$mod_gbm_13
      n.trees <- mod_gbm_13b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      benthos_reefs_temporal_summary_10 <- read_sampled_reefs_data_10_plot_2_
      ## ---- gbm_13b_pred_3
      newdata <- newdata_10b
      true_sum <- benthos_reefs_temporal_summary_10
      newdata <- readRDS(
        file = paste0(data_path, "synthetic/newdata_13b.rds")
      )
      gbm_13_pred <- pred_gbm(mod_gbm_13b,
        n.trees = n.trees,
        type = 3,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_13"
      )
      ## ----end
      gbm_13_pred
    }),
    tar_target(infl_gbm_13b_, {
      mod_gbm_13b <- mod_gbm_13b_$mod_gbm_13b
      n.trees <- mod_gbm_13b_$n.trees
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      ## ---- gbm_infl_13b
      ## gbm_3b_infl <- gbm::relative.influence(mod_gbm_3b, n.trees = n.trees, scale = TRUE, sort = TRUE)
      infl <- summary(mod_gbm_13b, n.trees =  n.trees, plot = FALSE)
      g <-
        infl |>
        as.data.frame() |>
        arrange(rel.inf) |>
        mutate(var = factor(var, levels = unique(var))) |>
        ggplot(aes(y = var, x = rel.inf)) +
        geom_bar(stat = "identity") +
        scale_y_discrete("") +
        scale_x_continuous("Relative influence") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_infl_mod_gbm_13b.png"
        ),
        g,
        width = 6, height = 4, dpi = 72
      )
      ## ----end
    }),
    tar_target(mse_1_mod_gbm_13b_, {
      mse_gbm <- mse_gbm_
      mod_gbm_13 <- mod_gbm_13b_$mod_gbm_13b
      n.trees <- mod_gbm_13b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_13b <- incomplete_spatial_newdata_13b_
      ## ---- gbm_13b_mse_1
      gbm_13_mse <- mse_gbm(mod_gbm_13,
        n.trees = n.trees,
        newdata = newdata_13b, type = 1, model_type = "covariates"
      )
      ## ----end
      gbm_13_mse
    }),
    tar_target(mse_2_mod_gbm_13b_, {
      mse_gbm <- mse_gbm_
      mod_gbm_13 <- mod_gbm_13b_$mod_gbm_13b
      n.trees <- mod_gbm_13b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      ## ---- gbm_13b_mse_2
      gbm_13_mse <- mse_gbm(mod_gbm_13,
        n.trees = n.trees,
        newdata = newdata_11b, type = 2, model_type = "covariates"
      )
      ## ----end
      gbm_13_mse
    }),
    tar_target(mse_3_mod_gbm_13b_, {
      mse_gbm <- mse_gbm_
      mod_gbm_13 <- mod_gbm_13b_$mod_gbm_13b
      n.trees <- mod_gbm_13b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      ## ---- gbm_13b_mse_3
      gbm_13_mse <- mse_gbm(mod_gbm_13,
        n.trees = n.trees,
        newdata = newdata_10b, type = 3, model_type = "covariates"
      )
      ## ----end
      gbm_13_mse
    }),

    ## dbarts --------------------------------------------------------
    tar_target(mod_dbarts_13_, {
      benthos_fixed_locs_obs_13 <- incomplete_spatial_data_prep_13_
      benthos_fixed_locs_obs_10 <- incomplete_spatial_data_prep_10_
      benthos_fixed_locs_obs_11 <- incomplete_spatial_data_prep_11_
      data_path <- incomplete_spatial_global_parameters_$data_path
      newdata_13 <- incomplete_spatial_newdata_13_
      newdata_13b <- incomplete_spatial_newdata_13b_
      newdata_13_a <- incomplete_spatial_data_prep_13_#incomplete_spatial_newdata_13b_
      newdata_10 <- incomplete_spatial_newdata_10_
      newdata_10_a <- incomplete_spatial_newdata_10b_
      newdata_11 <- incomplete_spatial_newdata_11_
      newdata_11_a <- incomplete_spatial_newdata_11b_
      ## ---- dbarts_13
      print(head(benthos_fixed_locs_obs_13))
      mod_dbarts_13 <- bart2(log(cover) ~ fYear,
        data =   benthos_fixed_locs_obs_13,
        keepTrees = TRUE
      )
      saveRDS(mod_dbarts_13,
        file = paste0(data_path, "synthetic/mod_dbarts_13.rds")
      ) 
      ## ----end
      ## Unfortunately, the next part must be in the same tar_target
      ## due to the way dbarts stores pointers - they cannot be stored
      ## ---- dbarts_pred_13
      preds <- predict(mod_dbarts_13, newdata_13, type = "ev") |>
        exp()
      saveRDS(preds,
        file = paste0(data_path, "synthetic/mod_dbarts_13_preds.rds")
      ) 
      preds_sum <- preds |>
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_13_preds_sum.rds")
      ) 

      newdata <- newdata_13_a |>
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
      preds_a <- predict(mod_dbarts_13, newdata, type = "ev") |>
        exp() 
      saveRDS(preds_a,
        file = paste0(data_path, "synthetic/mod_dbarts_13_a_preds.rds")
      ) 

      newdata_13b <- newdata_13b |>
        mutate(fYear = factor(Year))
      preds_13b <- predict(mod_dbarts_13, newdata_13b, type = "ev") |>
        exp() 
      saveRDS(preds_13b,
        file = paste0(data_path, "synthetic/mod_dbarts_13b_preds.rds")
      ) 
      
      preds_10 <- predict(mod_dbarts_13, newdata_10, type = "ev") |>
        exp()
      saveRDS(preds_10,
        file = paste0(data_path, "synthetic/mod_dbarts_13_10_preds.rds")
      ) 
      preds_10_sum <- preds_10 |>
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_10_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_13_10_preds_sum.rds")
      ) 
      newdata_10_a <- newdata_10_a |>
        ungroup() |>
        mutate(fYear = factor(Year))
      preds_10_a <- predict(mod_dbarts_13, newdata_10_a, type = "ev") |>
        exp()
      saveRDS(preds_10_a,
        file = paste0(data_path, "synthetic/mod_dbarts_13_10_a_preds.rds")
      ) 

      preds_11 <- predict(mod_dbarts_13, newdata_11, type = "ev") |>
        exp()
      saveRDS(preds_11,
        file = paste0(data_path, "synthetic/mod_dbarts_13_11_preds.rds")
      ) 
      preds_11_sum <- preds_11 |> 
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_11_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_11_preds_sum.rds")
      ) 
      newdata_11_a <- newdata_11_a |>
        ungroup() |>
        mutate(fYear = factor(Year))
      preds_11_a <- predict(mod_dbarts_13, newdata_11_a, type = "ev") |>
        exp()
      saveRDS(preds_11_a,
        file = paste0(data_path, "synthetic/mod_dbarts_13_11_a_preds.rds")
      ) 
      ## ----end
      list(
        mod_dbarts_13 = mod_dbarts_13,
        preds = preds,
        preds_sum = preds_sum,
        preds_a = preds_a,
        preds_10 = preds_10,
        preds_10_sum = preds_10_sum,
        preds_10_a = preds_10_a,
        preds_11 = preds_11,
        preds_11_sum = preds_11_sum,
        preds_11_a = preds_11_a
      )
    }),
    tar_target(pred_1_mod_dbarts_13_, {
      pred_dbarts <- pred_dbarts_
      preds <- mod_dbarts_13_$preds_sum
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_13 <- incomplete_spatial_newdata_13_
      mod_simple_13 <- mod_simple_13_ #sampled_simple_raw_means_
      ## ---- dbarts_13_pred_1
      newdata <- newdata_13
      true_sum <- mod_simple_13
      dbarts_13_pred <- pred_dbarts(
        preds,
        type = 1,
        model_type = "",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_13"
      )
      ## ----end
      dbarts_13_pred
    }),
    tar_target(pred_2_mod_dbarts_13_, {
      pred_dbarts <- pred_dbarts_
      preds <- mod_dbarts_13_$preds_11_sum
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11 <- incomplete_spatial_newdata_11_
      benthos_reefs_temporal_summary_11 <- read_sampled_reefs_data_11_plot_2_
      ## ---- dbarts_13_pred_2
      newdata <- newdata_11
      true_sum <- benthos_reefs_temporal_summary_11
      dbarts_13_pred <- pred_dbarts(
        preds,
        type = 2,
        model_type = "",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_13"
      )
      ## ----end
      dbarts_13_pred
    }),
    tar_target(pred_3_mod_dbarts_13_, {
      pred_dbarts <- pred_dbarts_
      preds <- mod_dbarts_13_$preds_10_sum
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10 <- incomplete_spatial_newdata_10_
      benthos_reefs_temporal_summary_10 <- read_sampled_reefs_data_10_plot_2_
      ## ---- dbarts_13_pred_3
      newdata <- newdata_10
      true_sum <- benthos_reefs_temporal_summary_10
      dbarts_13_pred <- pred_dbarts(
        preds,
        type = 3,
        model_type = "",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_13"
      )
      ## ----end
      dbarts_13_pred
    }),
    tar_target(mse_1_mod_dbarts_13_, {
      mse_dbarts <- mse_dbarts_
      preds <- mod_dbarts_13_$preds_a
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_13b <- incomplete_spatial_newdata_13b_
      ## ---- dbarts_13_mse_1
      newdata_13b <- newdata_13b |>
        mutate(fYear = factor(Year)) |> 
        ungroup()
      dbarts_13_mse <- mse_dbarts(preds,
        newdata = newdata_13b, type = 1, model_type = ""
      )
      ## ----end
      dbarts_13_mse
    }),
    tar_target(mse_2_mod_dbarts_13_, {
      mse_dbarts <- mse_dbarts_
      preds <- mod_dbarts_13_$preds_11_a
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_ #incomplete_spatial_newdata_13b_
      ## ---- dbarts_13_mse_2
      dbarts_13_mse <- mse_dbarts(preds,
        newdata = newdata_11b, type = 2, model_type = ""
      )
      ## ----end
      dbarts_13_mse
    }),
    tar_target(mse_3_mod_dbarts_13_, {
      mse_dbarts <- mse_dbarts_
      preds <- mod_dbarts_13_$preds_10_a
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_ #incomplete_spatial_newdata_13b_
      ## ---- dbarts_13_mse_3
      dbarts_13_mse <- mse_dbarts(preds,
        newdata = newdata_10b, type = 3, model_type = ""
      )
      ## ----end
      dbarts_13_mse
    }),

    ## dbarts + covariates --------------------------------------------
    tar_target(mod_dbarts_13b_, {
      benthos_fixed_locs_obs_13 <- incomplete_spatial_data_prep_13_
      benthos_fixed_locs_obs_10 <- incomplete_spatial_data_prep_10_
      benthos_fixed_locs_obs_11 <- incomplete_spatial_data_prep_11_
      data_path <- incomplete_spatial_global_parameters_$data_path
      newdata_13b <- incomplete_spatial_newdata_13b_
      newdata_10b <- incomplete_spatial_newdata_10b_
      newdata_11b <- incomplete_spatial_newdata_11b_
      ## ---- dbarts_13b
      print(head(benthos_fixed_locs_obs_13))
      mod_dbarts_13b <- bart2(log(cover) ~ fYear + Latitude + Longitude + CYC + DHW + OTHER,
        data =   benthos_fixed_locs_obs_13,
        keepTrees = TRUE
      )
      saveRDS(mod_dbarts_13b,
        file = paste0(data_path, "synthetic/mod_dbarts_13b.rds")
      ) 
      ## ----end
      ## Unfortunately, the next part must be in the same tar_target
      ## due to the way dbarts stores pointers - they cannot be stored
      ## ---- dbarts_pred_13
      newdata_13b <- newdata_13b |>
        mutate(fYear = factor(Year))
      preds_13b <- predict(mod_dbarts_13b, newdata_13b, type = "ev") |>
        exp() 
      saveRDS(preds_13b,
        file = paste0(data_path, "synthetic/mod_dbarts_13b_preds_13b.rds")
      ) 
      preds_13b_sum <- preds_13b |> summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_13b_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_13b_preds_13b_sum.rds")
      ) 
      newdata_10b <- newdata_10b |>
        mutate(fYear = factor(Year))
      preds_10b <- predict(mod_dbarts_13b, newdata_10b, type = "ev") |>
        exp() 
      saveRDS(preds_10b,
        file = paste0(data_path, "synthetic/mod_dbarts_13b_preds_10b.rds")
      ) 
      preds_10b_sum <- preds_10b |> 
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_10b_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_13b_preds_10b_sum.rds")
      ) 
      newdata_11b <- newdata_11b |>
        mutate(fYear = factor(Year))
      preds_11b <- predict(mod_dbarts_13b, newdata_11b, type = "ev") |>
        exp()
      saveRDS(preds_11b,
        file = paste0(data_path, "synthetic/mod_dbarts_13b_preds_11b.rds")
      ) 
      preds_11b_sum <- preds_11b |>
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_11b_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_13b_preds_11b_sum.rds")
      ) 
      ## ----end
      list(
        mod_dbarts_13 = mod_dbarts_13b,
        preds_13b = preds_13b,
        preds_13b_sum = preds_13b_sum,
        preds_10b = preds_10b,
        preds_10b_sum = preds_10b_sum,
        preds_11b = preds_11b,
        preds_11b_sum = preds_11b_sum
      )
    }),
    tar_target(pred_1_mod_dbarts_13b_, {
      pred_dbarts <- pred_dbarts_
      preds_13b <- mod_dbarts_13b_$preds_13b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_13b <- incomplete_spatial_newdata_13b_
      mod_simple_13 <- mod_simple_13_ #sampled_simple_raw_means_
      ## ---- dbarts_13b_pred_1
      true_sum <- mod_simple_13
      newdata <- newdata_13b
      dbarts_13_pred <- pred_dbarts(
        preds_13b,
        type = 1,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_13"
      )
      ## ----end
      dbarts_13_pred
    }),
    tar_target(pred_2_mod_dbarts_13b_, {
      pred_dbarts <- pred_dbarts_
      preds_11b <- mod_dbarts_13b_$preds_11b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      benthos_reefs_temporal_summary_11 <- read_sampled_reefs_data_11_plot_2_
      ## ---- dbarts_13b_pred_2
      newdata <- newdata_11b
      true_sum <- benthos_reefs_temporal_summary_11
      dbarts_13_pred <- pred_dbarts(
        preds_11b, 
        type = 2, 
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_13"
      )
      ## ----end
      dbarts_13_pred
    }),
    tar_target(pred_3_mod_dbarts_13b_, {
      pred_dbarts <- pred_dbarts_
      preds_10b <- mod_dbarts_13b_$preds_10b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      benthos_reefs_temporal_summary_10 <- read_sampled_reefs_data_10_plot_2_
      ## ---- dbarts_13b_pred_3
      newdata <- newdata_10b
      true_sum <- benthos_reefs_temporal_summary_10
      dbarts_13_pred <- pred_dbarts(
        preds_10b,
        type = 3,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_13"
      )
      ## ----end
      dbarts_13_pred
    }),
    tar_target(mse_1_mod_dbarts_13b_, {
      mse_dbarts <- mse_dbarts_
      preds_13b <- mod_dbarts_13b_$preds_13b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_13b <- incomplete_spatial_newdata_13b_
      ## ---- dbarts_13b_mse_1
      dbarts_13_mse <- mse_dbarts(preds_13b,
        newdata = newdata_13b, type = 1, model_type = "covariates"
      )
      ## ----end
      dbarts_13_mse
    }),
    tar_target(mse_2_mod_dbarts_13b_, {
      mse_dbarts <- mse_dbarts_
      preds_11b <- mod_dbarts_13b_$preds_11b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      ## ---- dbarts_13b_mse_2
      dbarts_13_mse <- mse_dbarts(preds_11b,
        newdata = newdata_11b, type = 2, model_type = "covariates"
      )
      ## ----end
      dbarts_13_mse
    }),
    tar_target(mse_3_mod_dbarts_13b_, {
      mse_dbarts <- mse_dbarts_
      preds_10b <- mod_dbarts_13b_$preds_10b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      ## ---- dbarts_13b_mse_3
      dbarts_13_mse <- mse_dbarts(preds_10b,
        newdata = newdata_10b, type = 3, model_type = "covariates"
      )
      ## ----end
      dbarts_13_mse
    }),
    
    ## xgboost + covariates -------------------------------------------
    tar_target(mod_xgboost_13b_prep_, {
      benthos_fixed_locs_obs_13 <- incomplete_spatial_data_prep_13_
      data_path <- incomplete_spatial_global_parameters_$data_path
      ## ---- xgboost_13_prep
      data_train_13b <- benthos_fixed_locs_obs_13 |>
        dplyr::select(cover, Year, Latitude, Longitude, CYC, DHW, OTHER) 
      ## ----end
      data_train_13b
    }),
    tar_target(mod_xgboost_13b_tune_, {
      data_train_13b <- mod_xgboost_13b_prep_
      data_path <- incomplete_spatial_global_parameters_$data_path
      ## ---- xgboost_13_tune
      ## Define the recipe
      tune_recipe <- recipe(cover ~ ., data = data_train_13b) |> 
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
        resamples = vfold_cv(data_train_13b, v = 5),
        grid = tune_grid_values)
      ## Get best set of parameters
      model_hyperparams <-
        select_best(tuned_results, metric = "rmse") |> 
        select(-".config") |> 
        as_tibble() |>
        mutate(nb_training = nrow(data_train_13b),
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
    tar_target(mod_xgboost_13b_fit_, {
      data_train_13b <- mod_xgboost_13b_prep_
      data_path <- incomplete_spatial_global_parameters_$data_path
      tune_recipe = mod_xgboost_13b_tune_$tune_recipe
      tune_model = mod_xgboost_13b_tune_$tune_model
      tune_workflow = mod_xgboost_13b_tune_$tune_workflow
      tune_grid_values = mod_xgboost_13b_tune_$tune_grid_values
      tuned_results = mod_xgboost_13b_tune_$tuned_results
      model_hyperparams = mod_xgboost_13b_tune_$model_hyperparams
      ## ---- xgboost_13_fit
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
      final_fitted_13b <- tune_workflow |>
        fit(data_train_13b)
      ## ----end
      final_fitted_13b
    }),
    tar_target(pred_1_mod_xgboost_13b_, {
      pred_xgboost <- pred_xgboost_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_13b <- incomplete_spatial_newdata_13b_
      final_fitted_13b <- mod_xgboost_13b_fit_
      mod_simple_13 <- mod_simple_13_ #sampled_simple_raw_means_
      ## ---- xgboost_13b_pred_1
      true_sum <- mod_simple_13
      newdata <- newdata_13b
      xgboost_13b_pred <- pred_xgboost(
        final_fitted_13b,
        type = 1,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_xgboost_13"
      )
      ## ----end
      xgboost_13b_pred
    }),
    tar_target(mse_1_mod_xgboost_13b_, {
      mse_xgboost <- mse_xgboost_
      final_fitted_13b <- mod_xgboost_13b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_13b <- incomplete_spatial_newdata_13b_
      ## ---- xgboost_13b_mse_1
      xgboost_13_mse <- mse_xgboost(final_fitted_13b,
        newdata = newdata_13b, type = 1, model_type = "covariates"
      )
      ## ----end
      xgboost_13_mse
    }),
    tar_target(pred_2_mod_xgboost_13b_, {
      pred_xgboost <- pred_xgboost_
      final_fitted_13b <- mod_xgboost_13b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      benthos_reefs_temporal_summary_11 <- read_sampled_reefs_data_11_plot_2_
      ## ---- xgboost_13b_pred_2
      newdata <- newdata_11b
      true_sum <- benthos_reefs_temporal_summary_11
      xgboost_13_pred <- pred_xgboost(final_fitted_13b,
        type = 2,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_xgboost_13"
      )
      ## ----end
      xgboost_13_pred
    }),
    tar_target(mse_2_mod_xgboost_13b_, {
      mse_xgboost <- mse_xgboost_
      final_fitted_13b <- mod_xgboost_13b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_11b <- incomplete_spatial_newdata_11b_
      ## ---- xgboost_13b_mse_2
      xgboost_13_mse <- mse_xgboost(final_fitted_13b,
        newdata = newdata_11b, type = 2, model_type = "covariates"
      )
      ## ----end
      xgboost_13_mse
    }),
    tar_target(pred_3_mod_xgboost_13b_, {
      pred_xgboost <- pred_xgboost_
      final_fitted_13b <- mod_xgboost_13b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      benthos_reefs_temporal_summary_10 <- read_sampled_reefs_data_10_plot_2_
      ## ---- xgboost_13b_pred_3
      newdata <- newdata_10b
      true_sum <- benthos_reefs_temporal_summary_10
      xgboost_13_pred <- pred_xgboost(final_fitted_13b,
        type = 3,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_xgboost_13"
      )
      ## ----end
      xgboost_13_pred
    }),
    tar_target(mse_3_mod_xgboost_13b_, {
      mse_xgboost <- mse_xgboost_
      final_fitted_13b <- mod_xgboost_13b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_10b <- incomplete_spatial_newdata_10b_
      ## ---- xgboost_13b_mse_3
      xgboost_13_mse <- mse_xgboost(final_fitted_13b,
        newdata = newdata_10b, type = 3, model_type = "covariates"
      )
      ## ----end
      xgboost_13_mse
    }),

    ## Comparisons ----------------------------------------------------
    tar_target(mse_1_mod_pymc_barts_13b_file_,
      paste0(
        incomplete_spatial_global_parameters_$data_path,
        "modelled/pymc_bart_13b_mse_1.csv"
      ), format = "file"),
    tar_target(mse_1_mod_pymc_barts_13b_, {
      read_csv(file = mse_1_mod_pymc_barts_13b_file_) 
    }),
    tar_target(mse_2_mod_pymc_barts_13b_file_,
      paste0(
        incomplete_spatial_global_parameters_$data_path,
        "modelled/pymc_bart_13b_mse_2.csv"
      ), format = "file"),
    tar_target(mse_2_mod_pymc_barts_13b_, {
      read_csv(file = mse_2_mod_pymc_barts_13b_file_) 
    }),
    tar_target(mse_3_mod_pymc_barts_13b_file_,
      paste0(
        incomplete_spatial_global_parameters_$data_path,
        "modelled/pymc_bart_13b_mse_3.csv"
      ), format = "file"),
    tar_target(mse_3_mod_pymc_barts_13b_, {
      read_csv(file = mse_3_mod_pymc_barts_13b_file_) 
    }),
    tar_target(comparisons_13_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      mse_1_mod_glmmTMB_13 <- mse_1_mod_glmmTMB_13_
      mse_2_mod_glmmTMB_13 <- mse_2_mod_glmmTMB_13_
      mse_3_mod_glmmTMB_13 <- mse_3_mod_glmmTMB_13_
      mse_1_mod_glmmTMB_13b <- mse_1_mod_glmmTMB_13b_
      mse_2_mod_glmmTMB_13b <- mse_2_mod_glmmTMB_13b_
      mse_3_mod_glmmTMB_13b <- mse_3_mod_glmmTMB_13b_

      mse_1_mod_brms_13 <- mse_1_mod_brms_13_
      mse_2_mod_brms_13 <- mse_2_mod_brms_13_
      mse_3_mod_brms_13 <- mse_3_mod_brms_13_
      mse_1_mod_brms_13b <- mse_1_mod_brms_13b_
      mse_2_mod_brms_13b <- mse_2_mod_brms_13b_
      mse_3_mod_brms_13b <- mse_3_mod_brms_13b_

      mse_1_mod_stan_13 <- mse_1_mod_stan_13_
      mse_2_mod_stan_13 <- mse_2_mod_stan_13_
      mse_3_mod_stan_13 <- mse_3_mod_stan_13_

      mse_1_mod_gbm_13 <- mse_1_mod_gbm_13_
      mse_2_mod_gbm_13 <- mse_2_mod_gbm_13_
      mse_3_mod_gbm_13 <- mse_3_mod_gbm_13_
      mse_1_mod_gbm_13b <- mse_1_mod_gbm_13b_
      mse_2_mod_gbm_13b <- mse_2_mod_gbm_13b_
      mse_3_mod_gbm_13b <- mse_3_mod_gbm_13b_

      mse_1_mod_dbarts_13 <- mse_1_mod_dbarts_13_
      mse_2_mod_dbarts_13 <- mse_2_mod_dbarts_13_
      mse_3_mod_dbarts_13 <- mse_3_mod_dbarts_13_
      mse_1_mod_dbarts_13b <- mse_1_mod_dbarts_13b_
      mse_2_mod_dbarts_13b <- mse_2_mod_dbarts_13b_
      mse_3_mod_dbarts_13b <- mse_3_mod_dbarts_13b_

      mse_1_mod_xgboost_13b <- mse_1_mod_xgboost_13b_
      mse_2_mod_xgboost_13b <- mse_2_mod_xgboost_13b_
      mse_3_mod_xgboost_13b <- mse_3_mod_xgboost_13b_

      mse_1_mod_pymc_barts_13b <- mse_1_mod_pymc_barts_13b_
      mse_2_mod_pymc_barts_13b <- mse_2_mod_pymc_barts_13b_
      mse_3_mod_pymc_barts_13b <- mse_3_mod_pymc_barts_13b_
      ## ---- comparisons_13
      mse_1_mod_pymc_barts_13b <- read_csv(
        file = paste0(data_path, "modelled/pymc_bart_13b_mse_1.csv")
      )   
      mse_2_mod_pymc_barts_13b <- read_csv(
        file = paste0(data_path, "modelled/pymc_bart_13b_mse_2.csv")
      )   
      mse_3_mod_pymc_barts_13b <- read_csv(
        file = paste0(data_path, "modelled/pymc_bart_13b_mse_3.csv")
      )   
      comparisons_13 <- bind_rows(
        mse_1_mod_glmmTMB_13,
        mse_2_mod_glmmTMB_13,
        mse_3_mod_glmmTMB_13,
        mse_1_mod_glmmTMB_13b,
        mse_2_mod_glmmTMB_13b,
        mse_3_mod_glmmTMB_13b,

        mse_1_mod_brms_13,
        mse_2_mod_brms_13,
        mse_3_mod_brms_13,
        mse_1_mod_brms_13b,
        mse_2_mod_brms_13b,
        mse_3_mod_brms_13b,
        
        mse_1_mod_stan_13,
        mse_2_mod_stan_13,
        mse_3_mod_stan_13,
        
        mse_1_mod_gbm_13,
        mse_2_mod_gbm_13,
        mse_3_mod_gbm_13,
        mse_1_mod_gbm_13b,
        mse_2_mod_gbm_13b,
        mse_3_mod_gbm_13b,

        mse_1_mod_dbarts_13,
        mse_2_mod_dbarts_13,
        mse_3_mod_dbarts_13,
        mse_1_mod_dbarts_13b,
        mse_2_mod_dbarts_13b,
        mse_3_mod_dbarts_13b,

        mse_1_mod_xgboost_13b,
        mse_2_mod_xgboost_13b,
        mse_3_mod_xgboost_13b,

        mse_1_mod_pymc_barts_13b,
        mse_2_mod_pymc_barts_13b,
        mse_3_mod_pymc_barts_13b
        )
      saveRDS(comparisons_13,
        file = paste0(data_path, "synthetic/comparisons_13.rds")
      ) 
      comps_13 <- 
        comparisons_13 |>
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
          type == 2 ~ "Predicting Southern reefs",
          type == 3 ~ "Predicting Northern reefs"
        )) 
      saveRDS(comps_13,
        file = paste0(data_path, "synthetic/comps_13.rds")
      ) 
      ## ----end
      comps_13
    }),
    tar_target(comparisons_13_plots_, {
      comps_13 <- comparisons_13_
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- comparisons_13
      g <-
        comps_13 |> 
        filter(stat == "mean", metric == "mse") |>
        ggplot(aes(x = value, y = model)) +
        geom_segment(aes(xend = 0, yend = model), color = "black") +
        geom_point() +
        facet_grid(model_type~type, scales = "free_y") +
        scale_y_discrete("") +
        scale_x_continuous("MSE") +
        theme_bw() +
        ggtitle("Mean Square Error of models built on Southern reefs data (2 reefs)")

      ggsave(
        filename = paste0(
          fig_path, "mse_13_1.png"
        ),
        g,
        width = 10, height = 6, dpi = 100
      )
      g <- 
        comps_13 |> 
        filter(stat == "mean", metric == "acc") |>
        ggplot(aes(x = value, y = model)) +
        geom_segment(aes(xend = 0, yend = model), color = "black") +
        geom_point() +
        facet_grid(model_type~type, scales = "free_y") +
        scale_y_discrete("") +
        scale_x_continuous("Mean inaccuracy (%)", labels = function(x) sprintf("%0.1f%%", x*100)) +
        theme_bw() +
        ggtitle("Mean accuracy of models built on Southern reefs data (2 reefs)")

      ggsave(
        filename = paste0(
          fig_path, "mse_13_2.png"
        ),
        g,
        width = 10, height = 6, dpi = 100
      )
      ## ----end
      comps_13
    }),

    ## Western subdomain (Sampled reefs) 
    ## simple ---------------------------------------------------------
    tar_target(mod_simple_16_, {
      benthos_fixed_locs_obs_16 <- incomplete_spatial_data_prep_16_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- simple_16
      mod_simple_16 <- benthos_fixed_locs_obs_16 |>
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
      mod_simple_16 <-
        benthos_fixed_locs_obs_16 |>
        group_by(Year) |>
        summarise(
          Mean = mean(cover),
          Median = median(cover)
        ) 
      saveRDS(mod_simple_16,
        file = paste0(data_path, "synthetic/mod_simple_16.rds")
      ) 
      ## ----end
      mod_simple_16
    }),
    
    ## glmmTMB --------------------------------------------------------
    tar_target(mod_glmmTMB_16_, {
      benthos_fixed_locs_obs_16 <- incomplete_spatial_data_prep_16_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_16
      mod_glmmTMB_16 <- glmmTMB(cover ~ fYear + (1 | Site) + (1 | Transect),
        data = benthos_fixed_locs_obs_16,
        family = "beta_family"
      )
      saveRDS(mod_glmmTMB_16,
        file = paste0(data_path, "synthetic/mod_glmmTMB_16.rds")
      ) 
      ## ----end
      mod_glmmTMB_16
    }),

    tar_target(mod_glmmTMB_16_sample_data_, {
      benthos_fixed_locs_obs_16 <- incomplete_spatial_data_prep_16_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_16 sample data
      benthos_fixed_locs_obs_16 
      benthos_fixed_locs_obs_16 |>
        group_by(Year) |>
        summarise(
          Mean = mean(cover),
          Median = median(cover)
          )
      ## ----end
      benthos_fixed_locs_obs_16
    }),
    tar_target(mod_glmmTMB_16_sample_data_summary_14_, {
      benthos_reefs_temporal_summary_14 <- read_sampled_reefs_data_14_plot_2_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_16 sample data summary 14
      benthos_reefs_temporal_summary_14 <- readRDS(
        file = paste0(data_path, "synthetic/benthos_reefs_temporal_summary_14.rds")
      )
      benthos_reefs_temporal_summary_14
      ## ----end
      benthos_reefs_temporal_summary_14
    }),
    tar_target(mod_glmmTMB_16_sample_data_summary_15_, {
      benthos_reefs_temporal_summary_15 <- read_sampled_reefs_data_15_plot_2_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_16 sample data summary 15
      benthos_reefs_temporal_summary_15 <- readRDS(
        file = paste0(data_path, "synthetic/benthos_reefs_temporal_summary_15.rds")
      )
      benthos_reefs_temporal_summary_15
      ## ----end
      benthos_reefs_temporal_summary_15
    }),

    tar_target(mod_glmmTMB_16_newdata_16_, {
      newdata_16 <- incomplete_spatial_newdata_16_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_16 newdata 16
      newdata_16 
      ## ----end
      newdata_16
    }),
    tar_target(mod_glmmTMB_16_newdata_16b_, {
      newdata_16b <- incomplete_spatial_newdata_16b_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_16 newdata 16b
      newdata_16b 
      ## ----end
      newdata_16b
    }),
    tar_target(mod_glmmTMB_16_newdata_17_, {
      newdata_17 <- incomplete_spatial_newdata_17_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_16 newdata 17
      newdata_17 
      ## ----end
      newdata_17
    }),
    tar_target(mod_glmmTMB_16_newdata_17b_, {
      newdata_17b <- incomplete_spatial_newdata_17_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_16 newdata 17b
      newdata_17b 
      ## ----end
      newdata_17b
    }),
    tar_target(mod_glmmTMB_16_newdata_14_, {
      newdata_14 <- incomplete_spatial_newdata_14_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_16 newdata 14
      newdata_14 
      ## ----end
      newdata_14
    }),
    tar_target(mod_glmmTMB_16_newdata_14b_, {
      newdata_14b <- incomplete_spatial_newdata_14b_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_16 newdata 14b
      newdata_14b 
      ## ----end
      newdata_14b
    }),
    tar_target(mod_glmmTMB_16_newdata_15_, {
      newdata_15 <- incomplete_spatial_newdata_15_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_16 newdata 15
      newdata_15 
      ## ----end
      newdata_15
    }),
    tar_target(mod_glmmTMB_16_newdata_15b_, {
      newdata_15b <- incomplete_spatial_newdata_15b_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_16 newdata 15b
      newdata_15b 
      ## ----end
      newdata_15b
    }),
    
    tar_target(dharma_mod_glmmTMB_16_, {
      DHARMa_glmmTMB <- DHARMa_glmmTMB_
      mod_glmmTMB_16 <- mod_glmmTMB_16_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- glmmTMB_16_dharma
      glmmTMB_16_dharma <- DHARMa_glmmTMB(mod_glmmTMB_16,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_16")
      ## ----end
      glmmTMB_16_dharma
    }),

    tar_target(pred_1_mod_glmmTMB_16_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_16 <- mod_glmmTMB_16_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_16 <- incomplete_spatial_newdata_16_
      mod_simple_16 <- mod_simple_16_ #sampled_simple_raw_means_
      ## ---- glmmTMB_16_pred_1
      newdata <- newdata_16
      true_sum <- mod_simple_16
      glmmTMB_16_pred <- pred_glmmTMB(mod_glmmTMB_16,
        type = 1,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_16"
      )
      ## ----end
      glmmTMB_16_pred
    }),
    tar_target(pred_2_mod_glmmTMB_16_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_16 <- mod_glmmTMB_16_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14 <- incomplete_spatial_newdata_14_
      benthos_reefs_temporal_summary_14 <- read_sampled_reefs_data_14_plot_2_
      ## ---- glmmTMB_16_pred_2
      newdata <- newdata_14
      true_sum <- benthos_reefs_temporal_summary_14
      glmmTMB_16_pred <- pred_glmmTMB(mod_glmmTMB_16,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_16"
      )
      ## ----end
      glmmTMB_16_pred
    }),
    tar_target(pred_3_mod_glmmTMB_16_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_16 <- mod_glmmTMB_16_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15 <- incomplete_spatial_newdata_15_
      benthos_reefs_temporal_summary_15 <- read_sampled_reefs_data_15_plot_2_
      ## ---- glmmTMB_16_pred_3
      newdata <- newdata_15
      true_sum <- benthos_reefs_temporal_summary_15
      glmmTMB_16_pred <- pred_glmmTMB(mod_glmmTMB_16,
        type = 3,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_16"
      )
      ## ----end
      glmmTMB_16_pred
    }),
    tar_target(mse_1_mod_glmmTMB_16_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_16 <- mod_glmmTMB_16_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_16b <- incomplete_spatial_newdata_16b_
      ## ---- glmmTMB_16_mse_1
      glmmTMB_16_mse <- mse_glmmTMB(mod_glmmTMB_16,
        newdata = newdata_16b, type = 1, model_type = ""
      )
      ## ----end
      glmmTMB_16_mse
    }),
    tar_target(mse_2_mod_glmmTMB_16_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_16 <- mod_glmmTMB_16_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      ## ---- glmmTMB_16_mse_2
      glmmTMB_16_mse <- mse_glmmTMB(mod_glmmTMB_16,
        newdata = newdata_14b, type = 2, model_type = ""
      )
      ## ----end
      glmmTMB_16_mse
    }),
    tar_target(mse_3_mod_glmmTMB_16_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_16 <- mod_glmmTMB_16_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      ## ---- glmmTMB_16_mse_3
      glmmTMB_16_mse <- mse_glmmTMB(mod_glmmTMB_16,
        newdata = newdata_15b, type = 3, model_type = ""
      )
      ## ----end
      glmmTMB_16_mse
    }),
    
    ## glmmTMB + covariates -------------------------------------------
    tar_target(mod_glmmTMB_16b_, {
      benthos_fixed_locs_obs_16 <- incomplete_spatial_data_prep_16_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_16b
      mod_glmmTMB_16b <- glmmTMB(cover ~ fYear + CYC + DHW + OTHER + Latitude + Longitude +
                                   (1 | Site) + (1 | Transect),
        data = benthos_fixed_locs_obs_16,
        family = "beta_family"
      )
      saveRDS(mod_glmmTMB_16b,
        file = paste0(data_path, "synthetic/mod_glmmTMB_16b.rds")
      ) 
      ## ----end
      mod_glmmTMB_16b
    }),
    tar_target(dharma_mod_glmmTMB_16b_, {
      DHARMa_glmmTMB <- DHARMa_glmmTMB_
      mod_glmmTMB_16b <- mod_glmmTMB_16b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- glmmTMB_16b_dharma
      glmmTMB_16b_dharma <- DHARMa_glmmTMB(mod_glmmTMB_16b,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_16b")
      ## ----end
      glmmTMB_16b_dharma
    }),
    tar_target(pred_1_mod_glmmTMB_16b_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_16b <- mod_glmmTMB_16b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_16b <- incomplete_spatial_newdata_16b_
      mod_simple_16 <- mod_simple_16_ #sampled_simple_raw_means_
      ## ---- glmmTMB_16b_pred_1
      true_sum <- mod_simple_16
      newdata <- newdata_16b
      glmmTMB_16b_pred <- pred_glmmTMB(mod_glmmTMB_16b,
        type = 1,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_16"
      )
      ## ----end
      glmmTMB_16b_pred
    }),
    tar_target(pred_2_mod_glmmTMB_16b_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_16 <- mod_glmmTMB_16b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      benthos_reefs_temporal_summary_14 <- read_sampled_reefs_data_14_plot_2_
      ## ---- glmmTMB_16b_pred_2
      newdata <- newdata_14b
      true_sum <- benthos_reefs_temporal_summary_14
      glmmTMB_16b_pred <- pred_glmmTMB(mod_glmmTMB_16,
        type = 2,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_16"
      ) 
      ## ----end
      glmmTMB_16b_pred
    }),
    tar_target(pred_3_mod_glmmTMB_16b_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_16 <- mod_glmmTMB_16b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      benthos_reefs_temporal_summary_15 <- read_sampled_reefs_data_15_plot_2_
      ## ---- glmmTMB_16b_pred_3
      newdata <- newdata_15b
      true_sum <- benthos_reefs_temporal_summary_15
      glmmTMB_16b_pred <- pred_glmmTMB(mod_glmmTMB_16,
        type = 3,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_16"
      )
      ## ----end
      glmmTMB_16b_pred
    }),
    tar_target(mse_1_mod_glmmTMB_16b_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_16 <- mod_glmmTMB_16b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_16b <- incomplete_spatial_newdata_16b_
      ## ---- glmmTMB_16b_mse_1
      glmmTMB_16_mse <- mse_glmmTMB(mod_glmmTMB_16,
        newdata = newdata_16b, type = 1, model_type = "covariates"
      )
      ## ----end
      glmmTMB_16_mse
    }),
    tar_target(mse_2_mod_glmmTMB_16b_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_16 <- mod_glmmTMB_16b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      ## ---- glmmTMB_16b_mse_2
      glmmTMB_16_mse <- mse_glmmTMB(mod_glmmTMB_16,
        newdata = newdata_14b, type = 2, model_type = "covariates"
      )
      ## ----end
      glmmTMB_16_mse
    }),
    tar_target(mse_3_mod_glmmTMB_16b_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_16 <- mod_glmmTMB_16b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      ## ---- glmmTMB_16b_mse_3
      glmmTMB_16_mse <- mse_glmmTMB(mod_glmmTMB_16,
        newdata = newdata_15b, type = 3, model_type = "covariates"
      )
      ## ----end
      glmmTMB_16_mse
    }),

    ## brms -----------------------------------------------------------
    tar_target(mod_brms_16_, {
      benthos_fixed_locs_obs_16 <- incomplete_spatial_data_prep_16_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- brms_pre_16
      benthos_fixed_locs_obs_16 |>
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
      ## ---- brms_16
      mod_brms_16 <- brm(mod_form,
        data = benthos_fixed_locs_obs_16,
        iter = 5000,
        warmup = 1000,
        chains = 3,
        cores = 3,
        prior = priors,
        thin =  5,
        control = list(adapt_delta = 0.99),
        backend = "cmdstanr"
      )
      saveRDS(mod_brms_16,
        file = paste0(data_path, "synthetic/mod_brms_16.rds")
      ) 
      ## ----end
      mod_brms_16
    }),
    tar_target(pred_1_mod_brms_16_, {
      pred_brms <- pred_brms_
      mod_brms_16 <- mod_brms_16_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_16 <- incomplete_spatial_newdata_16_
      mod_simple_16 <- mod_simple_16_ #sampled_simple_raw_means_
      ## ---- brms_16_pred_1
      newdata <- newdata_16
      true_sum <- mod_simple_16
      brms_16_pred <- pred_brms(mod_brms_16,
        type = 1,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_16"
      )
      ## ----end
      brms_16_pred
    }),
    tar_target(pred_2_mod_brms_16_, {
      pred_brms <- pred_brms_
      mod_brms_16 <- mod_brms_16_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14 <- incomplete_spatial_newdata_14_
      benthos_reefs_temporal_summary_14 <- read_sampled_reefs_data_14_plot_2_
      ## ---- brms_16_pred_2
      newdata <- newdata_14
      true_sum <- benthos_reefs_temporal_summary_14
      brms_16_pred <- pred_brms(mod_brms_16,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_16"
      )
      ## ----end
      brms_16_pred
    }),
    tar_target(pred_3_mod_brms_16_, {
      pred_brms <- pred_brms_
      mod_brms_16 <- mod_brms_16_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15 <- incomplete_spatial_newdata_15_
      benthos_reefs_temporal_summary_15 <- read_sampled_reefs_data_15_plot_2_
      ## ---- brms_16_pred_3
      newdata <- newdata_15
      true_sum <- benthos_reefs_temporal_summary_15
      brms_16_pred <- pred_brms(mod_brms_16,
        type = 3,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_16"
      )
      ## ----end
      brms_16_pred
    }),
    tar_target(brms_trace_16_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_16 <- mod_brms_16_
      ## ---- brms_trace_16
      mod_brms_16 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_16.rds")
      )
      vars <- mod_brms_16 |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_16$fit |> stan_trace(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_trace_16.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ac_16_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_16 <- mod_brms_16_
      ## ---- brms_ac_16
      mod_brms_16 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_16.rds")
      )
      vars <- mod_brms_16 |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_16$fit |> stan_ac(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ac_16.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_rhat_16_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_16 <- mod_brms_16_
      ## ---- brms_rhat_16
      mod_brms_16 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_16.rds")
      )
      g <- mod_brms_16$fit |> stan_rhat()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_rhat_16.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ess_16_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_16 <- mod_brms_16_
      ## ---- brms_ess_16
      mod_brms_16 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_16.rds")
      )
      g <- mod_brms_16$fit |> stan_ess()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ess_16.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ppc_16_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_16 <- mod_brms_16_
      ## ---- brms_ppc_16
      mod_brms_16 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_16.rds")
      )
      g <- mod_brms_16 |> pp_check( type='dens_overlay', ndraws=140)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ppc_16.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(mse_1_mod_brms_16_, {
      mse_brms <- mse_brms_
      mod_brms_16 <- mod_brms_16_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_16b <- incomplete_spatial_newdata_16b_
      ## ---- brms_16_mse_1
      brms_16_mse <- mse_brms(mod_brms_16,
        newdata = newdata_16b, type = 1, model_type = ""
      )
      ## ----end
      brms_16_mse
    }),
    tar_target(mse_2_mod_brms_16_, {
      mse_brms <- mse_brms_
      mod_brms_16 <- mod_brms_16_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      ## ---- brms_16_mse_2
      brms_16_mse <- mse_brms(mod_brms_16,
        newdata = newdata_14b, type = 2, model_type = ""
      )
      ## ----end
      brms_16_mse
    }),
    tar_target(mse_3_mod_brms_16_, {
      mse_brms <- mse_brms_
      mod_brms_16 <- mod_brms_16_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      ## ---- brms_16_mse_3
      brms_16_mse <- mse_brms(mod_brms_16,
        newdata = newdata_15b, type = 3, model_type = ""
      )
      ## ----end
      brms_16_mse
    }),
    
    ## brms + covariates ----------------------------------------------
    tar_target(mod_brms_16b_, {
      benthos_fixed_locs_obs_16 <- incomplete_spatial_data_prep_16_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- brms_pre_16b
      benthos_fixed_locs_obs_16 |>
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
      mod_form <- bf(cover ~ fYear + scale(Longitude) + scale(Latitude) +
                       scale(CYC) + scale(DHW) + scale(OTHER) +
                       (1 | Site) + (1 | Transect),
        family = "Beta"
      )
      ## ----end
      ## ---- brms_16
      mod_brms_16b <- brm(mod_form,
        data = benthos_fixed_locs_obs_16,
        iter = 5000,
        warmup = 1000,
        chains = 3,
        cores = 3,
        prior = priors,
        thin =  5,
        control = list(adapt_delta = 0.99),
        backend = "cmdstanr"
      )
      saveRDS(mod_brms_16b,
        file = paste0(data_path, "synthetic/mod_brms_16b.rds")
      ) 
      ## ----end
      mod_brms_16b
    }),
    tar_target(pred_1_mod_brms_16b_, {
      pred_brms <- pred_brms_
      mod_brms_16b <- mod_brms_16b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_16b <- incomplete_spatial_newdata_16b_
      mod_simple_16 <- mod_simple_16_ #sampled_simple_raw_means_
      ## ---- brms_16b_pred_1
      true_sum <- mod_simple_16
      newdata <- newdata_16b
      brms_16b_pred <- pred_brms(mod_brms_16b,
        type = 1,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_16"
      )
      ## ----end
      brms_16b_pred
    }),
    tar_target(pred_2_mod_brms_16b_, {
      pred_brms <- pred_brms_
      mod_brms_16 <- mod_brms_16b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      benthos_reefs_temporal_summary_14 <- read_sampled_reefs_data_14_plot_2_
      ## ---- brms_16b_pred_2
      newdata <- newdata_14b
      true_sum <- benthos_reefs_temporal_summary_14
      brms_16b_pred <- pred_brms(mod_brms_16,
        type = 2,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_16"
      )
      ## ----end
      brms_16b_pred
    }),
    tar_target(pred_3_mod_brms_16b_, {
      pred_brms <- pred_brms_
      mod_brms_16 <- mod_brms_16b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      benthos_reefs_temporal_summary_15 <- read_sampled_reefs_data_15_plot_2_
      ## ---- brms_16b_pred_3
      newdata <- newdata_15b
      true_sum <- benthos_reefs_temporal_summary_15
      brms_16b_pred <- pred_brms(mod_brms_16,
        type = 3,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_16"
      )
      ## ----end
      brms_16b_pred
    }),
    tar_target(brms_trace_16b_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_16b <- mod_brms_16b_
      ## ---- brms_trace_16
      mod_brms_16b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_16b.rds")
      )
      vars <- mod_brms_16b |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_16b$fit |> stan_trace(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_trace_16b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ac_16b_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_16b <- mod_brms_16b_
      ## ---- brms_ac_16b
      mod_brms_16b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_16b.rds")
      )
      vars <- mod_brms_16b |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_16b$fit |> stan_ac(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ac_16b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_rhat_16b_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_16b <- mod_brms_16b_
      ## ---- brms_rhat_16b
      mod_brms_16b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_16b.rds")
      )
      g <- mod_brms_16b$fit |> stan_rhat()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_rhat_16b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ess_16b_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_16b <- mod_brms_16b_
      ## ---- brms_ess_16b
      mod_brms_16b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_16b.rds")
      )
      g <- mod_brms_16b$fit |> stan_ess()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ess_16b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ppc_16b_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_16b <- mod_brms_16b_
      ## ---- brms_ppc_16b
      mod_brms_16b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_16b.rds")
      )
      g <- mod_brms_16b |> pp_check( type='dens_overlay', ndraws=100)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ppc_16b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(mse_1_mod_brms_16b_, {
      mse_brms <- mse_brms_
      mod_brms_16 <- mod_brms_16b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_16b <- incomplete_spatial_newdata_16b_
      ## ---- brms_16b_mse_1
      brms_16_mse <- mse_brms(mod_brms_16,
        newdata = newdata_16b, type = 1, model_type = "covariates"
      )
      ## ----end
      brms_16_mse
    }),
    tar_target(mse_2_mod_brms_16b_, {
      mse_brms <- mse_brms_
      mod_brms_16 <- mod_brms_16b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      ## ---- brms_16b_mse_2
      brms_16_mse <- mse_brms(mod_brms_16,
        newdata = newdata_14b, type = 2, model_type = "covariates"
      )
      ## ----end
      brms_16_mse
    }),
    tar_target(mse_3_mod_brms_16b_, {
      mse_brms <- mse_brms_
      mod_brms_16 <- mod_brms_16b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      ## ---- brms_16b_mse_3
      brms_16_mse <- mse_brms(mod_brms_16,
        newdata = newdata_15b, type = 3, model_type = "covariates"
      )
      ## ----end
      brms_16_mse
    }),

    ## stan -----------------------------------------------------------
    tar_target(mod_stan_16_, {
      source("model_functions.R")
      benthos_fixed_locs_obs_16 <- incomplete_spatial_data_prep_16_
      data_path <- incomplete_spatial_global_parameters_$data_path
      site_extra_functions_
      ## ---- stan_pre_16
      benthos_fixed_locs_obs_16 <-
        benthos_fixed_locs_obs_16 |>
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
      saveRDS(benthos_fixed_locs_obs_16,
        file = paste0(data_path, "synthetic/saveRDS(benthos_fixed_locs_obs_16_forstan.rds")
      ) 
      stan_data <- prepare_data_for_stan(benthos_fixed_locs_obs_16, yrs = NULL)
      model_stan <- cmdstanr::cmdstan_model(stan_file = "model1.stan")
      ## ----end
      ## ---- stan_16
      mod_stan_16 <- model_stan$sample(
        data = stan_data,
        seed = 163,
        iter_sampling = 5000,
        iter_warmup = 1000,
        thin = 5,
        chains = 3,
        parallel_chains = 3,
        adapt_delta = 0.99,
        output_dir = paste0(data_path, "synthetic/"),
      )
      saveRDS(mod_stan_16,
        file = paste0(data_path, "synthetic/mod_stan_16.rds")
      ) 
      ## ----end
      mod_stan_16
    }),
    tar_target(pred_1_mod_stan_16_, {
      pred_stan <- pred_stan_
      mod_stan_16 <- mod_stan_16_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_16 <- incomplete_spatial_newdata_16_
      mod_simple_16 <- mod_simple_16_ #sampled_simple_raw_means_
      ## ---- stan_16_pred_1
      newdata <- newdata_16
      true_sum <- mod_simple_16
      stan_16_pred <- pred_stan(mod_stan_16,
        type = 1,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_stan_16"
      )
      ## ----end
      stan_16_pred
    }),
    tar_target(pred_2_mod_stan_16_, {
      pred_stan <- pred_stan_
      mod_stan_16 <- mod_stan_16_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14 <- incomplete_spatial_newdata_14_
      benthos_reefs_temporal_summary_14 <- read_sampled_reefs_data_14_plot_2_
      ## ---- stan_16_pred_2
      newdata <- newdata_14
      true_sum <- benthos_reefs_temporal_summary_14
      stan_16_pred <- pred_stan(mod_stan_16,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_stan_16"
      )
      ## ----end
      stan_16_pred
    }),
    tar_target(pred_3_mod_stan_16_, {
      pred_stan <- pred_stan_
      mod_stan_16 <- mod_stan_16_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15 <- incomplete_spatial_newdata_15_
      benthos_reefs_temporal_summary_15 <- read_sampled_reefs_data_15_plot_2_
      ## ---- stan_16_pred_3
      newdata <- newdata_15
      true_sum <- benthos_reefs_temporal_summary_15
      stan_16_pred <- pred_stan(mod_stan_16,
        type = 3,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_stan_16"
      )
      ## ----end
      stan_16_pred
    }),
    tar_target(stan_trace_16_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_stan_16 <- mod_stan_16_
      ## ---- stan_trace_16
      mod_stan_16 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_16.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_16$draws(variables = c("beta", "phi", "sd_1", "sd_2", "sd_3")) |>
        mcmc_trace() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_trace_16.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ac_16_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_stan_16 <- mod_stan_16_
      ## ---- stan_ac_16
      mod_stan_16 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_16.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_16$draws(variables = c("beta", "phi", "sd_1", "sd_2", "sd_3")) |>
        mcmc_acf() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ac_16.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_rhat_16_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_stan_16 <- mod_stan_16_
      ## ---- stan_rhat_16
      mod_stan_16 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_16.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_16 |> bayesplot::rhat() |> 
        mcmc_rhat_hist() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_rhat_16.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ess_16_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_stan_16 <- mod_stan_16_
      ## ---- stan_ess_16
      mod_stan_16 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_16.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_16 |> bayesplot::neff_ratio() |> 
        mcmc_neff_hist() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ess_16.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ppc_16_, {
      benthos_fixed_locs_obs_16 <- incomplete_spatial_data_prep_16_
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_stan_16 <- mod_stan_16_
      ## ---- stan_ppc_16
      mod_stan_16 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_16.rds")
      )
      g <- 
        bayesplot::pp_check(
          benthos_fixed_locs_obs_16$cover,
          mod_stan_16$draws("ypred", format = "matrix")[1:100, ],
          ppc_dens_overlay
        ) +
        theme_classic()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ppc_16.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(mse_1_mod_stan_16_, {
      mse_stan <- mse_stan_
      mod_stan_16 <- mod_stan_16_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_16b <- incomplete_spatial_newdata_16b_
      ## ---- stan_16_mse_1
      stan_16_mse <- mse_stan(mod_stan_16,
        newdata = newdata_16b, type = 1, model_type = ""
      )
      ## ----end
      stan_16_mse
    }),
    tar_target(mse_2_mod_stan_16_, {
      mse_stan <- mse_stan_
      mod_stan_16 <- mod_stan_16_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      ## ---- stan_16_mse_2
      stan_16_mse <- mse_stan(mod_stan_16,
        newdata = newdata_14b, type = 2, model_type = ""
      )
      ## ----end
      stan_16_mse
    }),
    tar_target(mse_3_mod_stan_16_, {
      mse_stan <- mse_stan_
      mod_stan_16 <- mod_stan_16_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      ## ---- stan_16_mse_3
      stan_16_mse <- mse_stan(mod_stan_16,
        newdata = newdata_15b, type = 3, model_type = ""
      )
      ## ----end
      stan_16_mse
    }),

    ## gbm ------------------------------------------------------------
    tar_target(mod_gbm_16_, {
      benthos_fixed_locs_obs_16 <- incomplete_spatial_data_prep_16_
      data_path <- incomplete_spatial_global_parameters_$data_path
      ## ---- gbm_16
      mod_gbm_16 <- gbm(cover ~ fYear,
        data =  benthos_fixed_locs_obs_16,
        distribution = "gaussian",
        n.trees = 10000,
        interaction.depth = 5,
        shrinkage = 0.001,
        bag.fraction = 0.5,
        cv.folds = 5,
        verbose = TRUE
      )
      saveRDS(mod_gbm_16,
        file = paste0(data_path, "synthetic/mod_gbm_16.rds")
      ) 
      ## ----end
      ## ---- gbm_post_16
      n.trees <- gbm.perf(mod_gbm_16, method = "cv")
      ## ----end
      list(mod_gbm_16 = mod_gbm_16, n.trees = n.trees)
    }),
    tar_target(pred_1_mod_gbm_16_, {
      pred_gbm <- pred_gbm_
      mod_gbm_16 <- mod_gbm_16_$mod_gbm_16
      n.trees <- mod_gbm_16_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_16 <- incomplete_spatial_newdata_16_
      mod_simple_16 <- mod_simple_16_ #sampled_simple_raw_means_
      ## ---- gbm_16_pred_1
      newdata <- newdata_16
      true_sum <- mod_simple_16
      gbm_16_pred <- pred_gbm(mod_gbm_16,
        n.trees = n.trees,
        type = 1,
        model_type = "",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_16"
      )
      ## ----end
      gbm_16_pred
    }),
    tar_target(pred_2_mod_gbm_16_, {
      pred_gbm <- pred_gbm_
      mod_gbm_16 <- mod_gbm_16_$mod_gbm_16
      n.trees <- mod_gbm_16_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14 <- incomplete_spatial_newdata_14_
      benthos_reefs_temporal_summary_14 <- read_sampled_reefs_data_14_plot_2_
      ## ---- gbm_16_pred_2
      newdata <- newdata_14
      true_sum <- benthos_reefs_temporal_summary_14
      gbm_16_pred <- pred_gbm(mod_gbm_16,
        n.trees = n.trees,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_16"
      )
      ## ----end
      gbm_16_pred
    }),
    tar_target(pred_3_mod_gbm_16_, {
      pred_gbm <- pred_gbm_
      mod_gbm_16 <- mod_gbm_16_$mod_gbm_16
      n.trees <- mod_gbm_16_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15 <- incomplete_spatial_newdata_15_
      benthos_reefs_temporal_summary_15 <- read_sampled_reefs_data_15_plot_2_
      ## ---- gbm_16_pred_3
      newdata <- newdata_15
      true_sum <- benthos_reefs_temporal_summary_15
      gbm_16_pred <- pred_gbm(mod_gbm_16,
        n.trees = n.trees,
        type = 3,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_16"
      )
      ## ----end
      gbm_16_pred
    }),
    tar_target(mse_1_mod_gbm_16_, {
      mse_gbm <- mse_gbm_
      mod_gbm_16 <- mod_gbm_16_$mod_gbm_16
      n.trees <- mod_gbm_16_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_16b <- incomplete_spatial_newdata_16b_
      ## ---- gbm_16_mse_1
      gbm_16_mse <- mse_gbm(mod_gbm_16,
        n.trees = n.trees,
        newdata = newdata_16b, type = 1, model_type = ""
      )
      ## ----end
      gbm_16_mse
    }),
    tar_target(mse_2_mod_gbm_16_, {
      mse_gbm <- mse_gbm_
      mod_gbm_16 <- mod_gbm_16_$mod_gbm_16
      n.trees <- mod_gbm_16_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      ## ---- gbm_16_mse_2
      gbm_16_mse <- mse_gbm(mod_gbm_16,
        n.trees = n.trees,
        newdata = newdata_14b, type = 2, model_type = ""
      )
      ## ----end
      gbm_16_mse
    }),
    tar_target(mse_3_mod_gbm_16_, {
      mse_gbm <- mse_gbm_
      mod_gbm_16 <- mod_gbm_16_$mod_gbm_16
      n.trees <- mod_gbm_16_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      ## ---- gbm_16_mse_3
      gbm_16_mse <- mse_gbm(mod_gbm_16,
        n.trees = n.trees,
        newdata = newdata_15b, type = 3, model_type = ""
      )
      ## ----end
      gbm_16_mse
    }),
    
    ## gbm + covartiates ----------------------------------------------
    tar_target(mod_gbm_16b_, {
      benthos_fixed_locs_obs_16 <- incomplete_spatial_data_prep_16_
      data_path <- incomplete_spatial_global_parameters_$data_path
      ## ---- gbm_16b
      mod_gbm_16b <- gbm(cover ~ fYear + Latitude + Longitude + CYC + DHW + OTHER,
        data =  benthos_fixed_locs_obs_16,
        distribution = "gaussian",
        var.monotone = c(0, 0, 0, -1, -1, -1),
        n.trees = 10000,
        interaction.depth = 5,
        shrinkage = 0.001,
        bag.fraction = 0.5,
        cv.folds = 5, n.cores = 1,
        verbose = TRUE
      )
      saveRDS(mod_gbm_16b,
        file = paste0(data_path, "synthetic/mod_gbm_16b.rds")
      )
      ## ----end
      ## ---- gbm_post_16b
      n.trees <- gbm.perf(mod_gbm_16b, method = "cv")
      ## ----end
      list(mod_gbm_16b = mod_gbm_16b, n.trees = n.trees)
    }),
    tar_target(pred_1_mod_gbm_16b_, {
      pred_gbm <- pred_gbm_
      mod_gbm_16b <- mod_gbm_16b_$mod_gbm_16
      n.trees <- mod_gbm_16b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_16b <- incomplete_spatial_newdata_16b_
      mod_simple_16 <- mod_simple_16_ #sampled_simple_raw_means_
      ## ---- gbm_16b_pred_1
      true_sum <- mod_simple_16
      newdata <- newdata_16b
      gbm_16_pred <- pred_gbm(mod_gbm_16b,
        n.trees = n.trees,
        type = 1,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_16"
      )
      ## ----end
      gbm_16_pred
    }),
    tar_target(pred_2_mod_gbm_16b_, {
      pred_gbm <- pred_gbm_
      mod_gbm_16b <- mod_gbm_16b_$mod_gbm_16
      n.trees <- mod_gbm_16b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      benthos_reefs_temporal_summary_14 <- read_sampled_reefs_data_14_plot_2_
      ## ---- gbm_16b_pred_2
      newdata <- newdata_14b
      true_sum <- benthos_reefs_temporal_summary_14
      gbm_16_pred <- pred_gbm(mod_gbm_16b,
        n.trees = n.trees,
        type = 2,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_16"
      )
      ## ----end
      gbm_16_pred
    }),
    tar_target(pred_3_mod_gbm_16b_, {
      pred_gbm <- pred_gbm_
      mod_gbm_16b <- mod_gbm_16b_$mod_gbm_16
      n.trees <- mod_gbm_16b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      benthos_reefs_temporal_summary_15 <- read_sampled_reefs_data_15_plot_2_
      ## ---- gbm_16b_pred_3
      newdata <- newdata_15b
      true_sum <- benthos_reefs_temporal_summary_15
      gbm_16_pred <- pred_gbm(mod_gbm_16b,
        n.trees = n.trees,
        type = 3,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_16"
      )
      ## ----end
      gbm_16_pred
    }),
    tar_target(infl_gbm_16b_, {
      mod_gbm_16b <- mod_gbm_16b_$mod_gbm_16b
      n.trees <- mod_gbm_16b_$n.trees
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      ## ---- gbm_infl_16b
      ## gbm_3b_infl <- gbm::relative.influence(mod_gbm_3b, n.trees = n.trees, scale = TRUE, sort = TRUE)
      infl <- summary(mod_gbm_16b, n.trees =  n.trees, plot = FALSE)
      g <-
        infl |>
        as.data.frame() |>
        arrange(rel.inf) |>
        mutate(var = factor(var, levels = unique(var))) |>
        ggplot(aes(y = var, x = rel.inf)) +
        geom_bar(stat = "identity") +
        scale_y_discrete("") +
        scale_x_continuous("Relative influence") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_infl_mod_gbm_16b.png"
        ),
        g,
        width = 6, height = 4, dpi = 72
      )
      ## ----end
    }),
    tar_target(mse_1_mod_gbm_16b_, {
      mse_gbm <- mse_gbm_
      mod_gbm_16 <- mod_gbm_16b_$mod_gbm_16b
      n.trees <- mod_gbm_16b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_16b <- incomplete_spatial_newdata_16b_
      ## ---- gbm_16b_mse_1
      gbm_16_mse <- mse_gbm(mod_gbm_16,
        n.trees = n.trees,
        newdata = newdata_16b, type = 1, model_type = "covariates"
      )
      ## ----end
      gbm_16_mse
    }),
    tar_target(mse_2_mod_gbm_16b_, {
      mse_gbm <- mse_gbm_
      mod_gbm_16 <- mod_gbm_16b_$mod_gbm_16b
      n.trees <- mod_gbm_16b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      ## ---- gbm_16b_mse_2
      gbm_16_mse <- mse_gbm(mod_gbm_16,
        n.trees = n.trees,
        newdata = newdata_14b, type = 2, model_type = "covariates"
      )
      ## ----end
      gbm_16_mse
    }),
    tar_target(mse_3_mod_gbm_16b_, {
      mse_gbm <- mse_gbm_
      mod_gbm_16 <- mod_gbm_16b_$mod_gbm_16b
      n.trees <- mod_gbm_16b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      ## ---- gbm_16b_mse_3
      gbm_16_mse <- mse_gbm(mod_gbm_16,
        n.trees = n.trees,
        newdata = newdata_15b, type = 3, model_type = "covariates"
      )
      ## ----end
      gbm_16_mse
    }),
    
    ## dbarts --------------------------------------------------------
    tar_target(mod_dbarts_16_, {
      benthos_fixed_locs_obs_16 <- incomplete_spatial_data_prep_16_
      benthos_fixed_locs_obs_14 <- incomplete_spatial_data_prep_14_
      benthos_fixed_locs_obs_15 <- incomplete_spatial_data_prep_15_
      data_path <- incomplete_spatial_global_parameters_$data_path
      newdata_16 <- incomplete_spatial_newdata_16_
      newdata_16b <- incomplete_spatial_newdata_16b_
      newdata_16_a <- incomplete_spatial_data_prep_16_#incomplete_spatial_newdata_16b_
      newdata_14 <- incomplete_spatial_newdata_14_
      newdata_14_a <- incomplete_spatial_newdata_14b_
      newdata_15 <- incomplete_spatial_newdata_15_
      newdata_15_a <- incomplete_spatial_newdata_15b_
      ## ---- dbarts_16
      print(head(benthos_fixed_locs_obs_16))
      mod_dbarts_16 <- bart2(log(cover) ~ fYear,
        data =   benthos_fixed_locs_obs_16,
        keepTrees = TRUE
      )
      saveRDS(mod_dbarts_16,
        file = paste0(data_path, "synthetic/mod_dbarts_16.rds")
      ) 
      ## ----end
      ## Unfortunately, the next part must be in the same tar_target
      ## due to the way dbarts stores pointers - they cannot be stored
      ## ---- dbarts_pred_16
      preds <- predict(mod_dbarts_16, newdata_16, type = "ev") |>
        exp()
      saveRDS(preds,
        file = paste0(data_path, "synthetic/mod_dbarts_16_preds.rds")
      ) 
      preds_sum <- preds |>
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_16_preds_sum.rds")
      ) 

      newdata <- newdata_16_a |>
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
      preds_a <- predict(mod_dbarts_16, newdata, type = "ev") |>
        exp() 
      saveRDS(preds_a,
        file = paste0(data_path, "synthetic/mod_dbarts_16_a_preds.rds")
      ) 

      newdata_16b <- newdata_16b |>
        mutate(fYear = factor(Year))
      preds_16b <- predict(mod_dbarts_16, newdata_16b, type = "ev") |>
        exp() 
      saveRDS(preds_16b,
        file = paste0(data_path, "synthetic/mod_dbarts_16b_preds.rds")
      ) 


      
      preds_14 <- predict(mod_dbarts_16, newdata_14, type = "ev") |>
        exp()
      saveRDS(preds_14,
        file = paste0(data_path, "synthetic/mod_dbarts_14_preds.rds")
      ) 
      preds_14_sum <- preds_14 |>
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_14_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_14_preds_sum.rds")
      ) 
      newdata_14_a <- newdata_14_a |>
        ungroup() |>
        mutate(fYear = factor(Year))
      preds_14_a <- predict(mod_dbarts_16, newdata_14_a, type = "ev") |>
        exp()
      saveRDS(preds_14_a,
        file = paste0(data_path, "synthetic/mod_dbarts_14_a_preds.rds")
      ) 

      preds_15 <- predict(mod_dbarts_16, newdata_15, type = "ev") |>
        exp()
      saveRDS(preds_15,
        file = paste0(data_path, "synthetic/mod_dbarts_15_preds.rds")
      ) 
      preds_15_sum <- preds_15 |> 
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_15_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_15_preds_sum.rds")
      ) 
      newdata_15_a <- newdata_15_a |>
        ungroup() |>
        mutate(fYear = factor(Year))
      preds_15_a <- predict(mod_dbarts_16, newdata_15_a, type = "ev") |>
        exp()
      saveRDS(preds_15_a,
        file = paste0(data_path, "synthetic/mod_dbarts_15_a_preds.rds")
      ) 
      ## ----end
      list(
        mod_dbarts_16 = mod_dbarts_16,
        preds = preds,
        preds_sum = preds_sum,
        preds_a = preds_a,
        preds_14 = preds_14,
        preds_14_sum = preds_14_sum,
        preds_14_a = preds_14_a,
        preds_15 = preds_15,
        preds_15_sum = preds_15_sum,
        preds_15_a = preds_15_a
      )
    }),
    tar_target(pred_1_mod_dbarts_16_, {
      pred_dbarts <- pred_dbarts_
      preds <- mod_dbarts_16_$preds_sum
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_16 <- incomplete_spatial_newdata_16_
      mod_simple_16 <- mod_simple_16_ #sampled_simple_raw_means_
      ## ---- dbarts_16_pred_1
      newdata <- newdata_16
      true_sum <- mod_simple_16
      dbarts_16_pred <- pred_dbarts(
        preds,
        type = 1,
        model_type = "",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_16"
      )
      ## ----end
      dbarts_16_pred
    }),
    tar_target(pred_2_mod_dbarts_16_, {
      pred_dbarts <- pred_dbarts_
      preds <- mod_dbarts_16_$preds_14_sum
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14 <- incomplete_spatial_newdata_14_
      benthos_reefs_temporal_summary_14 <- read_sampled_reefs_data_14_plot_2_
      ## ---- dbarts_16_pred_2
      newdata <- newdata_14
      true_sum <- benthos_reefs_temporal_summary_14
      dbarts_16_pred <- pred_dbarts(
        preds,
        type = 2,
        model_type = "",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_16"
      )
      ## ----end
      dbarts_16_pred
    }),
    tar_target(pred_3_mod_dbarts_16_, {
      pred_dbarts <- pred_dbarts_
      preds <- mod_dbarts_16_$preds_15_sum
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15 <- incomplete_spatial_newdata_15_
      benthos_reefs_temporal_summary_15 <- read_sampled_reefs_data_15_plot_2_
      ## ---- dbarts_16_pred_3
      newdata <- newdata_15
      true_sum <- benthos_reefs_temporal_summary_15
      dbarts_16_pred <- pred_dbarts(
        preds,
        type = 3,
        model_type = "",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_16"
      )
      ## ----end
      dbarts_16_pred
    }),
    tar_target(mse_1_mod_dbarts_16_, {
      mse_dbarts <- mse_dbarts_
      preds <- mod_dbarts_16_$preds_a
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_16b <- incomplete_spatial_newdata_16b_
      ## ---- dbarts_16_mse_1
      newdata_16b <- newdata_16b |>
        mutate(fYear = factor(Year)) |> 
        ungroup()
      dbarts_16_mse <- mse_dbarts(preds,
        newdata = newdata_16b, type = 1, model_type = ""
      )
      ## ----end
      dbarts_16_mse
    }),
    tar_target(mse_2_mod_dbarts_16_, {
      mse_dbarts <- mse_dbarts_
      preds <- mod_dbarts_16_$preds_14_a
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_ #incomplete_spatial_newdata_16b_
      ## ---- dbarts_16_mse_2
      dbarts_16_mse <- mse_dbarts(preds,
        newdata = newdata_14b, type = 2, model_type = ""
      )
      ## ----end
      dbarts_16_mse
    }),
    tar_target(mse_3_mod_dbarts_16_, {
      mse_dbarts <- mse_dbarts_
      preds <- mod_dbarts_16_$preds_15_a
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_ #incomplete_spatial_newdata_16b_
      ## ---- dbarts_16_mse_3
      ## newdata <- newdata |>
      ##   ungroup() |>
      ##   mutate(fYear = factor(Year))
      dbarts_16_mse <- mse_dbarts(preds,
        newdata = newdata_15b, type = 3, model_type = ""
      )
      ## ----end
      dbarts_16_mse
    }),

    ## dbarts + covariates --------------------------------------------
    tar_target(mod_dbarts_16b_, {
      benthos_fixed_locs_obs_16 <- incomplete_spatial_data_prep_16_
      benthos_fixed_locs_obs_14 <- incomplete_spatial_data_prep_14_
      benthos_fixed_locs_obs_15 <- incomplete_spatial_data_prep_15_
      data_path <- incomplete_spatial_global_parameters_$data_path
      newdata_16b <- incomplete_spatial_newdata_16b_
      ## newdata_14 <- incomplete_spatial_newdata_14_
      newdata_14b <- incomplete_spatial_newdata_14b_
      newdata_15b <- incomplete_spatial_newdata_15b_
      ## ---- dbarts_16b
      print(head(benthos_fixed_locs_obs_16))
      mod_dbarts_16b <- bart2(log(cover) ~ fYear + Latitude + Longitude + CYC + DHW + OTHER,
        data =   benthos_fixed_locs_obs_16,
        keepTrees = TRUE
      )
      saveRDS(mod_dbarts_16b,
        file = paste0(data_path, "synthetic/mod_dbarts_16b.rds")
      ) 
      ## ----end
      ## Unfortunately, the next part must be in the same tar_target
      ## due to the way dbarts stores pointers - they cannot be stored
      ## ---- dbarts_pred_16b
      newdata_16b <- newdata_16b |>
        mutate(fYear = factor(Year))
      preds_16b <- predict(mod_dbarts_16b, newdata_16b, type = "ev") |>
        exp() 
      saveRDS(preds_16b,
        file = paste0(data_path, "synthetic/mod_dbarts_16b_preds_16b.rds")
      ) 
      preds_16b_sum <- preds_16b |> summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_16b_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_16b_preds_16b_sum.rds")
      )
      
      ## newdata_14 <- newdata_14 |>
      ##   mutate(fYear = factor(Year))
      ## preds_14 <- predict(mod_dbarts_16b, newdata_14, type = "ev") |>
      ##   exp() 
      ## saveRDS(preds_14,
      ##   file = paste0(data_path, "synthetic/mod_dbarts_14b_preds.rds")
      ## ) 
      ## preds_14_sum <- preds_14 |> 
      ##   summarise_draws(median, HDInterval::hdi)
      ## saveRDS(preds_14_sum,
      ##   file = paste0(data_path, "synthetic/mod_dbarts_14b_preds_sum.rds")
      ## ) 

      newdata_14b <- newdata_14b |>
        mutate(fYear = factor(Year))
      preds_14b <- predict(mod_dbarts_16b, newdata_14b, type = "ev") |>
        exp() 
      saveRDS(preds_14b,
        file = paste0(data_path, "synthetic/mod_dbarts_16b_preds_14b.rds")
      ) 
      preds_14b_sum <- preds_14b |> 
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_14b_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_16b_preds_14b_sum.rds")
      ) 

      newdata_15b <- newdata_15b |>
        mutate(fYear = factor(Year))
      preds_15b <- predict(mod_dbarts_16b, newdata_15b, type = "ev") |>
        exp()
      saveRDS(preds_15b,
        file = paste0(data_path, "synthetic/mod_dbarts_16b_preds_15b.rds")
      ) 
      preds_15b_sum <- preds_15b |>
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_15b_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_16b_preds_15b_sum.rds")
      ) 
      ## ----end
      list(
        mod_dbarts_16 = mod_dbarts_16b,
        preds_16b = preds_16b,
        preds_16b_sum = preds_16b_sum,
        preds_14b = preds_14b,
        preds_14b_sum = preds_14b_sum,
        preds_15b = preds_15b,
        preds_15b_sum = preds_15b_sum
      )
    }),
    tar_target(pred_1_mod_dbarts_16b_, {
      pred_dbarts <- pred_dbarts_
      preds_16b <- mod_dbarts_16b_$preds_16b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_16b <- incomplete_spatial_newdata_16b_
      mod_simple_16 <- mod_simple_16_ #sampled_simple_raw_means_
      ## ---- dbarts_16b_pred_1
      true_sum <- mod_simple_16
      newdata <- newdata_16b
      dbarts_16_pred <- pred_dbarts(
        preds_16b,
        type = 1,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_16"
      )
      ## ----end
      dbarts_16_pred
    }),
    tar_target(pred_2_mod_dbarts_16b_, {
      pred_dbarts <- pred_dbarts_
      preds_14b <- mod_dbarts_16b_$preds_14b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      benthos_reefs_temporal_summary_14 <- read_sampled_reefs_data_14_plot_2_
      ## ---- dbarts_16b_pred_2
      newdata <- newdata_14b
      true_sum <- benthos_reefs_temporal_summary_14
      dbarts_16_pred <- pred_dbarts(
        preds_14b, 
        type = 2, 
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_16"
      )
      ## ----end
      dbarts_16_pred
    }),
    tar_target(pred_3_mod_dbarts_16b_, {
      pred_dbarts <- pred_dbarts_
      preds_15b <- mod_dbarts_16b_$preds_15b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      benthos_reefs_temporal_summary_15 <- read_sampled_reefs_data_15_plot_2_
      ## ---- dbarts_16b_pred_3
      newdata <- newdata_15b
      true_sum <- benthos_reefs_temporal_summary_15
      dbarts_16_pred <- pred_dbarts(
        preds_15b,
        type = 3,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_16"
      )
      ## ----end
      dbarts_16_pred
    }),
    tar_target(mse_1_mod_dbarts_16b_, {
      mse_dbarts <- mse_dbarts_
      preds_16b <- mod_dbarts_16b_$preds_16b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata <- incomplete_spatial_data_prep_16_ #incomplete_spatial_newdata_16b_ #incomplete_spatial_newdata_16b_
      newdata_16b <- incomplete_spatial_newdata_16b_
      ## ---- dbarts_16b_mse_1
      dbarts_16_mse <- mse_dbarts(preds_16b,
        newdata = newdata_16b, type = 1, model_type = "covariates"
      )
      ## ----end
      dbarts_16_mse
    }),
    tar_target(mse_2_mod_dbarts_16b_, {
      mse_dbarts <- mse_dbarts_
      preds_14b <- mod_dbarts_16b_$preds_14b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      ## ---- dbarts_16b_mse_2
      dbarts_16_mse <- mse_dbarts(preds_14b,
        newdata = newdata_14b, type = 2, model_type = "covariates"
      )
      ## ----end
      dbarts_16_mse
    }),
    tar_target(mse_3_mod_dbarts_16b_, {
      mse_dbarts <- mse_dbarts_
      preds_15b <- mod_dbarts_16b_$preds_15b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      ## ---- dbarts_16b_mse_3
      dbarts_16_mse <- mse_dbarts(preds_15b,
        newdata = newdata_15b, type = 3, model_type = "covariates"
      )
      ## ----end
      dbarts_16_mse
    }),
    
    ## xgboost + covariates -------------------------------------------
    tar_target(mod_xgboost_16b_prep_, {
      benthos_fixed_locs_obs_16 <- incomplete_spatial_data_prep_16_
      data_path <- incomplete_spatial_global_parameters_$data_path
      ## ---- xgboost_16_prep
      data_train_16b <- benthos_fixed_locs_obs_16 |>
        dplyr::select(cover, Year, Latitude, Longitude, CYC, DHW, OTHER) 
      ## ----end
      data_train_16b
    }),
    tar_target(mod_xgboost_16b_tune_, {
      data_train_16b <- mod_xgboost_16b_prep_
      data_path <- incomplete_spatial_global_parameters_$data_path
      ## ---- xgboost_16_tune
      ## Define the recipe
      tune_recipe <- recipe(cover ~ ., data = data_train_16b) |> 
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
        resamples = vfold_cv(data_train_16b, v = 5),
        grid = tune_grid_values)
      ## Get best set of parameters
      model_hyperparams <-
        select_best(tuned_results, metric = "rmse") |> 
        select(-".config") |> 
        as_tibble() |>
        mutate(nb_training = nrow(data_train_16b),
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
    tar_target(mod_xgboost_16b_fit_, {
      data_train_16b <- mod_xgboost_16b_prep_
      data_path <- incomplete_spatial_global_parameters_$data_path
      tune_recipe = mod_xgboost_16b_tune_$tune_recipe
      tune_model = mod_xgboost_16b_tune_$tune_model
      tune_workflow = mod_xgboost_16b_tune_$tune_workflow
      tune_grid_values = mod_xgboost_16b_tune_$tune_grid_values
      tuned_results = mod_xgboost_16b_tune_$tuned_results
      model_hyperparams = mod_xgboost_16b_tune_$model_hyperparams
      ## ---- xgboost_16_fit
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
      final_fitted_16b <- tune_workflow |>
        fit(data_train_16b)
      ## ----end
      final_fitted_16b
    }),
    tar_target(pred_1_mod_xgboost_16b_, {
      pred_xgboost <- pred_xgboost_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_16b <- incomplete_spatial_newdata_16b_
      final_fitted_16b <- mod_xgboost_16b_fit_
      mod_simple_16 <- mod_simple_16_ #sampled_simple_raw_means_
      ## ---- xgboost_16b_pred_1
      true_sum <- mod_simple_16
      newdata <- newdata_16b
      xgboost_16b_pred <- pred_xgboost(
        final_fitted_16b,
        type = 1,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_xgboost_16"
      )
      ## ----end
      xgboost_16b_pred
    }),
    tar_target(mse_1_mod_xgboost_16b_, {
      mse_xgboost <- mse_xgboost_
      final_fitted_16b <- mod_xgboost_16b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_16b <- incomplete_spatial_newdata_16b_
      ## ---- xgboost_16b_mse_1
      xgboost_16_mse <- mse_xgboost(final_fitted_16b,
        newdata = newdata_16b, type = 1, model_type = "covariates"
      )
      ## ----end
      xgboost_16_mse
    }),
    tar_target(pred_2_mod_xgboost_16b_, {
      pred_xgboost <- pred_xgboost_
      final_fitted_16b <- mod_xgboost_16b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      benthos_reefs_temporal_summary_14 <- read_sampled_reefs_data_14_plot_2_
      ## ---- xgboost_16b_pred_2
      newdata <- newdata_14b
      true_sum <- benthos_reefs_temporal_summary_14
      xgboost_16_pred <- pred_xgboost(final_fitted_16b,
        type = 2,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_xgboost_16"
      )
      ## ----end
      xgboost_16_pred
    }),
    tar_target(mse_2_mod_xgboost_16b_, {
      mse_xgboost <- mse_xgboost_
      final_fitted_16b <- mod_xgboost_16b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      ## ---- xgboost_16b_mse_2
      xgboost_16_mse <- mse_xgboost(final_fitted_16b,
        newdata = newdata_14b, type = 2, model_type = "covariates"
      )
      ## ----end
      xgboost_16_mse
    }),
    tar_target(pred_3_mod_xgboost_16b_, {
      pred_xgboost <- pred_xgboost_
      final_fitted_16b <- mod_xgboost_16b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      benthos_reefs_temporal_summary_15 <- read_sampled_reefs_data_15_plot_2_
      ## ---- xgboost_16b_pred_3
      newdata <- newdata_15b
      true_sum <- benthos_reefs_temporal_summary_15
      xgboost_16_pred <- pred_xgboost(final_fitted_16b,
        type = 3,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_xgboost_16"
      )
      ## ----end
      xgboost_16_pred
    }),
    tar_target(mse_3_mod_xgboost_16b_, {
      mse_xgboost <- mse_xgboost_
      final_fitted_16b <- mod_xgboost_16b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      ## ---- xgboost_16b_mse_3
      xgboost_16_mse <- mse_xgboost(final_fitted_16b,
        newdata = newdata_15b, type = 3, model_type = "covariates"
      )
      ## ----end
      xgboost_16_mse
    }),

    ## Comparisons ----------------------------------------------------
    tar_target(mse_1_mod_pymc_barts_16b_file_,
      paste0(
        incomplete_spatial_global_parameters_$data_path,
        "modelled/pymc_bart_16b_mse_1.csv"
      ), format = "file"),
    tar_target(mse_1_mod_pymc_barts_16b_, {
      read_csv(file = mse_1_mod_pymc_barts_16b_file_) 
    }),
    tar_target(mse_2_mod_pymc_barts_16b_file_,
      paste0(
        incomplete_spatial_global_parameters_$data_path,
        "modelled/pymc_bart_16b_mse_2.csv"
      ), format = "file"),
    tar_target(mse_2_mod_pymc_barts_16b_, {
      read_csv(file = mse_2_mod_pymc_barts_16b_file_) 
    }),
    tar_target(mse_3_mod_pymc_barts_16b_file_,
      paste0(
        incomplete_spatial_global_parameters_$data_path,
        "modelled/pymc_bart_16b_mse_3.csv"
      ), format = "file"),
    tar_target(mse_3_mod_pymc_barts_16b_, {
      read_csv(file = mse_3_mod_pymc_barts_16b_file_) 
    }),
    tar_target(comparisons_16_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      mse_1_mod_glmmTMB_16 <- mse_1_mod_glmmTMB_16_
      mse_2_mod_glmmTMB_16 <- mse_2_mod_glmmTMB_16_
      mse_3_mod_glmmTMB_16 <- mse_3_mod_glmmTMB_16_
      mse_1_mod_glmmTMB_16b <- mse_1_mod_glmmTMB_16b_
      mse_2_mod_glmmTMB_16b <- mse_2_mod_glmmTMB_16b_
      mse_3_mod_glmmTMB_16b <- mse_3_mod_glmmTMB_16b_

      mse_1_mod_brms_16 <- mse_1_mod_brms_16_
      mse_2_mod_brms_16 <- mse_2_mod_brms_16_
      mse_3_mod_brms_16 <- mse_3_mod_brms_16_
      mse_1_mod_brms_16b <- mse_1_mod_brms_16b_
      mse_2_mod_brms_16b <- mse_2_mod_brms_16b_
      mse_3_mod_brms_16b <- mse_3_mod_brms_16b_

      mse_1_mod_stan_16 <- mse_1_mod_stan_16_
      mse_2_mod_stan_16 <- mse_2_mod_stan_16_
      mse_3_mod_stan_16 <- mse_3_mod_stan_16_

      mse_1_mod_gbm_16 <- mse_1_mod_gbm_16_
      mse_2_mod_gbm_16 <- mse_2_mod_gbm_16_
      mse_3_mod_gbm_16 <- mse_3_mod_gbm_16_
      mse_1_mod_gbm_16b <- mse_1_mod_gbm_16b_
      mse_2_mod_gbm_16b <- mse_2_mod_gbm_16b_
      mse_3_mod_gbm_16b <- mse_3_mod_gbm_16b_

      mse_1_mod_dbarts_16 <- mse_1_mod_dbarts_16_
      mse_2_mod_dbarts_16 <- mse_2_mod_dbarts_16_
      mse_3_mod_dbarts_16 <- mse_3_mod_dbarts_16_
      mse_1_mod_dbarts_16b <- mse_1_mod_dbarts_16b_
      mse_2_mod_dbarts_16b <- mse_2_mod_dbarts_16b_
      mse_3_mod_dbarts_16b <- mse_3_mod_dbarts_16b_

      mse_1_mod_xgboost_16b <- mse_1_mod_xgboost_16b_
      mse_2_mod_xgboost_16b <- mse_2_mod_xgboost_16b_
      mse_3_mod_xgboost_16b <- mse_3_mod_xgboost_16b_

      mse_1_mod_pymc_barts_16b <- mse_1_mod_pymc_barts_16b_
      mse_2_mod_pymc_barts_16b <- mse_2_mod_pymc_barts_16b_
      mse_3_mod_pymc_barts_16b <- mse_3_mod_pymc_barts_16b_
      ## ---- comparisons_16
      mse_1_mod_pymc_barts_16b <- read_csv(
        file = paste0(data_path, "modelled/pymc_bart_16b_mse_1.csv")
      )   
      mse_2_mod_pymc_barts_16b <- read_csv(
        file = paste0(data_path, "modelled/pymc_bart_16b_mse_2.csv")
      )   
      mse_3_mod_pymc_barts_16b <- read_csv(
        file = paste0(data_path, "modelled/pymc_bart_16b_mse_3.csv")
      )   
      comparisons_16 <- bind_rows(
        mse_1_mod_glmmTMB_16,
        mse_2_mod_glmmTMB_16,
        mse_3_mod_glmmTMB_16,
        mse_1_mod_glmmTMB_16b,
        mse_2_mod_glmmTMB_16b,
        mse_3_mod_glmmTMB_16b,

        mse_1_mod_brms_16,
        mse_2_mod_brms_16,
        mse_3_mod_brms_16,
        mse_1_mod_brms_16b,
        mse_2_mod_brms_16b,
        mse_3_mod_brms_16b,
        
        mse_1_mod_stan_16,
        mse_2_mod_stan_16,
        mse_3_mod_stan_16,
        
        mse_1_mod_gbm_16,
        mse_2_mod_gbm_16,
        mse_3_mod_gbm_16,
        mse_1_mod_gbm_16b,
        mse_2_mod_gbm_16b,
        mse_3_mod_gbm_16b,

        mse_1_mod_dbarts_16,
        mse_2_mod_dbarts_16,
        mse_3_mod_dbarts_16,
        mse_1_mod_dbarts_16b,
        mse_2_mod_dbarts_16b,
        mse_3_mod_dbarts_16b,

        mse_1_mod_xgboost_16b,
        mse_2_mod_xgboost_16b,
        mse_3_mod_xgboost_16b,

        mse_1_mod_pymc_barts_16b,
        mse_2_mod_pymc_barts_16b,
        mse_3_mod_pymc_barts_16b
        )
      saveRDS(comparisons_16,
        file = paste0(data_path, "synthetic/comparisons_16.rds")
      ) 
       
      comps_16 <- 
        comparisons_16 |>
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
          type == 2 ~ "Predicting Western reefs",
          type == 3 ~ "Predicting Eastern reefs"
        )) 
      saveRDS(comps_16,
        file = paste0(data_path, "synthetic/comps_16.rds")
      )  
      ## ----end
      comps_16
    }),
    tar_target(comparisons_16_plots_, {
      comps_16 <- comparisons_16_
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- comparisons_16_tab
      g <-
        comps_16 |> 
        filter(stat == "mean", metric == "mse") |>
        ggplot(aes(x = value, y = model)) +
        geom_segment(aes(xend = 0, yend = model), color = "black") +
        geom_point() +
        facet_grid(model_type~type, scales = "free_y") +
        scale_y_discrete("") +
        scale_x_continuous("MSE") +
        theme_bw() +
        ggtitle("Mean Square Error of models built on Western reefs data (25 reefs)")
 
      ggsave(
        filename = paste0(
          fig_path, "mse_16_1.png"
        ),
        g,
        width = 10, height = 6, dpi = 100
      )
      g <- 
        comps_16 |> 
        filter(stat == "mean", metric == "acc") |>
        ggplot(aes(x = value, y = model)) +
        geom_segment(aes(xend = 0, yend = model), color = "black") +
        geom_point() +
        facet_grid(model_type~type, scales = "free_y") +
        scale_y_discrete("") +
        scale_x_continuous("Mean inaccuracy (%)", labels = function(x) sprintf("%0.1f%%", x*100)) +
        theme_bw() +
        ggtitle("Mean inaccuracy (%) of models built on Western reefs data (2 reefs)")

      ggsave(
        filename = paste0(
          fig_path, "mse_16_2.png"
        ),
        g,
        width = 15, height = 6, dpi = 100
      )
      ## ----end
      comps_16
    }),

    ## Eastern subdomain (Sampled reefs) 
    ## simple ---------------------------------------------------------
    tar_target(mod_simple_17_, {
      benthos_fixed_locs_obs_17 <- incomplete_spatial_data_prep_17_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- simple_17
      mod_simple_17 <- benthos_fixed_locs_obs_17 |>
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
      mod_simple_17 <-
        benthos_fixed_locs_obs_17 |>
        group_by(Year) |>
        summarise(
          Mean = mean(cover),
          Median = median(cover)
        ) 
      saveRDS(mod_simple_17,
        file = paste0(data_path, "synthetic/mod_simple_17.rds")
      ) 
      ## ----end
      mod_simple_17
    }),
    
    ## glmmTMB --------------------------------------------------------
    tar_target(mod_glmmTMB_17_, {
      benthos_fixed_locs_obs_17 <- incomplete_spatial_data_prep_17_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_17
      mod_glmmTMB_17 <- glmmTMB(cover ~ fYear + (1 | Site) + (1 | Transect),
        data = benthos_fixed_locs_obs_17,
        family = "beta_family"
      )
      saveRDS(mod_glmmTMB_17,
        file = paste0(data_path, "synthetic/mod_glmmTMB_17.rds")
      ) 
      ## ----end
      mod_glmmTMB_17
    }),

    tar_target(mod_glmmTMB_17_sample_data_, {
      benthos_fixed_locs_obs_17 <- incomplete_spatial_data_prep_17_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_17 sample data
      benthos_fixed_locs_obs_17 
      benthos_fixed_locs_obs_17 |>
        group_by(Year) |>
        summarise(
          Mean = mean(cover),
          Median = median(cover)
          )
      ## ----end
      benthos_fixed_locs_obs_17
    }),
    tar_target(mod_glmmTMB_17_sample_data_summary_15_, {
      benthos_reefs_temporal_summary_15 <- read_sampled_reefs_data_15_plot_2_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_17 sample data summary 15
      benthos_reefs_temporal_summary_15 <- readRDS(
        file = paste0(data_path, "synthetic/benthos_reefs_temporal_summary_15.rds")
      )
      benthos_reefs_temporal_summary_15
      ## ----end
      benthos_reefs_temporal_summary_15
    }),
    tar_target(mod_glmmTMB_17_sample_data_summary_14_, {
      benthos_reefs_temporal_summary_14 <- read_sampled_reefs_data_14_plot_2_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_17 sample data summary 14
      benthos_reefs_temporal_summary_14 <- readRDS(
        file = paste0(data_path, "synthetic/benthos_reefs_temporal_summary_14.rds")
      )
      benthos_reefs_temporal_summary_14
      ## ----end
      benthos_reefs_temporal_summary_14
    }),
    
    tar_target(mod_glmmTMB_17_newdata_17_, {
      newdata_17 <- incomplete_spatial_newdata_17_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_17 newdata 17
      newdata_17 
      ## ----end
      newdata_17
    }),
    tar_target(mod_glmmTMB_17_newdata_17b_, {
      newdata_17b <- incomplete_spatial_newdata_17b_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_17 newdata 17b
      newdata_17b 
      ## ----end
      newdata_17b
    }),
    tar_target(mod_glmmTMB_17_newdata_16_, {
      newdata_16 <- incomplete_spatial_newdata_16_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_17 newdata 16
      newdata_16 
      ## ----end
      newdata_16
    }),
    tar_target(mod_glmmTMB_17_newdata_16b_, {
      newdata_16b <- incomplete_spatial_newdata_16_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_17 newdata 16b
      newdata_16b 
      ## ----end
      newdata_16b
    }),
    tar_target(mod_glmmTMB_17_newdata_15_, {
      newdata_15 <- incomplete_spatial_newdata_15_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_17 newdata 15
      newdata_15 
      ## ----end
      newdata_15
    }),
    tar_target(mod_glmmTMB_17_newdata_15b_, {
      newdata_15b <- incomplete_spatial_newdata_15b_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_17 newdata 15b
      newdata_15b 
      ## ----end
      newdata_15b
    }),
    tar_target(mod_glmmTMB_17_newdata_14_, {
      newdata_14 <- incomplete_spatial_newdata_14_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_17 newdata 14
      newdata_14 
      ## ----end
      newdata_14
    }),
    tar_target(mod_glmmTMB_17_newdata_14b_, {
      newdata_14b <- incomplete_spatial_newdata_14b_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_17 newdata 14b
      newdata_14b 
      ## ----end
      newdata_14b
    }),
    
    tar_target(dharma_mod_glmmTMB_17_, {
      DHARMa_glmmTMB <- DHARMa_glmmTMB_
      mod_glmmTMB_17 <- mod_glmmTMB_17_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- glmmTMB_17_dharma
      glmmTMB_17_dharma <- DHARMa_glmmTMB(mod_glmmTMB_17,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_17")
      ## ----end
      glmmTMB_17_dharma
    }),

    tar_target(pred_1_mod_glmmTMB_17_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_17 <- mod_glmmTMB_17_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_17 <- incomplete_spatial_newdata_17_
      mod_simple_17 <- mod_simple_17_ #sampled_simple_raw_means_
      ## ---- glmmTMB_17_pred_1
      newdata <- newdata_17
      true_sum <- mod_simple_17
      glmmTMB_17_pred <- pred_glmmTMB(mod_glmmTMB_17,
        type = 1,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_17"
      )
      ## ----end
      glmmTMB_17_pred
    }),
    tar_target(pred_2_mod_glmmTMB_17_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_17 <- mod_glmmTMB_17_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15 <- incomplete_spatial_newdata_15_
      benthos_reefs_temporal_summary_15 <- read_sampled_reefs_data_15_plot_2_
      ## ---- glmmTMB_17_pred_2
      newdata <- newdata_15
      true_sum <- benthos_reefs_temporal_summary_15
      glmmTMB_17_pred <- pred_glmmTMB(mod_glmmTMB_17,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_17"
      )
      ## ----end
      glmmTMB_17_pred
    }),
    tar_target(pred_3_mod_glmmTMB_17_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_17 <- mod_glmmTMB_17_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14 <- incomplete_spatial_newdata_14_
      benthos_reefs_temporal_summary_14 <- read_sampled_reefs_data_14_plot_2_
      ## ---- glmmTMB_17_pred_3
      newdata <- newdata_14
      true_sum <- benthos_reefs_temporal_summary_14
      glmmTMB_17_pred <- pred_glmmTMB(mod_glmmTMB_17,
        type = 3,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_17"
      )
      ## ----end
      glmmTMB_17_pred
    }),
    tar_target(mse_1_mod_glmmTMB_17_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_17 <- mod_glmmTMB_17_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_17b <- incomplete_spatial_newdata_17b_
      ## ---- glmmTMB_17_mse_1
      glmmTMB_17_mse <- mse_glmmTMB(mod_glmmTMB_17,
        newdata = newdata_17b, type = 1, model_type = ""
      )
      ## ----end
      glmmTMB_17_mse
    }),
    tar_target(mse_2_mod_glmmTMB_17_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_17 <- mod_glmmTMB_17_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      ## ---- glmmTMB_17_mse_2
      glmmTMB_17_mse <- mse_glmmTMB(mod_glmmTMB_17,
        newdata = newdata_15b, type = 2, model_type = ""
      )
      ## ----end
      glmmTMB_17_mse
    }),
    tar_target(mse_3_mod_glmmTMB_17_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_17 <- mod_glmmTMB_17_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      ## ---- glmmTMB_17_mse_3
      glmmTMB_17_mse <- mse_glmmTMB(mod_glmmTMB_17,
        newdata = newdata_14b, type = 3, model_type = ""
      )
      ## ----end
      glmmTMB_17_mse
    }),
    
    ## glmmTMB + covariates -------------------------------------------
    tar_target(mod_glmmTMB_17b_, {
      benthos_fixed_locs_obs_17 <- incomplete_spatial_data_prep_17_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_17b
      mod_glmmTMB_17b <- glmmTMB(cover ~ fYear + CYC + DHW + OTHER + Latitude + Longitude +
                                   (1 | Site) + (1 | Transect),
        data = benthos_fixed_locs_obs_17,
        family = "beta_family"
      )
      saveRDS(mod_glmmTMB_17b,
        file = paste0(data_path, "synthetic/mod_glmmTMB_17b.rds")
      ) 
      ## ----end
      mod_glmmTMB_17b
    }),
    tar_target(dharma_mod_glmmTMB_17b_, {
      DHARMa_glmmTMB <- DHARMa_glmmTMB_
      mod_glmmTMB_17b <- mod_glmmTMB_17b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- glmmTMB_17b_dharma
      glmmTMB_17b_dharma <- DHARMa_glmmTMB(mod_glmmTMB_17b,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_17b")
      ## ----end
      glmmTMB_17b_dharma
    }),
    tar_target(pred_1_mod_glmmTMB_17b_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_17b <- mod_glmmTMB_17b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_17b <- incomplete_spatial_newdata_17b_
      mod_simple_17 <- mod_simple_17_ #sampled_simple_raw_means_
      ## ---- glmmTMB_17b_pred_1
      true_sum <- mod_simple_17
      newdata <- newdata_17b
      glmmTMB_17b_pred <- pred_glmmTMB(mod_glmmTMB_17b,
        type = 1,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_17"
      )
      ## ----end
      glmmTMB_17b_pred
    }),
    tar_target(pred_2_mod_glmmTMB_17b_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_17 <- mod_glmmTMB_17b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      benthos_reefs_temporal_summary_15 <- read_sampled_reefs_data_15_plot_2_
      ## ---- glmmTMB_17b_pred_2
      newdata <- newdata_15b
      true_sum <- benthos_reefs_temporal_summary_15
      glmmTMB_17b_pred <- pred_glmmTMB(mod_glmmTMB_17,
        type = 2,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_17"
      ) 
      ## ----end
      glmmTMB_17b_pred
    }),
    tar_target(pred_3_mod_glmmTMB_17b_, {
      pred_glmmTMB <- pred_glmmTMB_
      mod_glmmTMB_17 <- mod_glmmTMB_17b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      benthos_reefs_temporal_summary_14 <- read_sampled_reefs_data_14_plot_2_
      ## ---- glmmTMB_17b_pred_3
      newdata <- newdata_14b
      true_sum <- benthos_reefs_temporal_summary_14
      glmmTMB_17b_pred <- pred_glmmTMB(mod_glmmTMB_17,
        type = 3,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_glmmTMB_17"
      )
      ## ----end
      glmmTMB_17b_pred
    }),
    tar_target(mse_1_mod_glmmTMB_17b_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_17 <- mod_glmmTMB_17b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_17b <- incomplete_spatial_newdata_17b_
      ## ---- glmmTMB_17b_mse_1
      glmmTMB_17_mse <- mse_glmmTMB(mod_glmmTMB_17,
        newdata = newdata_17b, type = 1, model_type = "covariates"
      )
      ## ----end
      glmmTMB_17_mse
    }),
    tar_target(mse_2_mod_glmmTMB_17b_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_17 <- mod_glmmTMB_17b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      ## ---- glmmTMB_17b_mse_2
      glmmTMB_17_mse <- mse_glmmTMB(mod_glmmTMB_17,
        newdata = newdata_15b, type = 2, model_type = "covariates"
      )
      ## ----end
      glmmTMB_17_mse
    }),
    tar_target(mse_3_mod_glmmTMB_17b_, {
      mse_glmmTMB <- mse_glmmTMB_
      mod_glmmTMB_17 <- mod_glmmTMB_17b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      ## ---- glmmTMB_17b_mse_3
      glmmTMB_17_mse <- mse_glmmTMB(mod_glmmTMB_17,
        newdata = newdata_14b, type = 3, model_type = "covariates"
      )
      ## ----end
      glmmTMB_17_mse
    }),
    
    ## brms -----------------------------------------------------------
    tar_target(mod_brms_17_, {
      benthos_fixed_locs_obs_17 <- incomplete_spatial_data_prep_17_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- brms_pre_17
      benthos_fixed_locs_obs_17 |>
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
      ## ---- brms_17
      mod_brms_17 <- brm(mod_form,
        data = benthos_fixed_locs_obs_17,
        iter = 5000,
        warmup = 1000,
        chains = 3,
        cores = 3,
        prior = priors,
        thin =  5,
        control = list(adapt_delta = 0.99),
        backend = "cmdstanr"
      )
      saveRDS(mod_brms_17,
        file = paste0(data_path, "synthetic/mod_brms_17.rds")
      ) 
      ## ----end
      mod_brms_17
    }),
    tar_target(pred_1_mod_brms_17_, {
      pred_brms <- pred_brms_
      mod_brms_17 <- mod_brms_17_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_17 <- incomplete_spatial_newdata_17_
      mod_simple_17 <- mod_simple_17_ #sampled_simple_raw_means_
      ## ---- brms_17_pred_1
      newdata <- newdata_17
      true_sum <- mod_simple_17
      brms_17_pred <- pred_brms(mod_brms_17,
        type = 1,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_17"
      )
      ## ----end
      brms_17_pred
    }),
    tar_target(pred_2_mod_brms_17_, {
      pred_brms <- pred_brms_
      mod_brms_17 <- mod_brms_17_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15 <- incomplete_spatial_newdata_15_
      benthos_reefs_temporal_summary_15 <- read_sampled_reefs_data_15_plot_2_
      ## ---- brms_17_pred_2
      newdata <- newdata_15
      true_sum <- benthos_reefs_temporal_summary_15
      brms_17_pred <- pred_brms(mod_brms_17,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_17"
      )
      ## ----end
      brms_17_pred
    }),
    tar_target(pred_3_mod_brms_17_, {
      pred_brms <- pred_brms_
      mod_brms_17 <- mod_brms_17_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14 <- incomplete_spatial_newdata_14_
      benthos_reefs_temporal_summary_14 <- read_sampled_reefs_data_14_plot_2_
      ## ---- brms_17_pred_3
      newdata <- newdata_14
      true_sum <- benthos_reefs_temporal_summary_14
      brms_17_pred <- pred_brms(mod_brms_17,
        type = 3,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_17"
      )
      ## ----end
      brms_17_pred
    }),
    tar_target(brms_trace_17_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_17 <- mod_brms_17_
      ## ---- brms_trace_17
      mod_brms_17 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_17.rds")
      )
      vars <- mod_brms_17 |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_17$fit |> stan_trace(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_trace_17.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ac_17_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_17 <- mod_brms_17_
      ## ---- brms_ac_17
      mod_brms_17 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_17.rds")
      )
      vars <- mod_brms_17 |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_17$fit |> stan_ac(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ac_17.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_rhat_17_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_17 <- mod_brms_17_
      ## ---- brms_rhat_17
      mod_brms_17 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_17.rds")
      )
      g <- mod_brms_17$fit |> stan_rhat()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_rhat_17.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ess_17_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_17 <- mod_brms_17_
      ## ---- brms_ess_17
      mod_brms_17 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_17.rds")
      )
      g <- mod_brms_17$fit |> stan_ess()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ess_17.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ppc_17_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_17 <- mod_brms_17_
      ## ---- brms_ppc_17
      mod_brms_17 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_17.rds")
      )
      g <- mod_brms_17 |> pp_check( type='dens_overlay', ndraws=140)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ppc_17.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(mse_1_mod_brms_17_, {
      mse_brms <- mse_brms_
      mod_brms_17 <- mod_brms_17_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_17b <- incomplete_spatial_newdata_17b_
      ## ---- brms_17_mse_1
      brms_17_mse <- mse_brms(mod_brms_17,
        newdata = newdata_17b, type = 1, model_type = ""
      )
      ## ----end
      brms_17_mse
    }),
    tar_target(mse_2_mod_brms_17_, {
      mse_brms <- mse_brms_
      mod_brms_17 <- mod_brms_17_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      ## ---- brms_17_mse_2
      brms_17_mse <- mse_brms(mod_brms_17,
        newdata = newdata_15b, type = 2, model_type = ""
      )
      ## ----end
      brms_17_mse
    }),
    tar_target(mse_3_mod_brms_17_, {
      mse_brms <- mse_brms_
      mod_brms_17 <- mod_brms_17_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      ## ---- brms_17_mse_3
      brms_17_mse <- mse_brms(mod_brms_17,
        newdata = newdata_14b, type = 3, model_type = ""
      )
      ## ----end
      brms_17_mse
    }),
    
    ## brms + covariates ----------------------------------------------
    tar_target(mod_brms_17b_, {
      benthos_fixed_locs_obs_17 <- incomplete_spatial_data_prep_17_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- brms_pre_17b
      benthos_fixed_locs_obs_17 |>
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
      mod_form <- bf(cover ~ fYear + scale(Longitude) + scale(Latitude) +
                       scale(CYC) + scale(DHW) + scale(OTHER) +
                       (1 | Site) + (1 | Transect),
        family = "Beta"
      )
      ## ----end
      ## ---- brms_17
      mod_brms_17b <- brm(mod_form,
        data = benthos_fixed_locs_obs_17,
        iter = 5000,
        warmup = 1000,
        chains = 3,
        cores = 3,
        prior = priors,
        thin =  5,
        control = list(adapt_delta = 0.99),
        backend = "cmdstanr"
      )
      saveRDS(mod_brms_17b,
        file = paste0(data_path, "synthetic/mod_brms_17b.rds")
      ) 
      ## ----end
      mod_brms_17b
    }),
    tar_target(pred_1_mod_brms_17b_, {
      pred_brms <- pred_brms_
      mod_brms_17b <- mod_brms_17b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_17b <- incomplete_spatial_newdata_17b_
      mod_simple_17 <- mod_simple_17_ #sampled_simple_raw_means_
      ## ---- brms_17b_pred_1
      true_sum <- mod_simple_17
      newdata <- newdata_17b
      brms_17b_pred <- pred_brms(mod_brms_17b,
        type = 1,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_17"
      )
      ## ----end
      brms_17b_pred
    }),
    tar_target(pred_2_mod_brms_17b_, {
      pred_brms <- pred_brms_
      mod_brms_17 <- mod_brms_17b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      benthos_reefs_temporal_summary_15 <- read_sampled_reefs_data_15_plot_2_
      ## ---- brms_17b_pred_2
      newdata <- newdata_15b
      true_sum <- benthos_reefs_temporal_summary_15
      brms_17b_pred <- pred_brms(mod_brms_17,
        type = 2,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_17"
      )
      ## ----end
      brms_17b_pred
    }),
    tar_target(pred_3_mod_brms_17b_, {
      pred_brms <- pred_brms_
      mod_brms_17 <- mod_brms_17b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      benthos_reefs_temporal_summary_14 <- read_sampled_reefs_data_14_plot_2_
      ## ---- brms_17b_pred_3
      newdata <- newdata_14b
      true_sum <- benthos_reefs_temporal_summary_14
      brms_17b_pred <- pred_brms(mod_brms_17,
        type = 3,
        model_type = "covariates",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_brms_17"
      )
      ## ----end
      brms_17b_pred
    }),
    tar_target(brms_trace_17b_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_17b <- mod_brms_17b_
      ## ---- brms_trace_17
      mod_brms_17b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_17b.rds")
      )
      vars <- mod_brms_17b |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_17b$fit |> stan_trace(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_trace_17b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ac_17b_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_17b <- mod_brms_17b_
      ## ---- brms_ac_17b
      mod_brms_17b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_17b.rds")
      )
      vars <- mod_brms_17b |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_17b$fit |> stan_ac(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ac_17b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_rhat_17b_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_17b <- mod_brms_17b_
      ## ---- brms_rhat_17b
      mod_brms_17b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_17b.rds")
      )
      g <- mod_brms_17b$fit |> stan_rhat()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_rhat_17b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ess_17b_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_17b <- mod_brms_17b_
      ## ---- brms_ess_17b
      mod_brms_17b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_17b.rds")
      )
      g <- mod_brms_17b$fit |> stan_ess()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ess_17b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ppc_17b_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_brms_17b <- mod_brms_17b_
      ## ---- brms_ppc_17b
      mod_brms_17b <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_17b.rds")
      )
      g <- mod_brms_17b |> pp_check( type='dens_overlay', ndraws=100)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ppc_17b.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(mse_1_mod_brms_17b_, {
      mse_brms <- mse_brms_
      mod_brms_17 <- mod_brms_17b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_17b <- incomplete_spatial_newdata_17b_
      ## ---- brms_17b_mse_1
      brms_17_mse <- mse_brms(mod_brms_17,
        newdata = newdata_17b, type = 1, model_type = "covariates"
      )
      ## ----end
      brms_17_mse
    }),
    tar_target(mse_2_mod_brms_17b_, {
      mse_brms <- mse_brms_
      mod_brms_17 <- mod_brms_17b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      ## ---- brms_17b_mse_2
      brms_17_mse <- mse_brms(mod_brms_17,
        newdata = newdata_15b, type = 2, model_type = "covariates"
      )
      ## ----end
      brms_17_mse
    }),
    tar_target(mse_3_mod_brms_17b_, {
      mse_brms <- mse_brms_
      mod_brms_17 <- mod_brms_17b_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      ## ---- brms_17b_mse_3
      brms_17_mse <- mse_brms(mod_brms_17,
        newdata = newdata_14b, type = 3, model_type = "covariates"
      )
      ## ----end
      brms_17_mse
    }),
    
    ## stan -----------------------------------------------------------
    tar_target(mod_stan_17_, {
      source("model_functions.R")
      benthos_fixed_locs_obs_17 <- incomplete_spatial_data_prep_17_
      data_path <- incomplete_spatial_global_parameters_$data_path
      site_extra_functions_
      ## ---- stan_pre_17
      benthos_fixed_locs_obs_17 <-
        benthos_fixed_locs_obs_17 |>
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
      saveRDS(benthos_fixed_locs_obs_17,
        file = paste0(data_path, "synthetic/saveRDS(benthos_fixed_locs_obs_17_forstan.rds")
      ) 
      stan_data <- prepare_data_for_stan(benthos_fixed_locs_obs_17, yrs = NULL)
      model_stan <- cmdstanr::cmdstan_model(stan_file = "model1.stan")
      ## ----end
      ## ---- stan_17
      mod_stan_17 <- model_stan$sample(
        data = stan_data,
        seed = 173,
        iter_sampling = 5000,
        iter_warmup = 1000,
        thin = 5,
        chains = 3,
        parallel_chains = 3,
        adapt_delta = 0.99,
        output_dir = paste0(data_path, "synthetic/"),
      )
      saveRDS(mod_stan_17,
        file = paste0(data_path, "synthetic/mod_stan_17.rds")
      ) 
      ## ----end
      mod_stan_17
    }),
    tar_target(pred_1_mod_stan_17_, {
      pred_stan <- pred_stan_
      mod_stan_17 <- mod_stan_17_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_17 <- incomplete_spatial_newdata_17_
      mod_simple_17 <- mod_simple_17_ #sampled_simple_raw_means_
      ## ---- stan_17_pred_1
      newdata <- newdata_17
      true_sum <- mod_simple_17
      stan_17_pred <- pred_stan(mod_stan_17,
        type = 1,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_stan_17"
      )
      ## ----end
      stan_17_pred
    }),
    tar_target(pred_2_mod_stan_17_, {
      pred_stan <- pred_stan_
      mod_stan_17 <- mod_stan_17_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15 <- incomplete_spatial_newdata_15_
      benthos_reefs_temporal_summary_15 <- read_sampled_reefs_data_15_plot_2_
      ## ---- stan_17_pred_2
      newdata <- newdata_15
      true_sum <- benthos_reefs_temporal_summary_15
      stan_17_pred <- pred_stan(mod_stan_17,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_stan_17"
      )
      ## ----end
      stan_17_pred
    }),
    tar_target(pred_3_mod_stan_17_, {
      pred_stan <- pred_stan_
      mod_stan_17 <- mod_stan_17_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14 <- incomplete_spatial_newdata_14_
      benthos_reefs_temporal_summary_14 <- read_sampled_reefs_data_14_plot_2_
      ## ---- stan_17_pred_3
      newdata <- newdata_14
      true_sum <- benthos_reefs_temporal_summary_14
      stan_17_pred <- pred_stan(mod_stan_17,
        type = 3,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_stan_17"
      )
      ## ----end
      stan_17_pred
    }),
    tar_target(stan_trace_17_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_stan_17 <- mod_stan_17_
      ## ---- stan_trace_17
      mod_stan_17 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_17.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_17$draws(variables = c("beta", "phi", "sd_1", "sd_2", "sd_3")) |>
        mcmc_trace() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_trace_17.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ac_17_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_stan_17 <- mod_stan_17_
      ## ---- stan_ac_17
      mod_stan_17 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_17.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_17$draws(variables = c("beta", "phi", "sd_1", "sd_2", "sd_3")) |>
        mcmc_acf() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ac_17.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_rhat_17_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_stan_17 <- mod_stan_17_
      ## ---- stan_rhat_17
      mod_stan_17 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_17.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_17 |> bayesplot::rhat() |> 
        mcmc_rhat_hist() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_rhat_17.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ess_17_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_stan_17 <- mod_stan_17_
      ## ---- stan_ess_17
      mod_stan_17 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_17.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_17 |> bayesplot::neff_ratio() |> 
        mcmc_neff_hist() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ess_17.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ppc_17_, {
      benthos_fixed_locs_obs_17 <- incomplete_spatial_data_prep_17_
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- incomplete_spatial_global_parameters_$fig_path
      mod_stan_17 <- mod_stan_17_
      ## ---- stan_ppc_17
      mod_stan_17 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_17.rds")
      )
      g <- 
        bayesplot::pp_check(
          benthos_fixed_locs_obs_17$cover,
          mod_stan_17$draws("ypred", format = "matrix")[1:100, ],
          ppc_dens_overlay
        ) +
        theme_classic()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ppc_17.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(mse_1_mod_stan_17_, {
      mse_stan <- mse_stan_
      mod_stan_17 <- mod_stan_17_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_17b <- incomplete_spatial_newdata_17b_
      ## ---- stan_17_mse_1
      stan_17_mse <- mse_stan(mod_stan_17,
        newdata = newdata_17b, type = 1, model_type = ""
      )
      ## ----end
      stan_17_mse
    }),
    tar_target(mse_2_mod_stan_17_, {
      mse_stan <- mse_stan_
      mod_stan_17 <- mod_stan_17_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      ## ---- stan_17_mse_2
      stan_17_mse <- mse_stan(mod_stan_17,
        newdata = newdata_15b, type = 2, model_type = ""
      )
      ## ----end
      stan_17_mse
    }),
    tar_target(mse_3_mod_stan_17_, {
      mse_stan <- mse_stan_
      mod_stan_17 <- mod_stan_17_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      ## ---- stan_17_mse_3
      stan_17_mse <- mse_stan(mod_stan_17,
        newdata = newdata_14b, type = 3, model_type = ""
      )
      ## ----end
      stan_17_mse
    }),
    
    ## gbm ------------------------------------------------------------
    tar_target(mod_gbm_17_, {
      benthos_fixed_locs_obs_17 <- incomplete_spatial_data_prep_17_
      data_path <- incomplete_spatial_global_parameters_$data_path
      ## ---- gbm_17
      mod_gbm_17 <- gbm(cover ~ fYear,
        data =  benthos_fixed_locs_obs_17,
        distribution = "gaussian",
        n.trees = 10000,
        interaction.depth = 5,
        shrinkage = 0.001,
        bag.fraction = 0.5,
        cv.folds = 5, n.cores = 1,
        verbose = TRUE
      )
      saveRDS(mod_gbm_17,
        file = paste0(data_path, "synthetic/mod_gbm_17.rds")
      ) 
      ## ----end
      ## ---- gbm_post_17
      n.trees <- gbm.perf(mod_gbm_17, method = "cv")
      ## ----end
      list(mod_gbm_17 = mod_gbm_17, n.trees = n.trees)
    }),
    tar_target(pred_1_mod_gbm_17_, {
      pred_gbm <- pred_gbm_
      mod_gbm_17 <- mod_gbm_17_$mod_gbm_17
      n.trees <- mod_gbm_17_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_17 <- incomplete_spatial_newdata_17_
      mod_simple_17 <- mod_simple_17_ #sampled_simple_raw_means_
      ## ---- gbm_17_pred_1
      newdata <- newdata_17
      true_sum <- mod_simple_17
      gbm_17_pred <- pred_gbm(mod_gbm_17,
        n.trees = n.trees,
        type = 1,
        model_type = "",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_17"
      )
      ## ----end
      gbm_17_pred
    }),
    tar_target(pred_2_mod_gbm_17_, {
      pred_gbm <- pred_gbm_
      mod_gbm_17 <- mod_gbm_17_$mod_gbm_17
      n.trees <- mod_gbm_17_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15 <- incomplete_spatial_newdata_15_
      benthos_reefs_temporal_summary_15 <- read_sampled_reefs_data_15_plot_2_
      ## ---- gbm_17_pred_2
      newdata <- newdata_15
      true_sum <- benthos_reefs_temporal_summary_15
      gbm_17_pred <- pred_gbm(mod_gbm_17,
        n.trees = n.trees,
        type = 2,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_17"
      )
      ## ----end
      gbm_17_pred
    }),
    tar_target(pred_3_mod_gbm_17_, {
      pred_gbm <- pred_gbm_
      mod_gbm_17 <- mod_gbm_17_$mod_gbm_17
      n.trees <- mod_gbm_17_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14 <- incomplete_spatial_newdata_14_
      benthos_reefs_temporal_summary_14 <- read_sampled_reefs_data_14_plot_2_
      ## ---- gbm_17_pred_3
      newdata <- newdata_14
      true_sum <- benthos_reefs_temporal_summary_14
      gbm_17_pred <- pred_gbm(mod_gbm_17,
        n.trees = n.trees,
        type = 3,
        model_type = "",
        newdata = newdata,
        true_sum = true_sum,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_17"
      )
      ## ----end
      gbm_17_pred
    }),
    tar_target(mse_1_mod_gbm_17_, {
      mse_gbm <- mse_gbm_
      mod_gbm_17 <- mod_gbm_17_$mod_gbm_17
      n.trees <- mod_gbm_17_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_17b <- incomplete_spatial_newdata_17b_
      ## ---- gbm_17_mse_1
      gbm_17_mse <- mse_gbm(mod_gbm_17,
        n.trees = n.trees,
        newdata = newdata_17b, type = 1, model_type = ""
      )
      ## ----end
      gbm_17_mse
    }),
    tar_target(mse_2_mod_gbm_17_, {
      mse_gbm <- mse_gbm_
      mod_gbm_17 <- mod_gbm_17_$mod_gbm_17
      n.trees <- mod_gbm_17_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      ## ---- gbm_17_mse_2
      gbm_17_mse <- mse_gbm(mod_gbm_17,
        n.trees = n.trees,
        newdata = newdata_15b, type = 2, model_type = ""
      )
      ## ----end
      gbm_17_mse
    }),
    tar_target(mse_3_mod_gbm_17_, {
      mse_gbm <- mse_gbm_
      mod_gbm_17 <- mod_gbm_17_$mod_gbm_17
      n.trees <- mod_gbm_17_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      ## ---- gbm_17_mse_3
      gbm_17_mse <- mse_gbm(mod_gbm_17,
        n.trees = n.trees,
        newdata = newdata_14b, type = 3, model_type = ""
      )
      ## ----end
      gbm_17_mse
    }),

    ## gbm + covartiates ----------------------------------------------
    tar_target(mod_gbm_17b_, {
      benthos_fixed_locs_obs_17 <- incomplete_spatial_data_prep_17_
      data_path <- incomplete_spatial_global_parameters_$data_path
      ## ---- gbm_17b
      mod_gbm_17b <- gbm(cover ~ fYear + Latitude + Longitude + CYC + DHW + OTHER,
        data =  benthos_fixed_locs_obs_17,
        distribution = "gaussian",
        var.monotone = c(0, 0, 0, -1, -1, -1),
        n.trees = 10000,
        interaction.depth = 5,
        shrinkage = 0.001,
        bag.fraction = 0.5,
        cv.folds = 5, n.cores = 1,
        verbose = TRUE
      )
      saveRDS(mod_gbm_17b,
        file = paste0(data_path, "synthetic/mod_gbm_17b.rds")
      )
      ## ----end
      ## ---- gbm_post_17b
      n.trees <- gbm.perf(mod_gbm_17b, method = "cv")
      ## ----end
      list(mod_gbm_17b = mod_gbm_17b, n.trees = n.trees)
    }),
    tar_target(pred_1_mod_gbm_17b_, {
      pred_gbm <- pred_gbm_
      mod_gbm_17b <- mod_gbm_17b_$mod_gbm_17
      n.trees <- mod_gbm_17b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_17b <- incomplete_spatial_newdata_17b_
      mod_simple_17 <- mod_simple_17_ #sampled_simple_raw_means_
      ## ---- gbm_17b_pred_1
      true_sum <- mod_simple_17
      newdata <- newdata_17b
      gbm_17_pred <- pred_gbm(mod_gbm_17b,
        n.trees = n.trees,
        type = 1,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_17"
      )
      ## ----end
      gbm_17_pred
    }),
    tar_target(pred_2_mod_gbm_17b_, {
      pred_gbm <- pred_gbm_
      mod_gbm_17b <- mod_gbm_17b_$mod_gbm_17
      n.trees <- mod_gbm_17b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      benthos_reefs_temporal_summary_15 <- read_sampled_reefs_data_15_plot_2_
      ## ---- gbm_17b_pred_2
      newdata <- newdata_15b
      true_sum <- benthos_reefs_temporal_summary_15
      gbm_17_pred <- pred_gbm(mod_gbm_17b,
        n.trees = n.trees,
        type = 2,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_17"
      )
      ## ----end
      gbm_17_pred
    }),
    tar_target(pred_3_mod_gbm_17b_, {
      pred_gbm <- pred_gbm_
      mod_gbm_17b <- mod_gbm_17b_$mod_gbm_17
      n.trees <- mod_gbm_17b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      benthos_reefs_temporal_summary_14 <- read_sampled_reefs_data_14_plot_2_
      ## ---- gbm_17b_pred_3
      newdata <- newdata_14b
      true_sum <- benthos_reefs_temporal_summary_14
      gbm_17_pred <- pred_gbm(mod_gbm_17b,
        n.trees = n.trees,
        type = 3,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_gbm_17"
      )
      ## ----end
      gbm_17_pred
    }),
    tar_target(infl_gbm_17b_, {
      mod_gbm_17b <- mod_gbm_17b_$mod_gbm_17b
      n.trees <- mod_gbm_17b_$n.trees
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      ## ---- gbm_infl_17b
      ## gbm_3b_infl <- gbm::relative.influence(mod_gbm_3b, n.trees = n.trees, scale = TRUE, sort = TRUE)
      infl <- summary(mod_gbm_17b, n.trees =  n.trees, plot = FALSE)
      g <-
        infl |>
        as.data.frame() |>
        arrange(rel.inf) |>
        mutate(var = factor(var, levels = unique(var))) |>
        ggplot(aes(y = var, x = rel.inf)) +
        geom_bar(stat = "identity") +
        scale_y_discrete("") +
        scale_x_continuous("Relative influence") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_infl_mod_gbm_17b.png"
        ),
        g,
        width = 6, height = 4, dpi = 72
      )
      ## ----end
    }),
    tar_target(mse_1_mod_gbm_17b_, {
      mse_gbm <- mse_gbm_
      mod_gbm_17 <- mod_gbm_17b_$mod_gbm_17b
      n.trees <- mod_gbm_17b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_17b <- incomplete_spatial_newdata_17b_
      ## ---- gbm_17b_mse_1
      gbm_17_mse <- mse_gbm(mod_gbm_17,
        n.trees = n.trees,
        newdata = newdata_17b, type = 1, model_type = "covariates"
      )
      ## ----end
      gbm_17_mse
    }),
    tar_target(mse_2_mod_gbm_17b_, {
      mse_gbm <- mse_gbm_
      mod_gbm_17 <- mod_gbm_17b_$mod_gbm_17b
      n.trees <- mod_gbm_17b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      ## ---- gbm_17b_mse_2
      gbm_17_mse <- mse_gbm(mod_gbm_17,
        n.trees = n.trees,
        newdata = newdata_15b, type = 2, model_type = "covariates"
      )
      ## ----end
      gbm_17_mse
    }),
    tar_target(mse_3_mod_gbm_17b_, {
      mse_gbm <- mse_gbm_
      mod_gbm_17 <- mod_gbm_17b_$mod_gbm_17b
      n.trees <- mod_gbm_17b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      ## ---- gbm_17b_mse_3
      gbm_17_mse <- mse_gbm(mod_gbm_17,
        n.trees = n.trees,
        newdata = newdata_14b, type = 3, model_type = "covariates"
      )
      ## ----end
      gbm_17_mse
    }),

    ## dbarts --------------------------------------------------------
    tar_target(mod_dbarts_17_, {
      benthos_fixed_locs_obs_17 <- incomplete_spatial_data_prep_17_
      benthos_fixed_locs_obs_14 <- incomplete_spatial_data_prep_14_
      benthos_fixed_locs_obs_15 <- incomplete_spatial_data_prep_15_
      data_path <- incomplete_spatial_global_parameters_$data_path
      newdata_17 <- incomplete_spatial_newdata_17_
      newdata_17b <- incomplete_spatial_newdata_17b_
      newdata_17_a <- incomplete_spatial_data_prep_17_#incomplete_spatial_newdata_17b_
      newdata_14 <- incomplete_spatial_newdata_14_
      newdata_14_a <- incomplete_spatial_newdata_14b_
      newdata_15 <- incomplete_spatial_newdata_15_
      newdata_15_a <- incomplete_spatial_newdata_15b_
      ## ---- dbarts_17
      print(head(benthos_fixed_locs_obs_17))
      mod_dbarts_17 <- bart2(log(cover) ~ fYear,
        data =   benthos_fixed_locs_obs_17,
        keepTrees = TRUE
      )
      saveRDS(mod_dbarts_17,
        file = paste0(data_path, "synthetic/mod_dbarts_17.rds")
      ) 
      ## ----end
      ## Unfortunately, the next part must be in the same tar_target
      ## due to the way dbarts stores pointers - they cannot be stored
      ## ---- dbarts_pred_17
      preds <- predict(mod_dbarts_17, newdata_17, type = "ev") |>
        exp()
      saveRDS(preds,
        file = paste0(data_path, "synthetic/mod_dbarts_17_preds.rds")
      ) 
      preds_sum <- preds |>
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_17_preds_sum.rds")
      ) 

      newdata <- newdata_17_a |>
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
      preds_a <- predict(mod_dbarts_17, newdata, type = "ev") |>
        exp() 
      saveRDS(preds_a,
        file = paste0(data_path, "synthetic/mod_dbarts_17_a_preds.rds")
      ) 

      newdata_17b <- newdata_17b |>
        mutate(fYear = factor(Year))
      preds_17b <- predict(mod_dbarts_17, newdata_17b, type = "ev") |>
        exp() 
      saveRDS(preds_17b,
        file = paste0(data_path, "synthetic/mod_dbarts_17b_preds.rds")
      ) 


      
      preds_14 <- predict(mod_dbarts_17, newdata_14, type = "ev") |>
        exp()
      saveRDS(preds_14,
        file = paste0(data_path, "synthetic/mod_dbarts_14_preds.rds")
      ) 
      preds_14_sum <- preds_14 |>
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_14_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_14_preds_sum.rds")
      ) 
      newdata_14_a <- newdata_14_a |>
        ungroup() |>
        mutate(fYear = factor(Year))
      preds_14_a <- predict(mod_dbarts_17, newdata_14_a, type = "ev") |>
        exp()
      saveRDS(preds_14_a,
        file = paste0(data_path, "synthetic/mod_dbarts_14_a_preds.rds")
      ) 

      preds_15 <- predict(mod_dbarts_17, newdata_15, type = "ev") |>
        exp()
      saveRDS(preds_15,
        file = paste0(data_path, "synthetic/mod_dbarts_15_preds.rds")
      ) 
      preds_15_sum <- preds_15 |> 
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_15_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_15_preds_sum.rds")
      ) 
      newdata_15_a <- newdata_15_a |>
        ungroup() |>
        mutate(fYear = factor(Year))
      preds_15_a <- predict(mod_dbarts_17, newdata_15_a, type = "ev") |>
        exp()
      saveRDS(preds_15_a,
        file = paste0(data_path, "synthetic/mod_dbarts_15_a_preds.rds")
      ) 
      ## ----end
      list(
        mod_dbarts_17 = mod_dbarts_17,
        preds = preds,
        preds_sum = preds_sum,
        preds_a = preds_a,
        preds_14 = preds_14,
        preds_14_sum = preds_14_sum,
        preds_14_a = preds_14_a,
        preds_15 = preds_15,
        preds_15_sum = preds_15_sum,
        preds_15_a = preds_15_a
      )
    }),
    tar_target(pred_1_mod_dbarts_17_, {
      pred_dbarts <- pred_dbarts_
      preds <- mod_dbarts_17_$preds_sum
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_17 <- incomplete_spatial_newdata_17_
      mod_simple_17 <- mod_simple_17_ #sampled_simple_raw_means_
      ## ---- dbarts_17_pred_1
      newdata <- newdata_17
      true_sum <- mod_simple_17
      dbarts_17_pred <- pred_dbarts(
        preds,
        type = 1,
        model_type = "",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_17"
      )
      ## ----end
      dbarts_17_pred
    }),
    tar_target(pred_2_mod_dbarts_17_, {
      pred_dbarts <- pred_dbarts_
      preds <- mod_dbarts_17_$preds_15_sum
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15 <- incomplete_spatial_newdata_15_
      benthos_reefs_temporal_summary_15 <- read_sampled_reefs_data_15_plot_2_
      ## ---- dbarts_17_pred_2
      newdata <- newdata_15
      true_sum <- benthos_reefs_temporal_summary_15
      dbarts_17_pred <- pred_dbarts(
        preds,
        type = 2,
        model_type = "",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_17"
      )
      ## ----end
      dbarts_17_pred
    }),
    tar_target(pred_3_mod_dbarts_17_, {
      pred_dbarts <- pred_dbarts_
      preds <- mod_dbarts_17_$preds_14_sum
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14 <- incomplete_spatial_newdata_14_
      benthos_reefs_temporal_summary_14 <- read_sampled_reefs_data_14_plot_2_
      ## ---- dbarts_17_pred_3
      newdata <- newdata_14
      true_sum <- benthos_reefs_temporal_summary_14
      dbarts_17_pred <- pred_dbarts(
        preds,
        type = 3,
        model_type = "",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_17"
      )
      ## ----end
      dbarts_17_pred
    }),
    tar_target(mse_1_mod_dbarts_17_, {
      mse_dbarts <- mse_dbarts_
      preds <- mod_dbarts_17_$preds_a
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_17b <- incomplete_spatial_newdata_17b_
      ## ---- dbarts_17_mse_1
      newdata_17b <- newdata_17b |>
        mutate(fYear = factor(Year)) |> 
        ungroup()
      dbarts_17_mse <- mse_dbarts(preds,
        newdata = newdata_17b, type = 1, model_type = ""
      )
      ## ----end
      dbarts_17_mse
    }),
    tar_target(mse_2_mod_dbarts_17_, {
      mse_dbarts <- mse_dbarts_
      preds <- mod_dbarts_17_$preds_15_a
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_ #incomplete_spatial_newdata_17b_
      ## ---- dbarts_17_mse_2
      dbarts_17_mse <- mse_dbarts(preds,
        newdata = newdata_15b, type = 2, model_type = ""
      )
      ## ----end
      dbarts_17_mse
    }),
    tar_target(mse_3_mod_dbarts_17_, {
      mse_dbarts <- mse_dbarts_
      preds <- mod_dbarts_17_$preds_14_a
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_ #incomplete_spatial_newdata_17b_
      ## ---- dbarts_17_mse_3
      ## newdata <- newdata |>
      ##   ungroup() |>
      ##   mutate(fYear = factor(Year))
      dbarts_17_mse <- mse_dbarts(preds,
        newdata = newdata_14b, type = 3, model_type = ""
      )
      ## ----end
      dbarts_17_mse
    }),
    
    ## dbarts + covariates --------------------------------------------
    tar_target(mod_dbarts_17b_, {
      benthos_fixed_locs_obs_17 <- incomplete_spatial_data_prep_17_
      benthos_fixed_locs_obs_14 <- incomplete_spatial_data_prep_14_
      benthos_fixed_locs_obs_15 <- incomplete_spatial_data_prep_15_
      data_path <- incomplete_spatial_global_parameters_$data_path
      newdata_17b <- incomplete_spatial_newdata_17b_
      ## newdata_14 <- incomplete_spatial_newdata_14_
      newdata_14b <- incomplete_spatial_newdata_14b_
      newdata_15b <- incomplete_spatial_newdata_15b_
      ## ---- dbarts_17b
      print(head(benthos_fixed_locs_obs_17))
      mod_dbarts_17b <- bart2(log(cover) ~ fYear + Latitude + Longitude + CYC + DHW + OTHER,
        data =   benthos_fixed_locs_obs_17,
        keepTrees = TRUE
      )
      saveRDS(mod_dbarts_17b,
        file = paste0(data_path, "synthetic/mod_dbarts_17b.rds")
      ) 
      ## ----end
      ## Unfortunately, the next part must be in the same tar_target
      ## due to the way dbarts stores pointers - they cannot be stored
      ## ---- dbarts_pred_17b
      newdata_17b <- newdata_17b |>
        mutate(fYear = factor(Year))
      preds_17b <- predict(mod_dbarts_17b, newdata_17b, type = "ev") |>
        exp() 
      saveRDS(preds_17b,
        file = paste0(data_path, "synthetic/mod_dbarts_17b_preds_17b.rds")
      ) 
      preds_17b_sum <- preds_17b |> summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_17b_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_17b_preds_17b_sum.rds")
      )
      
      ## newdata_14 <- newdata_14 |>
      ##   mutate(fYear = factor(Year))
      ## preds_14 <- predict(mod_dbarts_17b, newdata_14, type = "ev") |>
      ##   exp() 
      ## saveRDS(preds_14,
      ##   file = paste0(data_path, "synthetic/mod_dbarts_14b_preds.rds")
      ## ) 
      ## preds_14_sum <- preds_14 |> 
      ##   summarise_draws(median, HDInterval::hdi)
      ## saveRDS(preds_14_sum,
      ##   file = paste0(data_path, "synthetic/mod_dbarts_14b_preds_sum.rds")
      ## ) 

      newdata_14b <- newdata_14b |>
        mutate(fYear = factor(Year))
      preds_14b <- predict(mod_dbarts_17b, newdata_14b, type = "ev") |>
        exp() 
      saveRDS(preds_14b,
        file = paste0(data_path, "synthetic/mod_dbarts_17b_preds_14b.rds")
      ) 
      preds_14b_sum <- preds_14b |> 
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_14b_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_17b_preds_14b_sum.rds")
      ) 

      newdata_15b <- newdata_15b |>
        mutate(fYear = factor(Year))
      preds_15b <- predict(mod_dbarts_17b, newdata_15b, type = "ev") |>
        exp()
      saveRDS(preds_15b,
        file = paste0(data_path, "synthetic/mod_dbarts_17b_preds_15b.rds")
      ) 
      preds_15b_sum <- preds_15b |>
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds_15b_sum,
        file = paste0(data_path, "synthetic/mod_dbarts_17b_preds_15b_sum.rds")
      ) 
      ## ----end
      list(
        mod_dbarts_17 = mod_dbarts_17b,
        preds_17b = preds_17b,
        preds_17b_sum = preds_17b_sum,
        preds_14b = preds_14b,
        preds_14b_sum = preds_14b_sum,
        preds_15b = preds_15b,
        preds_15b_sum = preds_15b_sum
      )
    }),
    tar_target(pred_1_mod_dbarts_17b_, {
      pred_dbarts <- pred_dbarts_
      preds_17b <- mod_dbarts_17b_$preds_17b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_17b <- incomplete_spatial_newdata_17b_
      mod_simple_17 <- mod_simple_17_ #sampled_simple_raw_means_
      ## ---- dbarts_17b_pred_1
      true_sum <- mod_simple_17
      newdata <- newdata_17b
      dbarts_17_pred <- pred_dbarts(
        preds_17b,
        type = 1,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_17"
      )
      ## ----end
      dbarts_17_pred
    }),
    tar_target(pred_2_mod_dbarts_17b_, {
      pred_dbarts <- pred_dbarts_
      preds_15b <- mod_dbarts_17b_$preds_15b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      benthos_reefs_temporal_summary_15 <- read_sampled_reefs_data_15_plot_2_
      ## ---- dbarts_17b_pred_2
      newdata <- newdata_15b
      true_sum <- benthos_reefs_temporal_summary_15
      dbarts_17_pred <- pred_dbarts(
        preds_15b, 
        type = 2, 
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_17"
      )
      ## ----end
      dbarts_17_pred
    }),
    tar_target(pred_3_mod_dbarts_17b_, {
      pred_dbarts <- pred_dbarts_
      preds_14b <- mod_dbarts_17b_$preds_14b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      benthos_reefs_temporal_summary_14 <- read_sampled_reefs_data_14_plot_2_
      ## ---- dbarts_17b_pred_3
      newdata <- newdata_14b
      true_sum <- benthos_reefs_temporal_summary_14
      dbarts_17_pred <- pred_dbarts(
        preds_14b,
        type = 3,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_dbarts_17"
      )
      ## ----end
      dbarts_17_pred
    }),
    tar_target(mse_1_mod_dbarts_17b_, {
      mse_dbarts <- mse_dbarts_
      preds_17b <- mod_dbarts_17b_$preds_17b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata <- incomplete_spatial_data_prep_17_ #incomplete_spatial_newdata_17b_ #incomplete_spatial_newdata_17b_
      newdata_17b <- incomplete_spatial_newdata_17b_
      ## ---- dbarts_17b_mse_1
      dbarts_17_mse <- mse_dbarts(preds_17b,
        newdata = newdata_17b, type = 1, model_type = "covariates"
      )
      ## ----end
      dbarts_17_mse
    }),
    tar_target(mse_2_mod_dbarts_17b_, {
      mse_dbarts <- mse_dbarts_
      preds_15b <- mod_dbarts_17b_$preds_15b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      ## ---- dbarts_17b_mse_2
      dbarts_17_mse <- mse_dbarts(preds_15b,
        newdata = newdata_15b, type = 2, model_type = "covariates"
      )
      ## ----end
      dbarts_17_mse
    }),
    tar_target(mse_3_mod_dbarts_17b_, {
      mse_dbarts <- mse_dbarts_
      preds_14b <- mod_dbarts_17b_$preds_14b
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      ## ---- dbarts_17b_mse_3
      dbarts_17_mse <- mse_dbarts(preds_14b,
        newdata = newdata_14b, type = 3, model_type = "covariates"
      )
      ## ----end
      dbarts_17_mse
    }),
    
    ## xgboost + covariates -------------------------------------------
    tar_target(mod_xgboost_17b_prep_, {
      benthos_fixed_locs_obs_17 <- incomplete_spatial_data_prep_17_
      data_path <- incomplete_spatial_global_parameters_$data_path
      ## ---- xgboost_17_prep
      data_train_17b <- benthos_fixed_locs_obs_17 |>
        dplyr::select(cover, Year, Latitude, Longitude, CYC, DHW, OTHER) 
      ## ----end
      data_train_17b
    }),
    tar_target(mod_xgboost_17b_tune_, {
      data_train_17b <- mod_xgboost_17b_prep_
      data_path <- incomplete_spatial_global_parameters_$data_path
      ## ---- xgboost_17_tune
      ## Define the recipe
      tune_recipe <- recipe(cover ~ ., data = data_train_17b) |> 
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
        resamples = vfold_cv(data_train_17b, v = 5),
        grid = tune_grid_values)
      ## Get best set of parameters
      model_hyperparams <-
        select_best(tuned_results, metric = "rmse") |> 
        select(-".config") |> 
        as_tibble() |>
        mutate(nb_training = nrow(data_train_17b),
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
    tar_target(mod_xgboost_17b_fit_, {
      data_train_17b <- mod_xgboost_17b_prep_
      data_path <- incomplete_spatial_global_parameters_$data_path
      tune_recipe = mod_xgboost_17b_tune_$tune_recipe
      tune_model = mod_xgboost_17b_tune_$tune_model
      tune_workflow = mod_xgboost_17b_tune_$tune_workflow
      tune_grid_values = mod_xgboost_17b_tune_$tune_grid_values
      tuned_results = mod_xgboost_17b_tune_$tuned_results
      model_hyperparams = mod_xgboost_17b_tune_$model_hyperparams
      ## ---- xgboost_17_fit
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
      final_fitted_17b <- tune_workflow |>
        fit(data_train_17b)
      ## ----end
      final_fitted_17b
    }),
    tar_target(pred_1_mod_xgboost_17b_, {
      pred_xgboost <- pred_xgboost_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_17b <- incomplete_spatial_newdata_17b_
      final_fitted_17b <- mod_xgboost_17b_fit_
      mod_simple_17 <- mod_simple_17_ #sampled_simple_raw_means_
      ## ---- xgboost_17b_pred_1
      true_sum <- mod_simple_17
      newdata <- newdata_17b
      xgboost_17b_pred <- pred_xgboost(
        final_fitted_17b,
        type = 1,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_xgboost_17"
      )
      ## ----end
      xgboost_17b_pred
    }),
    tar_target(mse_1_mod_xgboost_17b_, {
      mse_xgboost <- mse_xgboost_
      final_fitted_17b <- mod_xgboost_17b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_17b <- incomplete_spatial_newdata_17b_
      ## ---- xgboost_17b_mse_1
      xgboost_17_mse <- mse_xgboost(final_fitted_17b,
        newdata = newdata_17b, type = 1, model_type = "covariates"
      )
      ## ----end
      xgboost_17_mse
    }),
    tar_target(pred_2_mod_xgboost_17b_, {
      pred_xgboost <- pred_xgboost_
      final_fitted_17b <- mod_xgboost_17b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      benthos_reefs_temporal_summary_15 <- read_sampled_reefs_data_15_plot_2_
      ## ---- xgboost_17b_pred_2
      newdata <- newdata_15b
      true_sum <- benthos_reefs_temporal_summary_15
      xgboost_17_pred <- pred_xgboost(final_fitted_17b,
        type = 2,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_xgboost_17"
      )
      ## ----end
      xgboost_17_pred
    }),
    tar_target(mse_2_mod_xgboost_17b_, {
      mse_xgboost <- mse_xgboost_
      final_fitted_17b <- mod_xgboost_17b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_15b <- incomplete_spatial_newdata_15b_
      ## ---- xgboost_17b_mse_2
      xgboost_17_mse <- mse_xgboost(final_fitted_17b,
        newdata = newdata_15b, type = 2, model_type = "covariates"
      )
      ## ----end
      xgboost_17_mse
    }),
    tar_target(pred_3_mod_xgboost_17b_, {
      pred_xgboost <- pred_xgboost_
      final_fitted_17b <- mod_xgboost_17b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      benthos_reefs_temporal_summary_14 <- read_sampled_reefs_data_14_plot_2_
      ## ---- xgboost_17b_pred_3
      newdata <- newdata_14b
      true_sum <- benthos_reefs_temporal_summary_14
      xgboost_17_pred <- pred_xgboost(final_fitted_17b,
        type = 3,
        model_type = "covariates",
        true_sum = true_sum,
        newdata = newdata,
        data_path = data_path,
        fig_path = fig_path,
        filename = "mod_xgboost_17"
      )
      ## ----end
      xgboost_17_pred
    }),
    tar_target(mse_3_mod_xgboost_17b_, {
      mse_xgboost <- mse_xgboost_
      final_fitted_17b <- mod_xgboost_17b_fit_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      newdata_14b <- incomplete_spatial_newdata_14b_
      ## ---- xgboost_17b_mse_3
      xgboost_17_mse <- mse_xgboost(final_fitted_17b,
        newdata = newdata_14b, type = 3, model_type = "covariates"
      )
      ## ----end
      xgboost_17_mse
    }),
    
    ## Comparisons ----------------------------------------------------
    tar_target(mse_1_mod_pymc_barts_17b_file_,
      paste0(
        incomplete_spatial_global_parameters_$data_path,
        "modelled/pymc_bart_17b_mse_1.csv"
      ), format = "file"),
    tar_target(mse_1_mod_pymc_barts_17b_, {
      read_csv(file = mse_1_mod_pymc_barts_17b_file_) 
    }),
    tar_target(mse_2_mod_pymc_barts_17b_file_,
      paste0(
        incomplete_spatial_global_parameters_$data_path,
        "modelled/pymc_bart_17b_mse_2.csv"
      ), format = "file"),
    tar_target(mse_2_mod_pymc_barts_17b_, {
      read_csv(file = mse_2_mod_pymc_barts_17b_file_) 
    }),
    tar_target(mse_3_mod_pymc_barts_17b_file_,
      paste0(
        incomplete_spatial_global_parameters_$data_path,
        "modelled/pymc_bart_17b_mse_3.csv"
      ), format = "file"),
    tar_target(mse_3_mod_pymc_barts_17b_, {
      read_csv(file = mse_3_mod_pymc_barts_17b_file_) 
    }),
    tar_target(comparisons_17_, {
      data_path <- incomplete_spatial_global_parameters_$data_path
      mse_1_mod_glmmTMB_17 <- mse_1_mod_glmmTMB_17_
      mse_2_mod_glmmTMB_17 <- mse_2_mod_glmmTMB_17_
      mse_3_mod_glmmTMB_17 <- mse_3_mod_glmmTMB_17_
      mse_1_mod_glmmTMB_17b <- mse_1_mod_glmmTMB_17b_
      mse_2_mod_glmmTMB_17b <- mse_2_mod_glmmTMB_17b_
      mse_3_mod_glmmTMB_17b <- mse_3_mod_glmmTMB_17b_

      mse_1_mod_brms_17 <- mse_1_mod_brms_17_
      mse_2_mod_brms_17 <- mse_2_mod_brms_17_
      mse_3_mod_brms_17 <- mse_3_mod_brms_17_
      mse_1_mod_brms_17b <- mse_1_mod_brms_17b_
      mse_2_mod_brms_17b <- mse_2_mod_brms_17b_
      mse_3_mod_brms_17b <- mse_3_mod_brms_17b_

      mse_1_mod_stan_17 <- mse_1_mod_stan_17_
      mse_2_mod_stan_17 <- mse_2_mod_stan_17_
      mse_3_mod_stan_17 <- mse_3_mod_stan_17_

      mse_1_mod_gbm_17 <- mse_1_mod_gbm_17_
      mse_2_mod_gbm_17 <- mse_2_mod_gbm_17_
      mse_3_mod_gbm_17 <- mse_3_mod_gbm_17_
      mse_1_mod_gbm_17b <- mse_1_mod_gbm_17b_
      mse_2_mod_gbm_17b <- mse_2_mod_gbm_17b_
      mse_3_mod_gbm_17b <- mse_3_mod_gbm_17b_

      mse_1_mod_dbarts_17 <- mse_1_mod_dbarts_17_
      mse_2_mod_dbarts_17 <- mse_2_mod_dbarts_17_
      mse_3_mod_dbarts_17 <- mse_3_mod_dbarts_17_
      mse_1_mod_dbarts_17b <- mse_1_mod_dbarts_17b_
      mse_2_mod_dbarts_17b <- mse_2_mod_dbarts_17b_
      mse_3_mod_dbarts_17b <- mse_3_mod_dbarts_17b_

      mse_1_mod_xgboost_17b <- mse_1_mod_xgboost_17b_
      mse_2_mod_xgboost_17b <- mse_2_mod_xgboost_17b_
      mse_3_mod_xgboost_17b <- mse_3_mod_xgboost_17b_

      mse_1_mod_pymc_barts_17b <- mse_1_mod_pymc_barts_17b_
      mse_2_mod_pymc_barts_17b <- mse_2_mod_pymc_barts_17b_
      mse_3_mod_pymc_barts_17b <- mse_3_mod_pymc_barts_17b_
      ## ---- comparisons_17
      mse_1_mod_pymc_barts_17b <- read_csv(
        file = paste0(data_path, "modelled/pymc_bart_17b_mse_1.csv")
      )   
      mse_2_mod_pymc_barts_17b <- read_csv(
        file = paste0(data_path, "modelled/pymc_bart_17b_mse_2.csv")
      )   
      mse_3_mod_pymc_barts_17b <- read_csv(
        file = paste0(data_path, "modelled/pymc_bart_17b_mse_3.csv")
      )   
      comparisons_17 <- bind_rows(
        mse_1_mod_glmmTMB_17,
        mse_2_mod_glmmTMB_17,
        mse_3_mod_glmmTMB_17,
        mse_1_mod_glmmTMB_17b,
        mse_2_mod_glmmTMB_17b,
        mse_3_mod_glmmTMB_17b,

        mse_1_mod_brms_17,
        mse_2_mod_brms_17,
        mse_3_mod_brms_17,
        mse_1_mod_brms_17b,
        mse_2_mod_brms_17b,
        mse_3_mod_brms_17b,
        
        mse_1_mod_stan_17,
        mse_2_mod_stan_17,
        mse_3_mod_stan_17,
        
        mse_1_mod_gbm_17,
        mse_2_mod_gbm_17,
        mse_3_mod_gbm_17,
        mse_1_mod_gbm_17b,
        mse_2_mod_gbm_17b,
        mse_3_mod_gbm_17b,

        mse_1_mod_dbarts_17,
        mse_2_mod_dbarts_17,
        mse_3_mod_dbarts_17,
        mse_1_mod_dbarts_17b,
        mse_2_mod_dbarts_17b,
        mse_3_mod_dbarts_17b,

        mse_1_mod_xgboost_17b,
        mse_2_mod_xgboost_17b,
        mse_3_mod_xgboost_17b,

        mse_1_mod_pymc_barts_17b,
        mse_2_mod_pymc_barts_17b,
        mse_3_mod_pymc_barts_17b
        )
      saveRDS(comparisons_17,
        file = paste0(data_path, "synthetic/comparisons_17.rds")
      ) 
       
      comps_17 <- 
        comparisons_17 |>
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
          type == 2 ~ "Predicting Western reefs",
          type == 3 ~ "Predicting Eastern reefs"
        )) 
      saveRDS(comps_17,
        file = paste0(data_path, "synthetic/comps_17.rds")
      )  
      ## ----end
      comps_17
    }),
    tar_target(comparisons_17_plots_, {
      comps_17 <- comparisons_17_
      data_path <- incomplete_spatial_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- comparisons_17_tab
      g <-
        comps_17 |> 
        filter(stat == "mean", metric == "mse") |>
        ggplot(aes(x = value, y = model)) +
        geom_segment(aes(xend = 0, yend = model), color = "black") +
        geom_point() +
        facet_grid(model_type~type, scales = "free_y") +
        scale_y_discrete("") +
        scale_x_continuous("MSE") +
        theme_bw() +
        ggtitle("Mean Square Error of models built on Western reefs data (25 reefs)")
 
      ggsave(
        filename = paste0(
          fig_path, "mse_17_1.png"
        ),
        g,
        width = 10, height = 6, dpi = 100
      )
      g <- 
        comps_17 |> 
        filter(stat == "mean", metric == "acc") |>
        ggplot(aes(x = value, y = model)) +
        geom_segment(aes(xend = 0, yend = model), color = "black") +
        geom_point() +
        facet_grid(model_type~type, scales = "free_y") +
        scale_y_discrete("") +
        scale_x_continuous("Mean inaccuracy (%)", labels = function(x) sprintf("%0.1f%%", x*100)) +
        theme_bw() +
        ggtitle("Mean inaccuracy (%) of models built on Western reefs data (2 reefs)")

      ggsave(
        filename = paste0(
          fig_path, "mse_17_2.png"
        ),
        g,
        width = 15, height = 6, dpi = 100
      )
      ## ----end
      comps_17
    })


  )

  return(targets)
}
