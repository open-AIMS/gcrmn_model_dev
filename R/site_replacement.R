 
source("model_functions.R")
## source("helper_functions.R")
site_replacement <- function() {
  targets <- list(
    tar_target(site_replacement_link,
      lnk <- synthetic_covariates_),
    # Target: Load raw data
    tar_target(
      site_replacement_libraries_,
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
      site_extra_functions_,
      {
        ## ---- site replacement functions
        source("model_functions.R")
        ## source("helper_functions.R")
        ## ----end
      }
    ),
    tar_target(
      site_replacement_global_parameters_,
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
    ## Import data ====================================================
    ## Full reef-level spatio-temporal data ---------------------------
    tar_target(read_all_reefs_data_, {
      data_path <- site_replacement_global_parameters_$data_path
      tmp <- synthetic_landscape_benthos_reefs_response_scale_
      ## ---- read all reefs data
      benthos_reefs_sf <- readRDS(file = paste0(
        data_path,
        "synthetic/benthos_reefs_sf.rds"
      ))
      ## ----end
      benthos_reefs_sf
    }),
    tar_target(read_all_temporal_summary_, {
      benthos_reefs_sf <- read_all_reefs_data_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- all reefs temporal summary
      benthos_reefs_temporal_summary <- benthos_reefs_sf |>
        st_drop_geometry() |>
        group_by(Year) |>
        summarise(Mean = mean(HCC),
          Median = median(HCC),
          SD = sd(HCC),
          Lower = quantile(HCC, 0.025),
          Upper = quantile(HCC, 0.975))
      benthos_reefs_temporal_summary |> 
        saveRDS(
          file = paste0(
            data_path,
            "synthetic/benthos_reefs_temporal_summary.rds"
          )
        )
      ## ----end
      benthos_reefs_temporal_summary 
    }),
    tar_target(read_all_temporal_summary_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_$benthos_reefs_temporal_summary
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- all reefs temporal summary plot
      benthos_reefs_temporal_summary <- 
        readRDS(
          file = paste0(
            data_path,
            "synthetic/benthos_reefs_temporal_summary.rds"
          )
        )
      g <- benthos_reefs_temporal_summary |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = Lower, ymax = Upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = Mean, colour = "mean")) +
        geom_line(aes(x = Year, y = Median, colour = "median")) +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_all_temporal_summary_plot.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
      ## ----end
    }),
    ## All sampled reefs ----------------------------------------------
    tar_target(read_sampled_reefs_data_, {
      data_path <- site_replacement_global_parameters_$data_path
      benthos_reefs_sf <- read_all_reefs_data_
      ## ---- read sampled reefs data
      benthos_fixed_locs_obs <- readRDS(file = paste0(
        data_path,
        "synthetic/benthos_fixed_locs_obs.rds"
      ))
      benthos_fixed_locs_obs <- benthos_fixed_locs_obs |>
        left_join(
          benthos_reefs_sf |>
            st_drop_geometry() |>
            dplyr::select(Year, Reef, CYC, DHW, OTHER) |>
            group_by(Year, Reef) |>
            summarise(across(c(CYC, DHW, OTHER), mean)),
          by = c("Year", "Reef")
        )
      ## ----end
      benthos_fixed_locs_obs
    }),
    tar_target(sampled_simple_raw_means_, {
      benthos_fixed_locs_obs <- read_sampled_reefs_data_
      ## ---- sampled simple raw means
      benthos_fixed_locs_obs <- benthos_fixed_locs_obs |>
        mutate(
          fYear = as.factor(Year),
          Reef = as.factor(Reef),
          Site = interaction(Reef, Site),
          Transect = interaction(Site, Transect),
          cover = HCC/100
        )
      print(benthos_fixed_locs_obs)
      ## Raw summary means
      all_sampled_sum <-
        benthos_fixed_locs_obs |>
        group_by(fYear) |>
        summarise(
          mean_mean_response = mean(cover),
          median_median_response = median(cover),
          mean_sd = sd(cover),
          mean_lower_lower = mean_mean_response - 1.96 * mean_sd,
          mean_upper_upper = mean_mean_response + 1.96 * mean_sd,
          median_lower_lower = quantile(cover, 0.025),
          median_upper_upper = quantile(cover, 0.975)
        ) |>
        dplyr::select(-mean_sd) |>
        pivot_longer(
          cols = c(
            mean_mean_response, median_median_response, mean_lower_lower,
            mean_upper_upper, median_lower_lower, median_upper_upper
          ),
          names_to = c("type", "variable", "stat"),
          names_sep = "_",
          values_to = "values"
        ) |>
        mutate(type = paste("all_sampled_", type)) |>
        dplyr::select(-variable) |>
        pivot_wider(
          names_from = stat,
          values_from = values
        ) |>
        mutate(Year = as.numeric(as.character(fYear)))
      print(all_sampled_sum)
      ## ----end
      all_sampled_sum
    }),
    tar_target(sampled_simple_raw_means_plot_, {
      all_sampled_sum <- sampled_simple_raw_means_
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- sampled simple raw means plot
      g <- all_sampled_sum |>
        ggplot(aes(x = Year, y = response, color = type, fill = type)) +
        geom_line() +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
        scale_y_continuous(
          name = "Hard Coral Cover (%)",
          limits = c(0, 1),
          labels = scales::percent_format(accuracy = 1)
        ) +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_full_simple_raw_means_plot.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
      ## ----end
    }),
    ## Replace a reef -------------------------------------------------
    tar_target(read_sampled_reefs_data_1_, {
      data_path <- site_replacement_global_parameters_$data_path
      benthos_reefs_sf <- read_all_reefs_data_
      ## ---- read sampled reefs data 1
      benthos_fixed_locs_obs_1 <- readRDS(
        file = paste0(
          data_path,
          "synthetic/benthos_fixed_locs_obs_1.rds"
        )
      )
      benthos_fixed_locs_obs_1 <- benthos_fixed_locs_obs_1 |>
        left_join(
          benthos_reefs_sf |>
            st_drop_geometry() |>
            dplyr::select(Year, Reef, CYC, DHW, OTHER) |>
            group_by(Year, Reef) |>
            summarise(across(c(CYC, DHW, OTHER), mean)),
          by = c("Year", "Reef")
        )
      ## ----end
      benthos_fixed_locs_obs_1 
    }),
    tar_target(sampled_reefs_data_1_plot_, {
      benthos_fixed_locs_obs_1 <- read_sampled_reefs_data_1_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- sampled reefs data 1 plot
      g <- benthos_fixed_locs_obs_1 |>
        ggplot() +
        geom_line(aes(
          y = HCC, x = Year, colour = Site,
          group = interaction(Site, Transect)
        )) +
        facet_wrap(~Reef) +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_sampled_reefs_1_plot.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
      ## ----end
    }),

    ## Replace a reef (V2) --------------------------------------------
    tar_target(read_sampled_reefs_data_2_, {
      data_path <- site_replacement_global_parameters_$data_path
      benthos_reefs_sf <- read_all_reefs_data_
      ## ---- read sampled reefs data 2
      benthos_fixed_locs_obs_2 <- readRDS(
        file = paste0(
          data_path,
          "synthetic/benthos_fixed_locs_obs_2.rds"
        )
      )
      benthos_fixed_locs_obs_2 <- benthos_fixed_locs_obs_2|>
        left_join(
          benthos_reefs_sf |>
            st_drop_geometry() |>
            dplyr::select(Year, Reef, CYC, DHW, OTHER) |>
            group_by(Year, Reef) |>
            summarise(across(c(CYC, DHW, OTHER), mean)),
          by = c("Year", "Reef")
        )
      ## ----end
      benthos_fixed_locs_obs_2 
    }),
    tar_target(sampled_reefs_data_2_plot_, {
      benthos_fixed_locs_obs_2 <- read_sampled_reefs_data_2_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- sampled reefs data 1 plot
      g <- benthos_fixed_locs_obs_2 |>
        ggplot() +
        geom_line(aes(
          y = HCC, x = Year, colour = Site,
          group = interaction(Site, Transect)
        )) +
        facet_wrap(~Reef) +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_sampled_reefs_2_plot.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
      ## ----end
    }),
    
    ## Data preparations ==============================================
    ## All sampled reefs ----------------------------------------------
    tar_target(site_replacements_data_prep_0_, {
      benthos_fixed_locs_obs <- read_sampled_reefs_data_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- sampled data prep 0
      benthos_fixed_locs_obs_0 <- benthos_fixed_locs_obs |>
        mutate(
          fYear = as.factor(Year),
          Reef = as.factor(Reef),
          Site = interaction(Reef, Site),
          Transect = interaction(Site, Transect),
          cover = HCC/100
        )
      ## ----end
      benthos_fixed_locs_obs_0 
    }),
    tar_target(site_replacements_newdata_0_, {
      benthos_fixed_locs_obs_0 <- site_replacements_data_prep_0_
      ## ---- newdata 0
      newdata_0 <-
        benthos_fixed_locs_obs_0 |>
        tidyr::expand(fYear)
      ## ----end
      newdata_0
    }),
    ## Replace a reef -------------------------------------------------
    tar_target(site_replacements_data_prep_1_, {
      benthos_fixed_locs_obs_1 <- read_sampled_reefs_data_1_
      ## ---- sampled data prep 1
      benthos_fixed_locs_obs_1 <- benthos_fixed_locs_obs_1 |>
        mutate(
          fYear = as.factor(Year),
          Reef = as.factor(Reef),
          Site = interaction(Reef, Site),
          Transect = interaction(Site, Transect),
          cover = HCC/100
        )
      ## ----end
      benthos_fixed_locs_obs_1 
    }),
    tar_target(site_replacements_newdata_1_, {
      benthos_fixed_locs_obs_1 <- site_replacements_data_prep_1_
      ## ---- newdata 1
      newdata_1 <-
        benthos_fixed_locs_obs_1 |>
        tidyr::expand(fYear)
      ## ----end
      newdata_1
    }),
    ## tar_target(site_replacements_newdata2_1_, {
    ##   benthos_fixed_locs_obs_1 <- site_replacements_data_prep_1_
    ##   ## ---- newdata 2
    ##   newdata2_1 <-
    ##     benthos_fixed_locs_obs_1 |>
    ##     tidyr::expand(fYear, Transect) |>
    ##     left_join(
    ##       benthos_fixed_locs_obs_1 |>
    ##         dplyr::select(Latitude, Longitude, Reef, Site, Transect) |>
    ##         distinct(),
    ##       by = c("Transect")
    ##     )
    ##   ## ----end
    ##   newdata2_1
    ## }),
    ## Replace a reef (V2) --------------------------------------------
    tar_target(site_replacements_data_prep_2_, {
      benthos_fixed_locs_obs_2 <- read_sampled_reefs_data_2_
      ## ---- sampled data prep 2
      benthos_fixed_locs_obs_2 <- benthos_fixed_locs_obs_2 |>
        mutate(
          fYear = as.factor(Year),
          Reef = as.factor(Reef),
          Site = interaction(Reef, Site),
          Transect = interaction(Site, Transect),
          cover = HCC/100
        )
      ## ----end
      benthos_fixed_locs_obs_2 
    }),
    tar_target(site_replacements_newdata_2_, {
      benthos_fixed_locs_obs_2 <- site_replacements_data_prep_2_
      ## ---- newdata 2
      newdata_2 <-
        benthos_fixed_locs_obs_2 |>
        tidyr::expand(fYear)
      ## ----end
      newdata_2
    }),
    ## tar_target(site_replacements_newdata2_2_, {
    ##   benthos_fixed_locs_obs_2 <- site_replacements_data_prep_2_
    ##   ## ---- newdata 2
    ##   newdata2_2 <-
    ##     benthos_fixed_locs_obs_2 |>
    ##     tidyr::expand(fYear, Transect) |>
    ##     left_join(
    ##       benthos_fixed_locs_obs_2 |>
    ##         dplyr::select(Latitude, Longitude, Reef, Site, Transect) |>
    ##         distinct(),
    ##       by = c("Transect")
    ##     )
    ##   ## ----end
    ##   newdata2_2
    ## }),
    ## Models =========================================================

    ## All sampled reefs ----------------------------------------------
    ## glmmTMB --------------------------------------------------------
    tar_target(mod_glmmTMB_0_, {
      benthos_fixed_locs_obs_0 <- site_replacements_data_prep_0_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_0
      mod_glmmTMB_0 <- glmmTMB(cover ~ fYear + (1 | Site) + (1 | Transect),
        data = benthos_fixed_locs_obs_0,
        family = "beta_family"
      )
      saveRDS(mod_glmmTMB_0,
        file = paste0(data_path, "synthetic/mod_glmmTMB_0.rds")
      ) 
      ## ----end
      mod_glmmTMB_0
    }),
    tar_target(dharma_mod_glmmTMB_0, {
      mod_glmmTMB_0 <- mod_glmmTMB_0_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_0_dharma
      glmmTMB_0_dharma <- mod_glmmTMB_0 |> 
        simulateResiduals(n = 1000)
      saveRDS(glmmTMB_0_dharma,
        file = paste0(data_path, "synthetic/glmmTMB_0_dharma.rds")
      ) 
      ## ----end
      glmmTMB_0_dharma
    }),
    tar_target(dharma_mod_glmmTMB_0_plot_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      glmmTMB_0_dharma <- dharma_mod_glmmTMB_0
      ## ---- glmmTMB_0_dharma plot
      glmmTMB_0_dharma <- readRDS(
        file = paste0(data_path, "synthetic/glmmTMB_0_dharma.rds")
      )
      g <- wrap_elements(~testUniformity(glmmTMB_0_dharma)) +
        wrap_elements(~plotResiduals(glmmTMB_0_dharma)) +
        wrap_elements(~testDispersion(glmmTMB_0_dharma))
      ggsave(
        filename = paste0(
          fig_path, "R_dharma_mod_glmmTMB_0.png"
        ),
        g,
        width = 10, height = 4, dpi = 72
      )
      ## ----end
    }),
    tar_target(emmeans_mod_glmmTMB_0, {
      mod_glmmTMB_0 <- mod_glmmTMB_0_
      newdata_0 <- site_replacements_newdata_0_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_0_emmeans
      glmmTMB_0_sum <- 
        mod_glmmTMB_0 |> emmeans(~fYear, at = newdata_0, type = "response") |>
        as.data.frame() |> 
        rename(median = response, lower = asymp.LCL, upper = asymp.UCL) |>
        mutate(Year = as.numeric(as.character(fYear))) |>
        mutate(type = "glmmTMB")
      saveRDS(glmmTMB_0_sum,
        file = paste0(data_path, "synthetic/glmmTMB_0_sum.rds")
      ) 
      ## ----end
      glmmTMB_0_sum 
    }),
    tar_target(emmeans_mod_glmmTMB_0_plot_, {
      glmmTMB_0_sum <- emmeans_mod_glmmTMB_0
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- glmmTMB_0_emmeans plot
      g <- glmmTMB_0_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "glmmTMB")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_glmmTMB_0.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
      ## ----end
    }),
    
    ## brms -----------------------------------------------------------
    tar_target(mod_brms_0_, {
      benthos_fixed_locs_obs_0 <- site_replacements_data_prep_0_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- brms_pre_0
      benthos_fixed_locs_obs_0 |>
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
      ## ---- brms_0
      mod_brms_0 <- brm(mod_form,
        data = benthos_fixed_locs_obs_0,
        iter = 5000,
        warmup = 1000,
        chains = 3,
        cores = 3,
        prior = priors,
        thin =  5,
        control = list(adapt_delta = 0.99),
        backend = "cmdstanr"
      )
      saveRDS(mod_brms_0,
        file = paste0(data_path, "synthetic/mod_brms_0.rds")
      ) 
      ## ----end
      mod_brms_0
    }),
    tar_target(emmeans_mod_brms_0_, {
      mod_brms_0 <- mod_brms_0_
      newdata_0 <- site_replacements_newdata_0_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- brms_0_emmeans
      brms_0_sum <- 
        mod_brms_0 |> emmeans(~fYear, at = newdata_0, type = "response") |>
        as.data.frame() |> 
        rename(median = response, lower = lower.HPD, upper = upper.HPD) |>
        mutate(Year = as.numeric(as.character(fYear))) |>
        mutate(type = "brms")
      saveRDS(brms_0_sum,
        file = paste0(data_path, "synthetic/brms_0_sum.rds")
      ) 
      ## ----end
      brms_0_sum 
    }),
    tar_target(emmeans_mod_brms_0_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      brms_0_sum <- emmeans_mod_brms_0_
      ## ---- brms_0_emmeans plot
      brms_0_sum <- readRDS(
        file = paste0(data_path, "synthetic/brms_0_sum.rds")
      )
      g <- 
        brms_0_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "brms")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_brms_0.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_trace_0_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_brms_0 <- mod_brms_0_
      ## ---- brms_trace_0
      mod_brms_0 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_0.rds")
      )
      vars <- mod_brms_0 |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_0$fit |> stan_trace(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_trace_0.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ac_0_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_brms_0 <- mod_brms_0_
      ## ---- brms_ac_0
      mod_brms_0 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_0.rds")
      )
      vars <- mod_brms_0 |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_0$fit |> stan_ac(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ac_0.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_rhat_0_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_brms_0 <- mod_brms_0_
      ## ---- brms_rhat_0
      mod_brms_0 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_0.rds")
      )
      g <- mod_brms_0$fit |> stan_rhat()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_rhat_0.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ess_0_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_brms_0 <- mod_brms_0_
      ## ---- brms_ess_0
      mod_brms_0 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_0.rds")
      )
      g <- mod_brms_0$fit |> stan_ess()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ess_0.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ppc_0_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_brms_0 <- mod_brms_0_
      ## ---- brms_ppc_0
      mod_brms_0 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_0.rds")
      )
      g <- mod_brms_0 |> pp_check( type='dens_overlay', ndraws=100)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ppc_0.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
      
    ## stan
    tar_target(mod_stan_0_, {
      source("model_functions.R")
      benthos_fixed_locs_obs_0 <- site_replacements_data_prep_0_
      data_path <- site_replacement_global_parameters_$data_path
      site_extra_functions_
      ## ---- stan_pre_0
      benthos_fixed_locs_obs_0 <-
        benthos_fixed_locs_obs_0 |>
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
      saveRDS(benthos_fixed_locs_obs_0,
        file = paste0(data_path, "synthetic/saveRDS(benthos_fixed_locs_obs_0_forstan.rds")
      ) 
      stan_data <- prepare_data_for_stan(benthos_fixed_locs_obs_0, yrs = NULL)
      model_stan <- cmdstanr::cmdstan_model(stan_file = "model1.stan")
      ## model_stan <- cmdstanr::cmdstan_model(stan_file = "mod1a.stan")
      ## model_stan <- cmdstanr::cmdstan_model(stan_file = "mod2a.stan")
      ## ----end
      ## ---- stan_0
      mod_stan_0 <- model_stan$sample(
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
      saveRDS(mod_stan_0,
        file = paste0(data_path, "synthetic/mod_stan_0.rds")
      ) 
      ## ----end
      mod_stan_0
    }),
    tar_target(pdp_mod_stan_0_, {
      mod_stan_0 <- mod_stan_0_
      newdata_0 <- site_replacements_newdata_0_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- stan_0_pdp
      stan_0_sum <-
        mod_stan_0$draws(variables = "cellmeans") |>
        posterior::as_draws_df() |>
        posterior::summarise_draws(
          median,
          HDInterval::hdi,
          ~ HDInterval::hdi(., credMass = c(0.9)),
          rhat,
          ess_bulk,
          ess_tail
        ) |>
        rename(lower_90 = V4, upper_90 = V5) |>
        bind_cols(newdata_0) |>
        mutate(Year = as.numeric(as.character(fYear)))
      saveRDS(stan_0_sum,
        file = paste0(data_path, "synthetic/stan_0_sum.rds")
      ) 
      ## ----end
      stan_0_sum 
    }),    
    tar_target(pdp_mod_stan_0_plot_, {
      stan_0_sum <- pdp_mod_stan_0_
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- stan_0_pdp plot
      stan_0_sum <- readRDS(
        file = paste0(data_path, "synthetic/stan_0_sum.rds")
      )
      g <- stan_0_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "stan")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_stan_0.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_trace_0_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_stan_0 <- mod_stan_0_
      ## ---- stan_trace_0
      mod_stan_0 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_0.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_0$draws(variables = c("beta", "phi", "sd_1", "sd_2", "sd_3")) |>
        mcmc_trace() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_trace_0.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ac_0_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_stan_0 <- mod_stan_0_
      ## ---- stan_ac_0
      mod_stan_0 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_0.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_0$draws(variables = c("beta", "phi", "sd_1", "sd_2", "sd_3")) |>
        mcmc_acf() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ac_0.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_rhat_0_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_stan_0 <- mod_stan_0_
      ## ---- stan_rhat_0
      mod_stan_0 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_0.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_0 |> bayesplot::rhat() |> 
        mcmc_rhat_hist() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_rhat_0.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ess_0_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_stan_0 <- mod_stan_0_
      ## ---- stan_ess_0
      mod_stan_0 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_0.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_0 |> bayesplot::neff_ratio() |> 
        mcmc_neff_hist() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ess_0.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ppc_0_, {
      benthos_fixed_locs_obs_0 <- site_replacements_data_prep_0_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_stan_0 <- mod_stan_0_
      ## ---- stan_ppc_0
      mod_stan_0 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_0.rds")
      )
      g <- 
        bayesplot::pp_check(
          benthos_fixed_locs_obs_0$cover,
          mod_stan_0$draws("ypred", format = "matrix")[1:100, ],
          ppc_dens_overlay
        ) +
        theme_classic()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ppc_0.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    
    ## gbm
    tar_target(mod_gbm_0_, {
      benthos_fixed_locs_obs_0 <- site_replacements_data_prep_0_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- gbm_0
      mod_gbm_0 <- gbm(cover ~ fYear,
        data =  benthos_fixed_locs_obs_0,
        distribution = "gaussian",
        n.trees = 10000,
        interaction.depth = 5,
        shrinkage = 0.001,
        bag.fraction = 0.5,
        cv.folds = 5,
        verbose = TRUE
      )
      saveRDS(mod_gbm_0,
        file = paste0(data_path, "synthetic/mod_gbm_0.rds")
      ) 
      ## ----end
      ## ---- gbm_post_0
      n.trees <- gbm.perf(mod_gbm_0, method = "cv")
      ## ----end
      list(mod_gbm_0 = mod_gbm_0, n.trees = n.trees)
    }),
    tar_target(pdp_gbm_0_, {
      mod_gbm_0 <- mod_gbm_0_$mod_gbm_0
      n.trees <- mod_gbm_0_$n.trees
      newdata_0 <- site_replacements_newdata_0_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- gbm_pdp_0
      gbm_0_sum <- newdata_0 |>
        mutate(median = predict(mod_gbm_0, newdata_0, n.trees = n.trees, type = "response")) |> 
        mutate(Year = as.numeric(as.character(fYear))) |>
        mutate(type = "gbm")
      saveRDS(gbm_0_sum,
        file = paste0(data_path, "synthetic/gbm_0_sum.rds")
      ) 
      ## ----end
      gbm_0_sum
    }),
    tar_target(pdp_gbm_0_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      gdm_0_sum <- pdp_gbm_0_
      ## ---- gbm_pdp_0
      gbm_0_sum <- readRDS(
        file = paste0(data_path, "synthetic/gbm_0_sum.rds")
      )
      g <-
        gbm_0_sum |>
        ggplot() +
        geom_line(aes(x = Year, y = median, color = "gbm")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_gbm_0.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
      ## ----end
    }),

    ## gbm + covartiates
    tar_target(mod_gbm_0b_, {
      benthos_fixed_locs_obs_0 <- site_replacements_data_prep_0_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- gbm_0b
      mod_gbm_0b <- gbm(cover ~ fYear + Latitude + Longitude + CYC + DHW + OTHER,
        data =  benthos_fixed_locs_obs_0,
        distribution = "gaussian",
        var.monotone = c(0, 0, 0, -1, -1, -1),
        n.trees = 10000,
        interaction.depth = 5,
        shrinkage = 0.001,
        bag.fraction = 0.5,
        cv.folds = 5,
        verbose = TRUE
      )
      saveRDS(mod_gbm_0b,
        file = paste0(data_path, "synthetic/mod_gbm_0b.rds")
      )
      ## ----end
      ## ---- gbm_post_0b
      n.trees <- gbm.perf(mod_gbm_0b, method = "cv")
      ## ----end
      list(mod_gbm_0b = mod_gbm_0b, n.trees = n.trees)
    }),
    tar_target(pdp_gbm_0b_, {
      benthos_fixed_locs_obs_0 <- site_replacements_data_prep_0_
      mod_gbm_0b <- mod_gbm_0b_$mod_gbm_0b
      n.trees <- mod_gbm_0b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- gbm_pdp_0b
      newdata_0b <- benthos_fixed_locs_obs_0
      gbm_0b_sum <- newdata_0b |>
        mutate(median = predict(mod_gbm_0b, newdata_0b, n.trees = n.trees, type = "response")) |>
        mutate(Year = as.numeric(as.character(fYear))) |>
        ## group_by(Year, Site, Transect) |>
        ## summarise(median = median(median)) |>
        ## group_by(Year, Site, .add = FALSE) |>
        ## summarise(median = median(median)) |>
        group_by(Year, .add = FALSE) |> 
        summarise(median = median(median)) |> 
        mutate(type = "gbm")
      saveRDS(gbm_0b_sum,
        file = paste0(data_path, "synthetic/gbm_0b_sum.rds")
      ) 
      ## ----end
      gbm_0b_sum
    }),
    tar_target(pdp_gbm_0b_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      gbm_0b_sum <- pdp_gbm_0b_
      ## ---- gbm_pdp_0b plot
      gbm_0b_sum <- readRDS(
        file = paste0(data_path, "synthetic/gbm_0b_sum.rds")
      )
      g <-
        gbm_0b_sum |>
        ggplot() +
        geom_line(aes(x = Year, y = median, color = "gbm")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_gbm_0b.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(infl_gbm_0b_, {
      mod_gbm_0b <- mod_gbm_0b_$mod_gbm_0b
      n.trees <- mod_gbm_0b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- gbm_infl_0b
      ## gbm_0b_infl <- gbm::relative.influence(mod_gbm_0b, n.trees = n.trees, scale = TRUE, sort = TRUE)
      infl <- summary(mod_gbm_0b, n.trees =  n.trees, plot = FALSE)
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
          fig_path, "R_infl_mod_gbm_0b.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(pdp_gbm_0c_, {
      benthos_reefs_sf <- read_all_reefs_data_
      mod_gbm_0b <- mod_gbm_0b_$mod_gbm_0b
      n.trees <- mod_gbm_0b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- gbm_pdp_0c
      newdata_0c <- benthos_reefs_sf |>
        mutate(
          Latitude = st_coordinates(geometry)[,2],
          Longitude = st_coordinates(geometry)[,1],
          fYear = as.factor(Year),
          )
      gbm_0c_sum <- newdata_0c |>
        mutate(median = predict(mod_gbm_0b, newdata_0c, n.trees = n.trees, type = "response")) |>
        mutate(Year = as.numeric(as.character(fYear))) |>
        ## group_by(Year, Site, Transect) |>
        ## summarise(median = median(median)) |>
        ## group_by(Year, Site, .add = FALSE) |>
        ## summarise(median = median(median)) |>
        group_by(Year, .add = FALSE) |> 
        summarise(median = median(median)) |> 
        mutate(type = "gbm")
      saveRDS(gbm_0c_sum,
        file = paste0(data_path, "synthetic/gbm_0c_sum.rds")
      ) 
      ## ----end
      gbm_0c_sum
    }),
    tar_target(pdp_gbm_0c_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      gbm_0c_sum <- pdp_gbm_0c_
      ## ---- gbm_pdp_0c plot
      gbm_0c_sum <- readRDS(
        file = paste0(data_path, "synthetic/gbm_0c_sum.rds")
      )
      g <-
        gbm_0c_sum |>
        ggplot() +
        geom_line(aes(x = Year, y = median, color = "gbm")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_gbm_0c.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(infl_gbm_0c_, {
      mod_gbm_0b <- mod_gbm_0b_$mod_gbm_0b
      n.trees <- mod_gbm_0b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- gbm_infl_0c
      ## gbm_0b_infl <- gbm::relative.influence(mod_gbm_0b, n.trees = n.trees, scale = TRUE, sort = TRUE)
      infl <- summary(mod_gbm_0b, n.trees =  n.trees, plot = FALSE)
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
          fig_path, "R_infl_mod_gbm_0c.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
      ## ----end
    }),

    ## dbarts --------------------------------------------------------
    tar_target(mod_dbarts_0_, {
      benthos_fixed_locs_obs_0 <- site_replacements_data_prep_0_
      data_path <- site_replacement_global_parameters_$data_path
      newdata_0 <- site_replacements_newdata_0_
      ## ---- dbarts_0
      print(head(benthos_fixed_locs_obs_0))
      mod_dbarts_0 <- bart2(log(cover) ~ fYear,
        data =   benthos_fixed_locs_obs_0,
        keepTrees = TRUE
      )
      saveRDS(mod_dbarts_0,
        file = paste0(data_path, "synthetic/mod_dbarts_0.rds")
      ) 
      ## ----end
      ## Unfortunately, the next part must be in the same tar_target
      ## due to the way dbarts stores pointers - they cannot be stored
      ## ---- dbarts_pred_0
      preds <- predict(mod_dbarts_0, newdata_0, type = "ev") |>
        exp() |> 
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds,
        file = paste0(data_path, "synthetic/mod_dbarts_0_preds.rds")
      ) 
      ## ----end
      list(
        mod_dbarts_0 = mod_dbarts_0,
        preds = preds
      )
    }),
    tar_target(dbarts_pdp_0, {
      data_path <- site_replacement_global_parameters_$data_path
      newdata_0 <- site_replacements_newdata_0_
      preds <- mod_dbarts_0_$preds
      ## ---- dbarts_pdp_0
      dbarts_0_sum <- newdata_0 |>
        bind_cols(preds) |> 
        mutate(Year = as.numeric(as.character(fYear))) |>
        mutate(type = "dbarts")
      saveRDS(dbarts_0_sum,
        file = paste0(data_path, "synthetic/dbarts_0_sum.rds")
      ) 
      # ----end
      dbarts_0_sum 
    }),
    tar_target(pdp_dbarts_0_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      dbarts_0_sum <- dbarts_pdp_0
      ## ---- dbarts_pdp_0 plot
      dbarts_0_sum <- readRDS(
        file = paste0(data_path, "synthetic/dbarts_0_sum.rds")
      )
      g <- 
        dbarts_0_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "dbarts")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_dbarts_0.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
      ## ----end
    }),
    
    ## dbarts + covariates ---------------------------------------------
    tar_target(mod_dbarts_0b_, {
      benthos_fixed_locs_obs_0 <- site_replacements_data_prep_0_
      benthos_reefs_sf <- read_all_reefs_data_
      data_path <- site_replacement_global_parameters_$data_path
      newdata_0 <- site_replacements_newdata_0_
      ## ---- dbarts_0b
      mod_dbarts_0b <- bart2(log(cover) ~ fYear + Latitude + Longitude + CYC + DHW + OTHER,
        data =   benthos_fixed_locs_obs_0,
        keepTrees = TRUE
      )
      saveRDS(mod_dbarts_0b,
        file = paste0(data_path, "synthetic/mod_dbarts_0b.rds")
      ) 
      ## ----end
      ## Unfortunately, the next part must be in the same tar_target
      ## due to the way dbarts stores pointers - they cannot be stored
      ## ---- dbarts_pred_0b
      newdata_0b <- benthos_fixed_locs_obs_0
      preds <- predict(mod_dbarts_0b, newdata_0b, type = "ev") |>
        exp()
      saveRDS(preds,
        file = paste0(data_path, "synthetic/mod_dbarts_0b_preds.rds")
      ) 
      ## ----end
      ## ---- dbarts_pred_0c
      newdata_0c <- benthos_reefs_sf |>
        mutate(
          Latitude = st_coordinates(geometry)[, 2],
          Longitude = st_coordinates(geometry)[, 1],
          fYear = as.factor(Year),
        ) |>
        st_drop_geometry() |>
        group_by(Year, fYear, Reef) |>
        summarise(across(c(Latitude, Longitude, CYC, DHW, OTHER), mean))
      preds_0c <- predict(mod_dbarts_0b, newdata_0c, type = "ev") |>
        exp()
      saveRDS(preds_0c,
        file = paste0(data_path, "synthetic/mod_dbarts_0c_preds.rds")
      ) 
      ## ----end
      list(
        mod_dbarts_0b = mod_dbarts_0b,
        preds = preds,
        preds_0c = preds_0c
      )
    }),
    tar_target(pdp_dbarts_0b_, {
      preds <- mod_dbarts_0b_$preds
      data_path <- site_replacement_global_parameters_$data_path
      benthos_fixed_locs_obs_0 <- site_replacements_data_prep_0_
      ## ---- dbarts_pdp_0
      dbarts_0b_sum <- preds |> 
        t() |>
        as_tibble(.name_repair = ~ paste0("V", seq_along(.))) |>
        bind_cols(benthos_fixed_locs_obs_0) |>
        pivot_longer(
          cols = matches("^V[0-9]*"),
          names_to = ".draw", values_to = "fit"
        ) |>
        group_by(Year, .draw) |>
        summarise(fit = mean(fit)) |>
        ## dplyr::select(Year, fit, .draw) |> 
        group_by(Year, .add = FALSE) |> 
        summarise_draws(median, HDInterval::hdi)
      saveRDS(dbarts_0b_sum,
        file = paste0(data_path, "synthetic/dbarts_0b_sum.rds")
      ) 
      ## ----end
      dbarts_0b_sum
    }),
    tar_target(pdp_dbarts_0b_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      dbarts_0b_sum <- pdp_dbarts_0b_
      ## ---- dbarts_pdp_0 plot
      dbarts_0b_sum <- readRDS(
        file = paste0(data_path, "synthetic/dbarts_0b_sum.rds")
      )
      g <- 
        dbarts_0b_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "dbarts")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_dbarts_0b.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(pdp_dbarts_0c_, {
      preds_0c <- mod_dbarts_0b_$preds_0c
      benthos_reefs_sf <- read_all_reefs_data_ |>
        mutate(
          Latitude = st_coordinates(geometry)[, 2],
          Longitude = st_coordinates(geometry)[, 1],
          fYear = as.factor(Year),
        ) |>
        st_drop_geometry() |>
        group_by(Year, fYear, Reef) |>
        summarise(across(c(Latitude, Longitude, CYC, DHW, OTHER), mean))
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- gbm_pdp_0c
      dbarts_0c_sum <- preds_0c |> 
        t() |>
        as_tibble(.name_repair = ~ paste0("V", seq_along(.))) |>
        bind_cols(benthos_reefs_sf) |>
        pivot_longer(
          cols = matches("^V[0-9]*"),
          names_to = ".draw", values_to = "fit"
        ) |>
        group_by(Year, .draw) |>
        summarise(fit = mean(fit)) |>
        ## dplyr::select(Year, fit, .draw) |> 
        group_by(Year, .add = FALSE) |> 
        summarise_draws(median, HDInterval::hdi)
      saveRDS(dbarts_0c_sum,
        file = paste0(data_path, "synthetic/dbarts_0c_sum.rds")
      ) 
      ## ----end
      dbarts_0c_sum
    }),
    tar_target(pdp_dbarts_0c_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      dbarts_0c_sum <- pdp_dbarts_0c_
      ## ---- dbarts_pdp_0 plot
      dbarts_0c_sum <- readRDS(
        file = paste0(data_path, "synthetic/dbarts_0c_sum.rds")
      )
      g <- 
        dbarts_0c_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "dbarts")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_dbarts_0c.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
      ## ----end
    }),
    
    ## Replace a reef (V1) --------------------------------------------
    ## simple ---------------------------------------------------------
    tar_target(mod_simple_1_, {
      benthos_fixed_locs_obs_1 <- site_replacements_data_prep_1_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- simple_1
      mod_simple_1 <- benthos_fixed_locs_obs_1 |>
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
      saveRDS(mod_simple_1,
        file = paste0(data_path, "synthetic/mod_simple_1.rds")
      ) 
      ## ----end
      mod_simple_1
    }),

    ## glmmTMB --------------------------------------------------------
    tar_target(mod_glmmTMB_1_, {
      benthos_fixed_locs_obs_1 <- site_replacements_data_prep_1_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_1
      mod_glmmTMB_1 <- glmmTMB(cover ~ fYear + (1 | Site) + (1 | Transect),
        data = benthos_fixed_locs_obs_1,
        family = "beta_family"
      )
      saveRDS(mod_glmmTMB_1,
        file = paste0(data_path, "synthetic/mod_glmmTMB_1.rds")
      ) 
      ## ----end
      mod_glmmTMB_1
    }),
    tar_target(dharma_mod_glmmTMB_1_, {
      mod_glmmTMB_1 <- mod_glmmTMB_1_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_1_dharma
      glmmTMB_1_dharma <- mod_glmmTMB_1 |> 
        simulateResiduals(n = 1000)
      saveRDS(glmmTMB_1_dharma,
        file = paste0(data_path, "synthetic/glmmTMB_1_dharma.rds")
      ) 
      ## ----end
      glmmTMB_1_dharma
    }),
    tar_target(dharma_mod_glmmTMB_1_plot_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      glmmTMB_1_dharma <- dharma_mod_glmmTMB_1_$glmmTMB_1_dharma
      ## ---- glmmTMB_1_dharma plot
      glmmTMB_1_dharma <- readRDS(
        file = paste0(data_path, "synthetic/glmmTMB_1_dharma.rds")
      )
      g <- wrap_elements(~testUniformity(glmmTMB_1_dharma)) +
        wrap_elements(~plotResiduals(glmmTMB_1_dharma)) +
        wrap_elements(~testDispersion(glmmTMB_1_dharma))
      ggsave(
        filename = paste0(
          fig_path, "R_dharma_mod_glmmTMB_1.png"
        ),
        g,
        width = 10, height = 4, dpi = 72
      )
      ## ----end
    }),
    tar_target(emmeans_mod_glmmTMB_1, {
      mod_glmmTMB_1 <- mod_glmmTMB_1_
      newdata_1 <- site_replacements_newdata_1_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_1_emmeans
      glmmTMB_1_sum <- 
        mod_glmmTMB_1 |> emmeans(~fYear, at = newdata_1, type = "response") |>
        as.data.frame() |> 
        rename(median = response, lower = asymp.LCL, upper = asymp.UCL) |>
        mutate(Year = as.numeric(as.character(fYear))) |>
        mutate(type = "glmmTMB")
      saveRDS(glmmTMB_1_sum,
        file = paste0(data_path, "synthetic/glmmTMB_1_sum.rds")
      ) 
      ## ----end
      glmmTMB_1_sum 
    }),
    tar_target(emmeans_mod_glmmTMB_1_plot_, {
      glmmTMB_1_sum <- emmeans_mod_glmmTMB_1
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_simple_1 <- mod_simple_1_
      ## ---- glmmTMB_1_emmeans plot
      g1 <- glmmTMB_1_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "glmmTMB")) +
        geom_line(data = mod_simple_1,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_1,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        ## geom_line(data = all_sampled_sum,
        ##   aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      g2 <- glmmTMB_1_sum |>
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
          fig_path, "R_pdp_mod_glmmTMB_1.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),

    ## brms -----------------------------------------------------------
    tar_target(mod_brms_1_, {
      benthos_fixed_locs_obs_1 <- site_replacements_data_prep_1_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- brms_pre_1
      benthos_fixed_locs_obs_1 |>
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
      ## ---- brms_1
      mod_brms_1 <- brm(mod_form,
        data = benthos_fixed_locs_obs_1,
        iter = 5000,
        warmup = 1000,
        chains = 3,
        cores = 3,
        prior = priors,
        thin =  5,
        control = list(adapt_delta = 0.99),
        backend = "cmdstanr"
      )
      saveRDS(mod_brms_1,
        file = paste0(data_path, "synthetic/mod_brms_1.rds")
      ) 
      ## ----end
      mod_brms_1
    }),
    tar_target(emmeans_mod_brms_1_, {
      mod_brms_1 <- mod_brms_1_
      newdata_1 <- site_replacements_newdata_1_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- brms_1_emmeans
      brms_1_sum <- 
        mod_brms_1 |> emmeans(~fYear, at = newdata_1, type = "response") |>
        as.data.frame() |> 
        rename(median = response, lower = lower.HPD, upper = upper.HPD) |>
        mutate(Year = as.numeric(as.character(fYear))) |>
        mutate(type = "brms")
      saveRDS(brms_1_sum,
        file = paste0(data_path, "synthetic/brms_1_sum.rds")
      ) 
      ## ----end
      brms_1_sum 
    }),
    tar_target(emmeans_mod_brms_1_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      brms_1_sum <- emmeans_mod_brms_1_
      mod_simple_1 <- mod_simple_1_
      ## ---- brms_1_emmeans plot
      brms_1_sum <- readRDS(
        file = paste0(data_path, "synthetic/brms_1_sum.rds")
      )
      g1 <- 
        brms_1_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "brms")) +
        geom_line(data = mod_simple_1,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_1,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <- 
        brms_1_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "brms")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_brms_1.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_trace_1_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_brms_1 <- mod_brms_1_
      ## ---- brms_trace_1
      mod_brms_1 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_1.rds")
      )
      vars <- mod_brms_1 |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_1$fit |> stan_trace(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_trace_1.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ac_1_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_brms_1 <- mod_brms_1_
      ## ---- brms_ac_1
      mod_brms_1 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_1.rds")
      )
      vars <- mod_brms_1 |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_1$fit |> stan_ac(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ac_1.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_rhat_1_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_brms_1 <- mod_brms_1_
      ## ---- brms_rhat_1
      mod_brms_1 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_1.rds")
      )
      g <- mod_brms_1$fit |> stan_rhat()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_rhat_1.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ess_1_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_brms_1 <- mod_brms_1_
      ## ---- brms_ess_1
      mod_brms_1 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_1.rds")
      )
      g <- mod_brms_1$fit |> stan_ess()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ess_1.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ppc_1_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_brms_1 <- mod_brms_1_
      ## ---- brms_ppc_1
      mod_brms_1 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_1.rds")
      )
      g <- mod_brms_1 |> pp_check( type='dens_overlay', ndraws=100)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ppc_1.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),

    ## stan -----------------------------------------------------------
    tar_target(mod_stan_1_, {
      source("model_functions.R")
      benthos_fixed_locs_obs_1 <- site_replacements_data_prep_1_
      data_path <- site_replacement_global_parameters_$data_path
      site_extra_functions_
      ## ---- stan_pre_1
      benthos_fixed_locs_obs_1 <-
        benthos_fixed_locs_obs_1 |>
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
      saveRDS(benthos_fixed_locs_obs_1,
        file = paste0(data_path, "synthetic/saveRDS(benthos_fixed_locs_obs_1_forstan.rds")
      ) 
      stan_data <- prepare_data_for_stan(benthos_fixed_locs_obs_1, yrs = NULL)
      model_stan <- cmdstanr::cmdstan_model(stan_file = "model1.stan")
      ## model_stan <- cmdstanr::cmdstan_model(stan_file = "mod1a.stan")
      ## model_stan <- cmdstanr::cmdstan_model(stan_file = "mod2a.stan")
      ## ----end
      ## ---- stan_1
      mod_stan_1 <- model_stan$sample(
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
      saveRDS(mod_stan_1,
        file = paste0(data_path, "synthetic/mod_stan_1.rds")
      ) 
      ## ----end
      mod_stan_1
    }),
    tar_target(pdp_mod_stan_1_, {
      mod_stan_1 <- mod_stan_1_
      newdata_1 <- site_replacements_newdata_1_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- stan_1_pdp
      stan_1_sum <-
        mod_stan_1$draws(variables = "cellmeans") |>
        posterior::as_draws_df() |>
        posterior::summarise_draws(
          median,
          HDInterval::hdi,
          ~ HDInterval::hdi(., credMass = c(0.9)),
          rhat,
          ess_bulk,
          ess_tail
        ) |>
        rename(lower_90 = V4, upper_90 = V5) |>
        bind_cols(newdata_1) |>
        mutate(Year = as.numeric(as.character(fYear)))
      saveRDS(stan_1_sum,
        file = paste0(data_path, "synthetic/stan_1_sum.rds")
      ) 
      ## ----end
      stan_1_sum 
    }),    
    tar_target(pdp_mod_stan_1_plot_, {
      stan_1_sum <- pdp_mod_stan_1_
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_simple_1 <- mod_simple_1_
      ## ---- stan_1_pdp plot
      stan_1_sum <- readRDS(
        file = paste0(data_path, "synthetic/stan_1_sum.rds")
      )
      g1 <- stan_1_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "stan")) +
        geom_line(data = mod_simple_1,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_1,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <- stan_1_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "stan")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_stan_1.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_trace_1_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_stan_1 <- mod_stan_1_
      ## ---- stan_trace_1
      mod_stan_1 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_1.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_1$draws(variables = c("beta", "phi", "sd_1", "sd_2", "sd_3")) |>
        mcmc_trace() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_trace_1.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ac_1_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_stan_1 <- mod_stan_1_
      ## ---- stan_ac_1
      mod_stan_1 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_1.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_1$draws(variables = c("beta", "phi", "sd_1", "sd_2", "sd_3")) |>
        mcmc_acf() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ac_1.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_rhat_1_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_stan_1 <- mod_stan_1_
      ## ---- stan_rhat_1
      mod_stan_1 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_1.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_1 |> bayesplot::rhat() |> 
        mcmc_rhat_hist() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_rhat_1.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ess_1_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_stan_1 <- mod_stan_1_
      ## ---- stan_ess_1
      mod_stan_1 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_1.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_1 |> bayesplot::neff_ratio() |> 
        mcmc_neff_hist() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ess_1.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ppc_1_, {
      benthos_fixed_locs_obs_1 <- site_replacements_data_prep_1_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_stan_1 <- mod_stan_1_
      ## ---- stan_ppc_1
      mod_stan_1 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_1.rds")
      )
      g <- 
        bayesplot::pp_check(
          benthos_fixed_locs_obs_1$cover,
          mod_stan_1$draws("ypred", format = "matrix")[1:100, ],
          ppc_dens_overlay
        ) +
        theme_classic()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ppc_1.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),

    ## gbm ------------------------------------------------------------
    tar_target(mod_gbm_1_, {
      benthos_fixed_locs_obs_1 <- site_replacements_data_prep_1_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- gbm_1
      mod_gbm_1 <- gbm(cover ~ fYear,
        data =  benthos_fixed_locs_obs_1,
        distribution = "gaussian",
        n.trees = 10000,
        interaction.depth = 5,
        shrinkage = 0.001,
        bag.fraction = 0.5,
        cv.folds = 5,
        verbose = TRUE
      )
      saveRDS(mod_gbm_1,
        file = paste0(data_path, "synthetic/mod_gbm_1.rds")
      ) 
      ## ----end
      ## ---- gbm_post_1
      n.trees <- gbm.perf(mod_gbm_1, method = "cv")
      ## ----end
      list(mod_gbm_1 = mod_gbm_1, n.trees = n.trees)
    }),
    tar_target(pdp_gbm_1_, {
      mod_gbm_1 <- mod_gbm_1_$mod_gbm_1
      n.trees <- mod_gbm_1_$n.trees
      newdata_1 <- site_replacements_newdata_1_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- gbm_pdp_1
      gbm_1_sum <- newdata_1 |>
        mutate(median = predict(mod_gbm_1, newdata_1, n.trees = n.trees, type = "response")) |> 
        mutate(Year = as.numeric(as.character(fYear))) |>
        mutate(type = "gbm")
      saveRDS(gbm_1_sum,
        file = paste0(data_path, "synthetic/gbm_1_sum.rds")
      ) 
      ## ----end
      gbm_1_sum
    }),
    tar_target(pdp_gbm_1_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      gdm_1_sum <- pdp_gbm_1_
      mod_simple_1 <- mod_simple_1_
      ## ---- gbm_pdp_1
      gbm_1_sum <- readRDS(
        file = paste0(data_path, "synthetic/gbm_1_sum.rds")
      )
      g1 <-
        gbm_1_sum |>
        ggplot() +
        geom_line(aes(x = Year, y = median, color = "gbm")) +
        geom_line(data = mod_simple_1,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_1,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <-
        gbm_1_sum |>
        ggplot() +
        geom_line(aes(x = Year, y = median, color = "gbm")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_gbm_1.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    
    ## gbm + covartiates ----------------------------------------------
    tar_target(mod_gbm_1b_, {
      benthos_fixed_locs_obs_1 <- site_replacements_data_prep_1_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- gbm_1b
      mod_gbm_1b <- gbm(cover ~ fYear + Latitude + Longitude + CYC + DHW + OTHER,
        data =  benthos_fixed_locs_obs_1,
        distribution = "gaussian",
        var.monotone = c(0, 0, 0, -1, -1, -1),
        n.trees = 10000,
        interaction.depth = 5,
        shrinkage = 0.001,
        bag.fraction = 0.5,
        cv.folds = 5,
        verbose = TRUE
      )
      saveRDS(mod_gbm_1b,
        file = paste0(data_path, "synthetic/mod_gbm_1b.rds")
      )
      ## ----end
      ## ---- gbm_post_1b
      n.trees <- gbm.perf(mod_gbm_1b, method = "cv")
      ## ----end
      list(mod_gbm_1b = mod_gbm_1b, n.trees = n.trees)
    }),
    tar_target(pdp_gbm_1b_, {
      benthos_fixed_locs_obs_1 <- site_replacements_data_prep_1_
      mod_gbm_1b <- mod_gbm_1b_$mod_gbm_1b
      n.trees <- mod_gbm_1b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- gbm_pdp_1b
      newdata_1b <- benthos_fixed_locs_obs_1
      gbm_1b_sum <- newdata_1b |>
        mutate(median = predict(mod_gbm_1b, newdata_1b, n.trees = n.trees, type = "response")) |>
        mutate(Year = as.numeric(as.character(fYear))) |>
        ## group_by(Year, Site, Transect) |>
        ## summarise(median = median(median)) |>
        ## group_by(Year, Site, .add = FALSE) |>
        ## summarise(median = median(median)) |>
        group_by(Year, .add = FALSE) |> 
        summarise(median = median(median)) |> 
        mutate(type = "gbm")
      saveRDS(gbm_1b_sum,
        file = paste0(data_path, "synthetic/gbm_1b_sum.rds")
      ) 
      ## ----end
      gbm_1b_sum
    }),
    tar_target(pdp_gbm_1b_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      gbm_1b_sum <- pdp_gbm_1b_
      mod_simple_1 <- mod_simple_1_
      ## ---- gbm_pdp_1b plot
      gbm_1b_sum <- readRDS(
        file = paste0(data_path, "synthetic/gbm_1b_sum.rds")
      )
      g1 <-
        gbm_1b_sum |>
        ggplot() +
        geom_line(aes(x = Year, y = median, color = "gbm")) +
        geom_line(data = mod_simple_1,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_1,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <-
        gbm_1b_sum |>
        ggplot() +
        geom_line(aes(x = Year, y = median, color = "gbm")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_gbm_1b.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(infl_gbm_1b_, {
      mod_gbm_1b <- mod_gbm_1b_$mod_gbm_1b
      n.trees <- mod_gbm_1b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- gbm_infl_1b
      ## gbm_1b_infl <- gbm::relative.influence(mod_gbm_1b, n.trees = n.trees, scale = TRUE, sort = TRUE)
      infl <- summary(mod_gbm_1b, n.trees =  n.trees, plot = FALSE)
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
          fig_path, "R_infl_mod_gbm_1b.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(pdp_gbm_1c_, {
      benthos_reefs_sf <- read_all_reefs_data_
      mod_gbm_1b <- mod_gbm_1b_$mod_gbm_1b
      n.trees <- mod_gbm_1b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- gbm_pdp_1c
      newdata_1c <- benthos_reefs_sf |>
        mutate(
          Latitude = st_coordinates(geometry)[,2],
          Longitude = st_coordinates(geometry)[,1],
          fYear = as.factor(Year),
          )
      gbm_1c_sum <- newdata_1c |>
        mutate(median = predict(mod_gbm_1b, newdata_1c, n.trees = n.trees, type = "response")) |>
        mutate(Year = as.numeric(as.character(fYear))) |>
        ## group_by(Year, Site, Transect) |>
        ## summarise(median = median(median)) |>
        ## group_by(Year, Site, .add = FALSE) |>
        ## summarise(median = median(median)) |>
        group_by(Year, .add = FALSE) |> 
        summarise(median = median(median)) |> 
        mutate(type = "gbm")
      saveRDS(gbm_1c_sum,
        file = paste0(data_path, "synthetic/gbm_1c_sum.rds")
      ) 
      ## ----end
      gbm_1c_sum
    }),
    tar_target(pdp_gbm_1c_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      gbm_1c_sum <- pdp_gbm_1c_
      mod_simple_1 <- mod_simple_1_
      ## ---- gbm_pdp_1c plot
      gbm_1c_sum <- readRDS(
        file = paste0(data_path, "synthetic/gbm_1c_sum.rds")
      )
      g1 <-
        gbm_1c_sum |>
        ggplot() +
        geom_line(aes(x = Year, y = median, color = "gbm")) +
        geom_line(data = mod_simple_1,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_1,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <-
        gbm_1c_sum |>
        ggplot() +
        geom_line(aes(x = Year, y = median, color = "gbm")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_gbm_1c.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(infl_gbm_1c_, {
      mod_gbm_1b <- mod_gbm_1b_$mod_gbm_1b
      n.trees <- mod_gbm_1b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- gbm_infl_1c
      ## gbm_1b_infl <- gbm::relative.influence(mod_gbm_1b, n.trees = n.trees, scale = TRUE, sort = TRUE)
      infl <- summary(mod_gbm_1b, n.trees =  n.trees, plot = FALSE)
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
          fig_path, "R_infl_mod_gbm_1c.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
      ## ----end
    }),

    ## dbarts --------------------------------------------------------
    tar_target(mod_dbarts_1_, {
      benthos_fixed_locs_obs_1 <- site_replacements_data_prep_1_
      data_path <- site_replacement_global_parameters_$data_path
      newdata_1 <- site_replacements_newdata_1_
      ## ---- dbarts_1
      print(head(benthos_fixed_locs_obs_1))
      mod_dbarts_1 <- bart2(log(cover) ~ fYear,
        data =   benthos_fixed_locs_obs_1,
        keepTrees = TRUE
      )
      saveRDS(mod_dbarts_1,
        file = paste0(data_path, "synthetic/mod_dbarts_1.rds")
      ) 
      ## ----end
      ## Unfortunately, the next part must be in the same tar_target
      ## due to the way dbarts stores pointers - they cannot be stored
      ## ---- dbarts_pred_1
      preds <- predict(mod_dbarts_1, newdata_1, type = "ev") |>
        exp() |> 
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds,
        file = paste0(data_path, "synthetic/mod_dbarts_1_preds.rds")
      ) 
      ## ----end
      list(
        mod_dbarts_1 = mod_dbarts_1,
        preds = preds
      )
    }),
    tar_target(dbarts_pdp_1, {
      data_path <- site_replacement_global_parameters_$data_path
      newdata_1 <- site_replacements_newdata_1_
      preds <- mod_dbarts_1_$preds
      ## ---- dbarts_pdp_1
      dbarts_1_sum <- newdata_1 |>
        bind_cols(preds) |> 
        mutate(Year = as.numeric(as.character(fYear))) |>
        mutate(type = "dbarts")
      saveRDS(dbarts_1_sum,
        file = paste0(data_path, "synthetic/dbarts_1_sum.rds")
      ) 
      # ----end
      dbarts_1_sum 
    }),
    tar_target(pdp_dbarts_1_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      dbarts_1_sum <- dbarts_pdp_1
      mod_simple_1 <- mod_simple_1_
      ## ---- dbarts_pdp_1 plot
      dbarts_1_sum <- readRDS(
        file = paste0(data_path, "synthetic/dbarts_1_sum.rds")
      )
      g1 <- 
        dbarts_1_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "dbarts")) +
        geom_line(data = mod_simple_1,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_1,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <- 
        dbarts_1_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "dbarts")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_dbarts_1.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    
    ## dbarts + covariates ---------------------------------------------
    tar_target(mod_dbarts_1b_, {
      benthos_fixed_locs_obs_1 <- site_replacements_data_prep_1_
      benthos_reefs_sf <- read_all_reefs_data_
      data_path <- site_replacement_global_parameters_$data_path
      newdata_1 <- site_replacements_newdata_1_
      ## ---- dbarts_1b
      mod_dbarts_1b <- bart2(log(cover) ~ fYear + Latitude + Longitude + CYC + DHW + OTHER,
        data =   benthos_fixed_locs_obs_1,
        keepTrees = TRUE
      )
      saveRDS(mod_dbarts_1b,
        file = paste0(data_path, "synthetic/mod_dbarts_1b.rds")
      ) 
      ## ----end
      ## Unfortunately, the next part must be in the same tar_target
      ## due to the way dbarts stores pointers - they cannot be stored
      ## ---- dbarts_pred_1b
      newdata_1b <- benthos_fixed_locs_obs_1
      preds <- predict(mod_dbarts_1b, newdata_1b, type = "ev") |>
        exp()
      saveRDS(preds,
        file = paste0(data_path, "synthetic/mod_dbarts_1b_preds.rds")
      ) 
      ## ----end
      ## ---- dbarts_pred_1c
      newdata_1c <- benthos_reefs_sf |>
        mutate(
          Latitude = st_coordinates(geometry)[, 2],
          Longitude = st_coordinates(geometry)[, 1],
          fYear = as.factor(Year),
        ) |>
        st_drop_geometry() |>
        group_by(Year, fYear, Reef) |>
        summarise(across(c(Latitude, Longitude, CYC, DHW, OTHER), mean))
      preds_1c <- predict(mod_dbarts_1b, newdata_1c, type = "ev") |>
        exp()
      saveRDS(preds_1c,
        file = paste0(data_path, "synthetic/mod_dbarts_1c_preds.rds")
      ) 
      ## ----end
      list(
        mod_dbarts_1b = mod_dbarts_1b,
        preds = preds,
        preds_1c = preds_1c
      )
    }),
    tar_target(pdp_dbarts_1b_, {
      preds <- mod_dbarts_1b_$preds
      data_path <- site_replacement_global_parameters_$data_path
      benthos_fixed_locs_obs_1 <- site_replacements_data_prep_1_
      ## ---- dbarts_pdp_1
      dbarts_1b_sum <- preds |> 
        t() |>
        as_tibble(.name_repair = ~ paste0("V", seq_along(.))) |>
        bind_cols(benthos_fixed_locs_obs_1) |>
        pivot_longer(
          cols = matches("^V[0-9]*"),
          names_to = ".draw", values_to = "fit"
        ) |>
        group_by(Year, .draw) |>
        summarise(fit = mean(fit)) |>
        ## dplyr::select(Year, fit, .draw) |> 
        group_by(Year, .add = FALSE) |> 
        summarise_draws(median, HDInterval::hdi)
      saveRDS(dbarts_1b_sum,
        file = paste0(data_path, "synthetic/dbarts_1b_sum.rds")
      ) 
      ## ----end
      dbarts_1b_sum
    }),
    tar_target(pdp_dbarts_1b_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      dbarts_1b_sum <- pdp_dbarts_1b_
      mod_simple_1 <- mod_simple_1_
      ## ---- dbarts_pdp_1 plot
      dbarts_1b_sum <- readRDS(
        file = paste0(data_path, "synthetic/dbarts_1b_sum.rds")
      )
      g1 <- 
        dbarts_1b_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "dbarts")) +
        geom_line(data = mod_simple_1,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_1,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <- 
        dbarts_1b_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "dbarts")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_dbarts_1b.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(infl_dbarts_1b_plot_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_dbarts_1b <- mod_dbarts_1b_$mod_dbarts_1b
      ## ---- dbarts_infl_1 plot
      rel_inf <- apply(mod_dbarts_1b$varcount, c(2, 3), sum) |>
        as.data.frame() |>
        mutate(.draw = 1:n()) |>
        pivot_longer(
          cols = -.draw,
          names_to = "var",
          values_to = "count"
        ) |>
        mutate(variable = str_remove(var, "\\.[0-9]*")) |>
        group_by(.draw, variable) |>
        summarise(count = sum(count)) |>
        mutate(rel_inf = count / sum(count)) |>
        ungroup() |>
        dplyr::select(-count) |>
        group_by(variable) |>
        summarise(
          mean = mean(rel_inf),
          median = median(rel_inf),
          lower = quantile(rel_inf, 0.025),
          upper = quantile(rel_inf, 0.975)
        ) |>
        arrange(median) |>
        mutate(variable = factor(variable, levels = unique(variable)))
      g <-
        rel_inf |>
        ggplot(aes(y = variable, x = median)) +
        geom_pointrange(aes(xmin = lower, xmax = upper)) +
        scale_y_discrete("") +
        scale_x_continuous("Relative influence") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_infl_mod_dbarts_1b.png"
        ),
        g,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    
    tar_target(pdp_dbarts_1c_, {
      preds_1c <- mod_dbarts_1b_$preds_1c
      benthos_reefs_sf <- read_all_reefs_data_ |>
        mutate(
          Latitude = st_coordinates(geometry)[, 2],
          Longitude = st_coordinates(geometry)[, 1],
          fYear = as.factor(Year),
        ) |>
        st_drop_geometry() |>
        group_by(Year, fYear, Reef) |>
        summarise(across(c(Latitude, Longitude, CYC, DHW, OTHER), mean))
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- gbm_pdp_1c
      dbarts_1c_sum <- preds_1c |> 
        t() |>
        as_tibble(.name_repair = ~ paste0("V", seq_along(.))) |>
        bind_cols(benthos_reefs_sf) |>
        pivot_longer(
          cols = matches("^V[0-9]*"),
          names_to = ".draw", values_to = "fit"
        ) |>
        group_by(Year, .draw) |>
        summarise(fit = mean(fit)) |>
        ## dplyr::select(Year, fit, .draw) |> 
        group_by(Year, .add = FALSE) |> 
        summarise_draws(median, HDInterval::hdi)
      saveRDS(dbarts_1c_sum,
        file = paste0(data_path, "synthetic/dbarts_1c_sum.rds")
      ) 
      ## ----end
      dbarts_1c_sum
    }),
    tar_target(pdp_dbarts_1c_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      dbarts_1c_sum <- pdp_dbarts_1c_
      mod_simple_1 <- mod_simple_1_
      ## ---- dbarts_pdp_1 plot
      dbarts_1c_sum <- readRDS(
        file = paste0(data_path, "synthetic/dbarts_1c_sum.rds")
      )
      g1 <- 
        dbarts_1c_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "dbarts")) +
        geom_line(data = mod_simple_1,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_1,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <- 
        dbarts_1c_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "dbarts")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_dbarts_1c.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    
    ## Replace a reef (V2) --------------------------------------------
    ## simple ---------------------------------------------------------
    tar_target(mod_simple_2_, {
      benthos_fixed_locs_obs_2 <- site_replacements_data_prep_2_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- simple_2
      mod_simple_2 <- benthos_fixed_locs_obs_2 |>
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
      saveRDS(mod_simple_2,
        file = paste0(data_path, "synthetic/mod_simple_2.rds")
      ) 
      ## ----end
      mod_simple_2
    }),
    ## glmmTMB --------------------------------------------------------
    tar_target(mod_glmmTMB_2_, {
      benthos_fixed_locs_obs_2 <- site_replacements_data_prep_2_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_2
      mod_glmmTMB_2 <- glmmTMB(cover ~ fYear + (1 | Site) + (1 | Transect),
        data = benthos_fixed_locs_obs_2,
        family = "beta_family"
      )
      saveRDS(mod_glmmTMB_2,
        file = paste0(data_path, "synthetic/mod_glmmTMB_2.rds")
      ) 
      ## ----end
      mod_glmmTMB_2
    }),
    tar_target(dharma_mod_glmmTMB_2_, {
      mod_glmmTMB_2 <- mod_glmmTMB_2_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_2_dharma
      glmmTMB_2_dharma <- mod_glmmTMB_2 |> 
        simulateResiduals(n = 1000)
      saveRDS(glmmTMB_2_dharma,
        file = paste0(data_path, "synthetic/glmmTMB_2_dharma.rds")
      ) 
      ## ----end
      glmmTMB_2_dharma
    }),
    tar_target(dharma_mod_glmmTMB_2_plot_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      glmmTMB_2_dharma <- dharma_mod_glmmTMB_2_$glmmTMB_2_dharma
      ## ---- glmmTMB_2_dharma plot
      glmmTMB_2_dharma <- readRDS(
        file = paste0(data_path, "synthetic/glmmTMB_2_dharma.rds")
      )
      g <- wrap_elements(~testUniformity(glmmTMB_2_dharma)) +
        wrap_elements(~plotResiduals(glmmTMB_2_dharma)) +
        wrap_elements(~testDispersion(glmmTMB_2_dharma))
      ggsave(
        filename = paste0(
          fig_path, "R_dharma_mod_glmmTMB_2.png"
        ),
        g,
        width = 10, height = 4, dpi = 72
      )
      ## ----end
    }),
    tar_target(emmeans_mod_glmmTMB_2, {
      mod_glmmTMB_2 <- mod_glmmTMB_2_
      newdata_2 <- site_replacements_newdata_2_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_2_emmeans
      glmmTMB_2_sum <- 
        mod_glmmTMB_2 |> emmeans(~fYear, at = newdata_2, type = "response") |>
        as.data.frame() |> 
        rename(median = response, lower = asymp.LCL, upper = asymp.UCL) |>
        mutate(Year = as.numeric(as.character(fYear))) |>
        mutate(type = "glmmTMB")
      saveRDS(glmmTMB_2_sum,
        file = paste0(data_path, "synthetic/glmmTMB_2_sum.rds")
      ) 
      ## ----end
      glmmTMB_2_sum 
    }),
    tar_target(emmeans_mod_glmmTMB_2_plot_, {
      glmmTMB_2_sum <- emmeans_mod_glmmTMB_2
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_simple_2 <- mod_simple_2_
      ## ---- glmmTMB_2_emmeans plot
      g1 <- glmmTMB_2_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "glmmTMB")) +
        geom_line(data = mod_simple_2,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_2,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <- glmmTMB_2_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "glmmTMB")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_glmmTMB_2.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    
    ## brms -----------------------------------------------------------
    tar_target(mod_brms_2_, {
      benthos_fixed_locs_obs_2 <- site_replacements_data_prep_2_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- brms_pre_2
      benthos_fixed_locs_obs_2 |>
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
      ## ---- brms_2
      mod_brms_2 <- brm(mod_form,
        data = benthos_fixed_locs_obs_2,
        iter = 5000,
        warmup = 1000,
        chains = 3,
        cores = 3,
        prior = priors,
        thin =  5,
        control = list(adapt_delta = 0.99),
        backend = "cmdstanr"
      )
      saveRDS(mod_brms_2,
        file = paste0(data_path, "synthetic/mod_brms_2.rds")
      ) 
      ## ----end
      mod_brms_2
    }),
    tar_target(emmeans_mod_brms_2_, {
      mod_brms_2 <- mod_brms_2_
      newdata_2 <- site_replacements_newdata_2_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- brms_2_emmeans
      brms_2_sum <- 
        mod_brms_2 |> emmeans(~fYear, at = newdata_2, type = "response") |>
        as.data.frame() |> 
        rename(median = response, lower = lower.HPD, upper = upper.HPD) |>
        mutate(Year = as.numeric(as.character(fYear))) |>
        mutate(type = "brms")
      saveRDS(brms_2_sum,
        file = paste0(data_path, "synthetic/brms_2_sum.rds")
      ) 
      ## ----end
      brms_2_sum 
    }),
    tar_target(emmeans_mod_brms_2_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      brms_2_sum <- emmeans_mod_brms_2_
      mod_simple_2 <- mod_simple_2_
      ## ---- brms_2_emmeans plot
      brms_2_sum <- readRDS(
        file = paste0(data_path, "synthetic/brms_2_sum.rds")
      )
      g1 <- 
        brms_2_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "brms")) +
        geom_line(data = mod_simple_2,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_2,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <- 
        brms_2_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "brms")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_brms_2.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_trace_2_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_brms_2 <- mod_brms_2_
      ## ---- brms_trace_2
      mod_brms_2 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_2.rds")
      )
      vars <- mod_brms_2 |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_2$fit |> stan_trace(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_trace_2.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ac_2_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_brms_2 <- mod_brms_2_
      ## ---- brms_ac_2
      mod_brms_2 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_2.rds")
      )
      vars <- mod_brms_2 |>
        brms::variables() |>
        str_subset("^b.*")
      g <- mod_brms_2$fit |> stan_ac(pars = vars)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ac_2.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_rhat_2_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_brms_2 <- mod_brms_2_
      ## ---- brms_rhat_2
      mod_brms_2 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_2.rds")
      )
      g <- mod_brms_2$fit |> stan_rhat()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_rhat_2.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ess_2_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_brms_2 <- mod_brms_2_
      ## ---- brms_ess_2
      mod_brms_2 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_2.rds")
      )
      g <- mod_brms_2$fit |> stan_ess()
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ess_2.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(brms_ppc_2_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_brms_2 <- mod_brms_2_
      ## ---- brms_ppc_2
      mod_brms_2 <- readRDS(
        file = paste0(data_path, "synthetic/mod_brms_2.rds")
      )
      g <- mod_brms_2 |> pp_check( type='dens_overlay', ndraws=100)
      ggsave(
        filename = paste0(
          fig_path, "R_brms_ppc_2.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),

    ## stan -----------------------------------------------------------
    tar_target(mod_stan_2_, {
      source("model_functions.R")
      benthos_fixed_locs_obs_2 <- site_replacements_data_prep_2_
      data_path <- site_replacement_global_parameters_$data_path
      site_extra_functions_
      ## ---- stan_pre_2
      benthos_fixed_locs_obs_2 <-
        benthos_fixed_locs_obs_2 |>
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
      saveRDS(benthos_fixed_locs_obs_2,
        file = paste0(data_path, "synthetic/saveRDS(benthos_fixed_locs_obs_2_forstan.rds")
      ) 
      stan_data <- prepare_data_for_stan(benthos_fixed_locs_obs_2, yrs = NULL)
      model_stan <- cmdstanr::cmdstan_model(stan_file = "model1.stan")
      ## model_stan <- cmdstanr::cmdstan_model(stan_file = "mod1a.stan")
      ## model_stan <- cmdstanr::cmdstan_model(stan_file = "mod2a.stan")
      ## ----end
      ## ---- stan_2
      mod_stan_2 <- model_stan$sample(
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
      saveRDS(mod_stan_2,
        file = paste0(data_path, "synthetic/mod_stan_2.rds")
      ) 
      ## ----end
      mod_stan_2
    }),
    tar_target(pdp_mod_stan_2_, {
      mod_stan_2 <- mod_stan_2_
      newdata_2 <- site_replacements_newdata_2_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- stan_2_pdp
      stan_2_sum <-
        mod_stan_2$draws(variables = "cellmeans") |>
        posterior::as_draws_df() |>
        posterior::summarise_draws(
          median,
          HDInterval::hdi,
          ~ HDInterval::hdi(., credMass = c(0.9)),
          rhat,
          ess_bulk,
          ess_tail
        ) |>
        rename(lower_90 = V4, upper_90 = V5) |>
        bind_cols(newdata_2) |>
        mutate(Year = as.numeric(as.character(fYear)))
      saveRDS(stan_2_sum,
        file = paste0(data_path, "synthetic/stan_2_sum.rds")
      ) 
      ## ----end
      stan_2_sum 
    }),    
    tar_target(pdp_mod_stan_2_plot_, {
      stan_2_sum <- pdp_mod_stan_2_
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_simple_2 <- mod_simple_2_
      ## ---- stan_2_pdp plot
      stan_2_sum <- readRDS(
        file = paste0(data_path, "synthetic/stan_2_sum.rds")
      )
      g1 <- stan_2_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "stan")) +
        geom_line(data = mod_simple_2,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_2,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <- stan_2_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "stan")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_stan_2.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_trace_2_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_stan_2 <- mod_stan_2_
      ## ---- stan_trace_2
      mod_stan_2 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_2.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_2$draws(variables = c("beta", "phi", "sd_2", "sd_2", "sd_3")) |>
        mcmc_trace() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_trace_2.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ac_2_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_stan_2 <- mod_stan_2_
      ## ---- stan_ac_2
      mod_stan_2 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_2.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_2$draws(variables = c("beta", "phi", "sd_2", "sd_2", "sd_3")) |>
        mcmc_acf() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ac_2.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_rhat_2_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_stan_2 <- mod_stan_2_
      ## ---- stan_rhat_2
      mod_stan_2 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_2.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_2 |> bayesplot::rhat() |> 
        mcmc_rhat_hist() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_rhat_2.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ess_2_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_stan_2 <- mod_stan_2_
      ## ---- stan_ess_2
      mod_stan_2 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_2.rds")
      )
      color_scheme_set("viridis")
      g <-
        mod_stan_2 |> bayesplot::neff_ratio() |> 
        mcmc_neff_hist() +
        theme_minimal()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ess_2.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    tar_target(stan_ppc_2_, {
      benthos_fixed_locs_obs_2 <- site_replacements_data_prep_2_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_stan_2 <- mod_stan_2_
      ## ---- stan_ppc_2
      mod_stan_2 <- readRDS(
        file = paste0(data_path, "synthetic/mod_stan_2.rds")
      )
      g <- 
        bayesplot::pp_check(
          benthos_fixed_locs_obs_2$cover,
          mod_stan_2$draws("ypred", format = "matrix")[1:100, ],
          ppc_dens_overlay
        ) +
        theme_classic()
      ggsave(
        filename = paste0(
          fig_path, "R_stan_ppc_2.png"
        ),
        g,
        width = 10, height = 8, dpi = 72
      )
      ## ----end
    }),
    
    ## gbm ------------------------------------------------------------
    tar_target(mod_gbm_2_, {
      benthos_fixed_locs_obs_2 <- site_replacements_data_prep_2_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- gbm_2
      mod_gbm_2 <- gbm(cover ~ fYear,
        data =  benthos_fixed_locs_obs_2,
        distribution = "gaussian",
        n.trees = 10000,
        interaction.depth = 5,
        shrinkage = 0.001,
        bag.fraction = 0.5,
        cv.folds = 5,
        verbose = TRUE
      )
      saveRDS(mod_gbm_2,
        file = paste0(data_path, "synthetic/mod_gbm_2.rds")
      ) 
      ## ----end
      ## ---- gbm_post_2
      n.trees <- gbm.perf(mod_gbm_2, method = "cv")
      ## ----end
      list(mod_gbm_2 = mod_gbm_2, n.trees = n.trees)
    }),
    tar_target(pdp_gbm_2_, {
      mod_gbm_2 <- mod_gbm_2_$mod_gbm_2
      n.trees <- mod_gbm_2_$n.trees
      newdata_2 <- site_replacements_newdata_2_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- gbm_pdp_2
      gbm_2_sum <- newdata_2 |>
        mutate(median = predict(mod_gbm_2, newdata_2, n.trees = n.trees, type = "response")) |> 
        mutate(Year = as.numeric(as.character(fYear))) |>
        mutate(type = "gbm")
      saveRDS(gbm_2_sum,
        file = paste0(data_path, "synthetic/gbm_2_sum.rds")
      ) 
      ## ----end
      gbm_2_sum
    }),
    tar_target(pdp_gbm_2_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      gdm_2_sum <- pdp_gbm_2_
      mod_simple_2 <- mod_simple_2_
      ## ---- gbm_pdp_2
      gbm_2_sum <- readRDS(
        file = paste0(data_path, "synthetic/gbm_2_sum.rds")
      )
      g1 <-
        gbm_2_sum |>
        ggplot() +
        geom_line(aes(x = Year, y = median, color = "gbm")) +
        geom_line(data = mod_simple_2,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_2,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <-
        gbm_2_sum |>
        ggplot() +
        geom_line(aes(x = Year, y = median, color = "gbm")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_gbm_2.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    
    ## gbm + covartiates ----------------------------------------------
    tar_target(mod_gbm_2b_, {
      benthos_fixed_locs_obs_2 <- site_replacements_data_prep_2_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- gbm_2b
      mod_gbm_2b <- gbm(cover ~ fYear + Latitude + Longitude + CYC + DHW + OTHER,
        data =  benthos_fixed_locs_obs_2,
        distribution = "gaussian",
        var.monotone = c(0, 0, 0, -1, -1, -1),
        n.trees = 10000,
        interaction.depth = 5,
        shrinkage = 0.001,
        bag.fraction = 0.5,
        cv.folds = 5,
        verbose = TRUE
      )
      saveRDS(mod_gbm_2b,
        file = paste0(data_path, "synthetic/mod_gbm_2b.rds")
      )
      ## ----end
      ## ---- gbm_post_2b
      n.trees <- gbm.perf(mod_gbm_2b, method = "cv")
      ## ----end
      list(mod_gbm_2b = mod_gbm_2b, n.trees = n.trees)
    }),
    tar_target(pdp_gbm_2b_, {
      benthos_fixed_locs_obs_2 <- site_replacements_data_prep_2_
      mod_gbm_2b <- mod_gbm_2b_$mod_gbm_2b
      n.trees <- mod_gbm_2b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- gbm_pdp_2b
      newdata_2b <- benthos_fixed_locs_obs_2
      gbm_2b_sum <- newdata_2b |>
        mutate(median = predict(mod_gbm_2b, newdata_2b, n.trees = n.trees, type = "response")) |>
        mutate(Year = as.numeric(as.character(fYear))) |>
        ## group_by(Year, Site, Transect) |>
        ## summarise(median = median(median)) |>
        ## group_by(Year, Site, .add = FALSE) |>
        ## summarise(median = median(median)) |>
        group_by(Year, .add = FALSE) |> 
        summarise(median = median(median)) |> 
        mutate(type = "gbm")
      saveRDS(gbm_2b_sum,
        file = paste0(data_path, "synthetic/gbm_2b_sum.rds")
      ) 
      ## ----end
      gbm_2b_sum
    }),
    tar_target(pdp_gbm_2b_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      gbm_2b_sum <- pdp_gbm_2b_
      mod_simple_2 <- mod_simple_2_
      ## ---- gbm_pdp_2b plot
      gbm_2b_sum <- readRDS(
        file = paste0(data_path, "synthetic/gbm_2b_sum.rds")
      )
      g1 <-
        gbm_2b_sum |>
        ggplot() +
        geom_line(aes(x = Year, y = median, color = "gbm")) +
        geom_line(data = mod_simple_2,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_2,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <-
        gbm_2b_sum |>
        ggplot() +
        geom_line(aes(x = Year, y = median, color = "gbm")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_gbm_2b.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(infl_gbm_2b_, {
      mod_gbm_2b <- mod_gbm_2b_$mod_gbm_2b
      n.trees <- mod_gbm_2b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- gbm_infl_2b
      ## gbm_2b_infl <- gbm::relative.influence(mod_gbm_2b, n.trees = n.trees, scale = TRUE, sort = TRUE)
      infl <- summary(mod_gbm_2b, n.trees =  n.trees, plot = FALSE)
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
          fig_path, "R_infl_mod_gbm_2b.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(pdp_gbm_2c_, {
      benthos_reefs_sf <- read_all_reefs_data_
      mod_gbm_2b <- mod_gbm_2b_$mod_gbm_2b
      n.trees <- mod_gbm_2b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- gbm_pdp_2c
      newdata_2c <- benthos_reefs_sf |>
        mutate(
          Latitude = st_coordinates(geometry)[,2],
          Longitude = st_coordinates(geometry)[,1],
          fYear = as.factor(Year),
          )
      gbm_2c_sum <- newdata_2c |>
        mutate(median = predict(mod_gbm_2b, newdata_2c, n.trees = n.trees, type = "response")) |>
        mutate(Year = as.numeric(as.character(fYear))) |>
        ## group_by(Year, Site, Transect) |>
        ## summarise(median = median(median)) |>
        ## group_by(Year, Site, .add = FALSE) |>
        ## summarise(median = median(median)) |>
        group_by(Year, .add = FALSE) |> 
        summarise(median = median(median)) |> 
        mutate(type = "gbm")
      saveRDS(gbm_2c_sum,
        file = paste0(data_path, "synthetic/gbm_2c_sum.rds")
      ) 
      ## ----end
      gbm_2c_sum
    }),
    tar_target(pdp_gbm_2c_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      gbm_2c_sum <- pdp_gbm_2c_
      mod_simple_2 <- mod_simple_2_
      ## ---- gbm_pdp_2c plot
      gbm_2c_sum <- readRDS(
        file = paste0(data_path, "synthetic/gbm_2c_sum.rds")
      )
      g1 <-
        gbm_2c_sum |>
        ggplot() +
        geom_line(aes(x = Year, y = median, color = "gbm")) +
        geom_line(data = mod_simple_2,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_2,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <-
        gbm_2c_sum |>
        ggplot() +
        geom_line(aes(x = Year, y = median, color = "gbm")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_gbm_2c.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(infl_gbm_2c_, {
      mod_gbm_2b <- mod_gbm_2b_$mod_gbm_2b
      n.trees <- mod_gbm_2b_$n.trees
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      ## ---- gbm_infl_2c
      ## gbm_2c_infl <- gbm::relative.influence(mod_gbm_2c, n.trees = n.trees, scale = TRUE, sort = TRUE)
      infl <- summary(mod_gbm_2b, n.trees =  n.trees, plot = FALSE)
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
          fig_path, "R_infl_mod_gbm_2c.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
    }),

    ## dbarts --------------------------------------------------------
    tar_target(mod_dbarts_2_, {
      benthos_fixed_locs_obs_2 <- site_replacements_data_prep_2_
      data_path <- site_replacement_global_parameters_$data_path
      newdata_2 <- site_replacements_newdata_2_
      ## ---- dbarts_2
      print(head(benthos_fixed_locs_obs_2))
      mod_dbarts_2 <- bart2(log(cover) ~ fYear,
        data =   benthos_fixed_locs_obs_2,
        keepTrees = TRUE
      )
      saveRDS(mod_dbarts_2,
        file = paste0(data_path, "synthetic/mod_dbarts_2.rds")
      ) 
      ## ----end
      ## Unfortunately, the next part must be in the same tar_target
      ## due to the way dbarts stores pointers - they cannot be stored
      ## ---- dbarts_pred_2
      preds <- predict(mod_dbarts_2, newdata_2, type = "ev") |>
        exp() |> 
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds,
        file = paste0(data_path, "synthetic/mod_dbarts_2_preds.rds")
      ) 
      ## ----end
      list(
        mod_dbarts_2 = mod_dbarts_2,
        preds = preds
      )
    }),
    tar_target(dbarts_pdp_2, {
      data_path <- site_replacement_global_parameters_$data_path
      newdata_2 <- site_replacements_newdata_2_
      preds <- mod_dbarts_2_$preds
      ## ---- dbarts_pdp_2
      dbarts_2_sum <- newdata_2 |>
        bind_cols(preds) |> 
        mutate(Year = as.numeric(as.character(fYear))) |>
        mutate(type = "dbarts")
      saveRDS(dbarts_2_sum,
        file = paste0(data_path, "synthetic/dbarts_2_sum.rds")
      ) 
      # ----end
      dbarts_2_sum 
    }),
    tar_target(pdp_dbarts_2_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      dbarts_2_sum <- dbarts_pdp_2
      mod_simple_2 <- mod_simple_2_
      ## ---- dbarts_pdp_2 plot
      dbarts_2_sum <- readRDS(
        file = paste0(data_path, "synthetic/dbarts_2_sum.rds")
      )
      g1 <- 
        dbarts_2_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "dbarts")) +
        geom_line(data = mod_simple_2,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_2,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <- 
        dbarts_2_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "dbarts")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_dbarts_2.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    
    ## dbarts + covariates ---------------------------------------------
    tar_target(mod_dbarts_2b_, {
      benthos_fixed_locs_obs_2 <- site_replacements_data_prep_2_
      benthos_reefs_sf <- read_all_reefs_data_
      data_path <- site_replacement_global_parameters_$data_path
      newdata_2 <- site_replacements_newdata_2_
      ## ---- dbarts_2b
      mod_dbarts_2b <- bart2(log(cover) ~ fYear + Latitude + Longitude + CYC + DHW + OTHER,
        data =   benthos_fixed_locs_obs_2,
        keepTrees = TRUE
      )
      saveRDS(mod_dbarts_2b,
        file = paste0(data_path, "synthetic/mod_dbarts_2b.rds")
      ) 
      ## ----end
      ## Unfortunately, the next part must be in the same tar_target
      ## due to the way dbarts stores pointers - they cannot be stored
      ## ---- dbarts_pred_2b
      newdata_2b <- benthos_fixed_locs_obs_2
      preds <- predict(mod_dbarts_2b, newdata_2b, type = "ev") |>
        exp()
      saveRDS(preds,
        file = paste0(data_path, "synthetic/mod_dbarts_2b_preds.rds")
      ) 
      ## ----end
      ## ---- dbarts_pred_2c
      newdata_2c <- benthos_reefs_sf |>
        mutate(
          Latitude = st_coordinates(geometry)[, 2],
          Longitude = st_coordinates(geometry)[, 1],
          fYear = as.factor(Year),
        ) |>
        st_drop_geometry() |>
        group_by(Year, fYear, Reef) |>
        summarise(across(c(Latitude, Longitude, CYC, DHW, OTHER), mean))
      preds_2c <- predict(mod_dbarts_2b, newdata_2c, type = "ev") |>
        exp()
      saveRDS(preds_2c,
        file = paste0(data_path, "synthetic/mod_dbarts_2c_preds.rds")
      ) 
      ## ----end
      list(
        mod_dbarts_2b = mod_dbarts_2b,
        preds = preds,
        preds_2c = preds_2c
      )
    }),
    tar_target(pdp_dbarts_2b_, {
      preds <- mod_dbarts_2b_$preds
      data_path <- site_replacement_global_parameters_$data_path
      benthos_fixed_locs_obs_2 <- site_replacements_data_prep_2_
      ## ---- dbarts_pdp_2c
      dbarts_2b_sum <- preds |> 
        t() |>
        as_tibble(.name_repair = ~ paste0("V", seq_along(.))) |>
        bind_cols(benthos_fixed_locs_obs_2) |>
        pivot_longer(
          cols = matches("^V[0-9]*"),
          names_to = ".draw", values_to = "fit"
        ) |>
        group_by(Year, .draw) |>
        summarise(fit = mean(fit)) |>
        ## dplyr::select(Year, fit, .draw) |> 
        group_by(Year, .add = FALSE) |> 
        summarise_draws(median, HDInterval::hdi)
      saveRDS(dbarts_2b_sum,
        file = paste0(data_path, "synthetic/dbarts_2b_sum.rds")
      ) 
      ## ----end
      dbarts_2b_sum
    }),
    tar_target(pdp_dbarts_2b_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      dbarts_2b_sum <- pdp_dbarts_2b_
      mod_simple_2 <- mod_simple_2_
      ## ---- dbarts_pdp_2c plot
      dbarts_2b_sum <- readRDS(
        file = paste0(data_path, "synthetic/dbarts_2b_sum.rds")
      )
      g1 <- 
        dbarts_2b_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "dbarts")) +
        geom_line(data = mod_simple_2,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_2,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <- 
        dbarts_2b_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "dbarts")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_dbarts_2b.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(pdp_dbarts_2c_, {
      preds_2c <- mod_dbarts_2b_$preds_2c
      benthos_reefs_sf <- read_all_reefs_data_ |>
        mutate(
          Latitude = st_coordinates(geometry)[, 2],
          Longitude = st_coordinates(geometry)[, 1],
          fYear = as.factor(Year),
        ) |>
        st_drop_geometry() |>
        group_by(Year, fYear, Reef) |>
        summarise(across(c(Latitude, Longitude, CYC, DHW, OTHER), mean))
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- gbm_pdp_2c
      dbarts_2c_sum <- preds_2c |> 
        t() |>
        as_tibble(.name_repair = ~ paste0("V", seq_along(.))) |>
        bind_cols(benthos_reefs_sf) |>
        pivot_longer(
          cols = matches("^V[0-9]*"),
          names_to = ".draw", values_to = "fit"
        ) |>
        group_by(Year, .draw) |>
        summarise(fit = mean(fit)) |>
        ## dplyr::select(Year, fit, .draw) |> 
        group_by(Year, .add = FALSE) |> 
        summarise_draws(median, HDInterval::hdi)
      saveRDS(dbarts_2c_sum,
        file = paste0(data_path, "synthetic/dbarts_2c_sum.rds")
      ) 
      ## ----end
      dbarts_2c_sum
    }),
    tar_target(pdp_dbarts_2c_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      dbarts_2c_sum <- pdp_dbarts_2c_
      mod_simple_2 <- mod_simple_2_
      ## ---- dbarts_pdp_2 plot
      dbarts_2c_sum <- readRDS(
        file = paste0(data_path, "synthetic/dbarts_2c_sum.rds")
      )
      g1 <- 
        dbarts_2c_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "dbarts")) +
        geom_line(data = mod_simple_2,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_2,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <- 
        dbarts_2c_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "dbarts")) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_dbarts_2c.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(infl_dbarts_2b_plot_, {
      data_path <- site_replacement_global_parameters_$data_path
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_dbarts_2b <- mod_dbarts_2b_$mod_dbarts_2b
      ## ---- dbarts_infl_1 plot
      rel_inf <- apply(mod_dbarts_2b$varcount, c(2, 3), sum) |>
        as.data.frame() |>
        mutate(.draw = 1:n()) |>
        pivot_longer(
          cols = -.draw,
          names_to = "var",
          values_to = "count"
        ) |>
        mutate(variable = str_remove(var, "\\.[0-9]*")) |>
        group_by(.draw, variable) |>
        summarise(count = sum(count)) |>
        mutate(rel_inf = count / sum(count)) |>
        ungroup() |>
        dplyr::select(-count) |>
        group_by(variable) |>
        summarise(
          mean = mean(rel_inf),
          median = median(rel_inf),
          lower = quantile(rel_inf, 0.025),
          upper = quantile(rel_inf, 0.975)
        ) |>
        arrange(median) |>
        mutate(variable = factor(variable, levels = unique(variable)))
      g <-
        rel_inf |>
        ggplot(aes(y = variable, x = median)) +
        geom_pointrange(aes(xmin = lower, xmax = upper)) +
        scale_y_discrete("") +
        scale_x_continuous("Relative influence") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_infl_mod_dbarts_2b.png"
        ),
        g,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    })
    

    ## tar_target(mod_stan_1_, {
    ##   benthos_fixed_locs_obs_1 <- site_replacements_data_prep_1_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   site_extra_functions_
    ##   ## ---- stan_pre_1
    ##   benthos_fixed_locs_obs_1 <-
    ##     benthos_fixed_locs_obs_1 |>
    ##     mutate(
    ##       cYear = fYear,
    ##       grid_id = factor(Reef),
    ##       cSite = factor(interaction(Reef, Site)),
    ##       cReplicate = ifelse(is.na(Transect),
    ##         interaction(Site, Year),
    ##         interaction(cSite, Transect)),
    ##       Cover = cover,
    ##       area = 1,
    ##       sum = 1
    ##     )
    ##   saveRDS(benthos_fixed_locs_obs_1,
    ##     file = paste0(data_path, "synthetic/saveRDS(benthos_fixed_locs_obs_1_forstan.rds")
    ##   ) 
    ##   stan_data <- prepare_data_for_stan(benthos_fixed_locs_obs_1, yrs = NULL)
    ##   model_stan <- cmdstanr::cmdstan_model(stan_file = "model1.stan")
    ##   ## model_stan <- cmdstanr::cmdstan_model(stan_file = "mod1a.stan")
    ##   ## model_stan <- cmdstanr::cmdstan_model(stan_file = "mod2a.stan")
    ##   ## ----end
    ##   ## ---- stan_1
    ##   mod_stan_1 <- model_stan$sample(
    ##     data = stan_data,
    ##     seed = 123,
    ##     iter_sampling = 5000,
    ##     iter_warmup = 1000,
    ##     thin = 5,
    ##     chains = 3,
    ##     parallel_chains = 3,
    ##     adapt_delta = 0.99,
    ##     output_dir = paste0(data_path, "synthetic/"),
    ##   )
    ##   saveRDS(mod_stan_1,
    ##     file = paste0(data_path, "synthetic/mod_stan_1.rds")
    ##   ) 
    ##   ## ----end
    ##   mod_stan_1
    ## }),
    ## tar_target(pdp_mod_stan_1, {
    ##   mod_stan_1 <- mod_stan_1_
    ##   newdata_1 <- site_replacements_newdata_1_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- stan_1_pdp
    ##   stan_1_sum <-
    ##     mod_stan_1$draws(variables = "cellmeans") |>
    ##     posterior::as_draws_df() |>
    ##     posterior::summarise_draws(
    ##       median,
    ##       HDInterval::hdi,
    ##       ~ HDInterval::hdi(., credMass = c(0.9)),
    ##       rhat,
    ##       ess_bulk,
    ##       ess_tail
    ##     ) |>
    ##     rename(lower_90 = V4, upper_90 = V5) |>
    ##     bind_cols(newdata_1) |>
    ##     mutate(Year = as.numeric(as.character(fYear)))
    ##   saveRDS(stan_1_sum,
    ##     file = paste0(data_path, "synthetic/stan_1_sum.rds")
    ##   ) 
    ##   ## ----end
    ##   stan_1_sum 
    ## }),    
    
    ## ## gbm
    ## tar_target(mod_gbm_1_, {
    ##   benthos_fixed_locs_obs_1 <- site_replacements_data_prep_1_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- gbm_1
    ##   mod_gbm_1 <- gbm(cover ~ fYear,
    ##     data =  benthos_fixed_locs_obs_1,
    ##     distribution = "gaussian",
    ##     n.trees = 10000,
    ##     interaction.depth = 5,
    ##     shrinkage = 0.001,
    ##     bag.fraction = 0.5,
    ##     cv.folds = 5,
    ##     verbose = TRUE
    ##   )
    ##   saveRDS(mod_gbm_1,
    ##     file = paste0(data_path, "synthetic/mod_gbm_1.rds")
    ##   ) 
    ##   ## ----end
    ##   ## ---- gbm_post_1
    ##   n.trees <- gbm.perf(mod_gbm_1, method = "cv")
    ##   ## ----end
    ##   list(mod_gbm_1 = mod_gbm_1, n.trees = n.trees)
    ## }),
    ## tar_target(pdp_gbm_1_, {
    ##   mod_gbm_1 <- mod_gbm_1_$mod_gbm_1
    ##   n.trees <- mod_gbm_1_$n.trees
    ##   newdata_1 <- site_replacements_newdata_1_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- gbm_pdp_1
    ##   gbm_1_sum <- newdata_1 |>
    ##     mutate(median = predict(mod_gbm_1, newdata_1, n.trees = n.trees, type = "response")) |> 
    ##     mutate(Year = as.numeric(as.character(fYear))) |>
    ##     mutate(type = "gbm")
    ##   saveRDS(gbm_1_sum,
    ##     file = paste0(data_path, "synthetic/gbm_1_sum.rds")
    ##   ) 
    ##   ## ----end
    ## }),
    ## ## gbm 2
    ## tar_target(mod_gbm2_1_, {
    ##   benthos_fixed_locs_obs_1 <- site_replacements_data_prep_1_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- gbm2_1
    ##   mod_gbm2_1 <- gbm(cover ~ fYear + Latitude + Longitude + Reef + Site + Transect,
    ##     data =  benthos_fixed_locs_obs_1,
    ##     distribution = "gaussian",
    ##     n.trees = 10000,
    ##     interaction.depth = 5,
    ##     shrinkage = 0.001,
    ##     bag.fraction = 0.5,
    ##     cv.folds = 5,
    ##     verbose = TRUE
    ##   )
    ##   saveRDS(mod_gbm2_1,
    ##     file = paste0(data_path, "synthetic/mod_gbm2_1.rds")
    ##   ) 
    ##   ## ----end
    ##   ## ---- gbm2_post_1
    ##   n.trees <- gbm.perf(mod_gbm2_1, method = "cv")
    ##   ## ----end
    ##   list(mod_gbm2_1 = mod_gbm2_1, n.trees = n.trees)
    ## }),
    ## tar_target(pdp_gbm2_1_, {
    ##   mod_gbm2_1 <- mod_gbm2_1_$mod_gbm2_1
    ##   n.trees <- mod_gbm2_1_$n.trees
    ##   newdata2_1 <- site_replacements_newdata2_1_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- gbm2_pdp_1
    ##   gbm2_1_sum <- newdata2_1 |>
    ##     mutate(median = predict(mod_gbm2_1, newdata2_1, n.trees = n.trees, type = "response")) |>
    ##     group_by(fYear) |>
    ##     summarise(
    ##       median = median(median)
    ##     ) |>
    ##     ungroup() |> 
    ##     mutate(Year = as.numeric(as.character(fYear))) |>
    ##     mutate(type = "gbm2")
    ##   saveRDS(gbm2_1_sum,
    ##     file = paste0(data_path, "synthetic/gbm2_1_sum.rds")
    ##   ) 
    ##   ## ----end
    ## }),


    ## ## dbarts
    ## tar_target(mod_dbarts_1_, {
    ##   benthos_fixed_locs_obs_1 <- site_replacements_data_prep_1_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- dbarts_1
    ##   mod_dbarts_1 <- bart2(log(cover) ~ fYear,
    ##     data =   benthos_fixed_locs_obs_1,
    ##     keepTrees = TRUE
    ##   )
    ##   saveRDS(mod_dbarts_1,
    ##     file = paste0(data_path, "synthetic/mod_dbarts_1.rds")
    ##   ) 
    ##   ## ----end
    ##   mod_dbarts_1
    ## }),
    ## tar_target(pdp_dbarts_1_, {
    ##   mod_dbarts_1 <- mod_dbarts_1_
    ##   newdata_1 <- site_replacements_newdata_1_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- dbarts_pdp_1
    ##   preds <- predict(mod_dbarts_1, newdata_1, type = "ev") |>
    ##     exp() |> 
    ##     summarise_draws(median, HDInterval::hdi)
    ##   dbarts_1_sum <- newdata_1 |>
    ##     bind_cols(preds) |> 
    ##     mutate(Year = as.numeric(as.character(fYear))) |>
    ##     mutate(type = "dbarts")
    ##   saveRDS(dbarts_1_sum,
    ##     file = paste0(data_path, "synthetic/dbarts_1_sum.rds")
    ##   ) 
    ##   ## ----end
    ## }),

    ## ## dbarts 2
    ## tar_target(mod_dbarts2_1_, {
    ##   benthos_fixed_locs_obs_1 <- site_replacements_data_prep_1_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- dbarts2_1
    ##   mod_dbarts2_1 <- bart2(log(cover) ~ fYear + Latitude + Longitude + Reef + Site + Transect,
    ##     data =   benthos_fixed_locs_obs_1,
    ##     keepTrees = TRUE
    ##   )
    ##   saveRDS(mod_dbarts2_1,
    ##     file = paste0(data_path, "synthetic/mod_dbarts2_1.rds")
    ##   ) 
    ##   ## ----end
    ##   mod_dbarts2_1
    ## }),
    ## tar_target(pdp_dbarts2_1_, {
    ##   mod_dbarts2_1 <- mod_dbarts2_1_
    ##   newdata2_1 <- site_replacements_newdata2_1_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- dbarts2_pdp_1
    ##   preds <- predict(mod_dbarts2_1, newdata2_1, type = "ev") |>
    ##     exp()
    ##   ## summarise_draws(median, HDInterval::hdi)
    ##   dbarts2_1_sum <-
    ##     preds |>
    ##     as_tibble() |>
    ##     mutate(.draw = 1:n()) |>
    ##     pivot_longer(cols = -.draw, names_to = "Variable", values_to = "values") |>
    ##     left_join(
    ##       newdata2_1 |>
    ##         mutate(Variable = paste0("V", 1:n())),
    ##       by = c("Variable")
    ##     ) |>
    ##     group_by(.draw, fYear) |>
    ##     summarise(
    ##       median = median(values)
    ##     ) |>
    ##     ungroup() |>
    ##     group_by(fYear) |>
    ##     tidybayes::summarise_draws(median, HDInterval::hdi) |> 
    ##     mutate(Year = as.numeric(as.character(fYear))) |>
    ##     mutate(type = "dbarts2")
    ##   saveRDS(dbarts2_1_sum,
    ##     file = paste0(data_path, "synthetic/dbarts2_1_sum.rds")
    ##   ) 
    ##   ## ----end
    ## }),

    
    ## ## Type 2 ---------------------------------------------------------
    ## ## glmmTMB --------------------------------------------------------
    ## tar_target(mod_simple_2_, {
    ##   benthos_fixed_locs_obs_2 <- site_replacements_data_prep_2_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- simple_2
    ##   mod_simple_2 <- benthos_fixed_locs_obs_2 |>
    ##     group_by(Year, Site) |>
    ##     summarise(
    ##       Mean = mean(cover),
    ##       Median = median(cover)
    ##     ) |>
    ##     ungroup() |>
    ##     group_by(Year) |> 
    ##     summarise(
    ##       Mean = mean(Mean),
    ##       Median = median(Median)
    ##     ) 
    ##   saveRDS(mod_simple_2,
    ##     file = paste0(data_path, "synthetic/mod_simple_2.rds")
    ##   ) 
    ##   ## ----end
    ##   mod_simple_2
    ## }),
    ## tar_target(mod_glmmTMB_2_, {
    ##   benthos_fixed_locs_obs_2 <- site_replacements_data_prep_2_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- glmmTMB_2
    ##   mod_glmmTMB_2 <- glmmTMB(cover ~ fYear + (1 | Site) + (1 | Transect),
    ##     data = benthos_fixed_locs_obs_2,
    ##     family = "beta_family"
    ##   )
    ##   saveRDS(mod_glmmTMB_2,
    ##     file = paste0(data_path, "synthetic/mod_glmmTMB_2.rds")
    ##   ) 
    ##   ## ----end
    ##   mod_glmmTMB_2
    ## }),
    ## tar_target(dharma_mod_glmmTMB_2, {
    ##   mod_glmmTMB_2 <- mod_glmmTMB_2_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- glmmTMB_2_dharma
    ##   glmmTMB_2_dharma <- mod_glmmTMB_2 |> 
    ##     simulateResiduals(n = 1000)
    ##   saveRDS(glmmTMB_2_dharma,
    ##     file = paste0(data_path, "synthetic/glmmTMB_2_dharma.rds")
    ##   ) 
    ##   ## ----end
    ##   glmmTMB_2_dharma
    ## }),
    ## tar_target(emmeans_mod_glmmTMB_2, {
    ##   mod_glmmTMB_2 <- mod_glmmTMB_2_
    ##   newdata_2 <- site_replacements_newdata_2_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- glmmTMB_2_emmeans
    ##   glmmTMB_2_sum <- 
    ##     mod_glmmTMB_2 |> emmeans(~fYear, at = newdata_2, type = "response") |>
    ##     as.data.frame() |> 
    ##     rename(median = response, lower = asymp.LCL, upper = asymp.UCL) |>
    ##     mutate(Year = as.numeric(as.character(fYear))) |>
    ##     mutate(type = "glmmTMB")
    ##   saveRDS(glmmTMB_2_sum,
    ##     file = paste0(data_path, "synthetic/glmmTMB_2_sum.rds")
    ##   ) 
    ##   ## ----end
    ##   glmmTMB_2_sum 
    ## }),

    ## ## brms -----------------------------------------------------------
    ## tar_target(mod_brms_2_, {
    ##   benthos_fixed_locs_obs_2 <- site_replacements_data_prep_2_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- brms_pre_2
    ##   benthos_fixed_locs_obs_2 |>
    ##     group_by(fYear) |>
    ##     summarise(mean_abundance = qlogis(mean(cover)),
    ##       median_abundance = qlogis(median(cover)),
    ##       sd_abundance = sd(qlogis(cover)) ,
    ##       mad_abundance = mad(qlogis(cover)) 
    ##     ) 
    ##   priors <- c(
    ##     prior(normal(0, 0.5), class = "b"),
    ##     prior(normal(-1.4, 0.4), class = "Intercept"),
    ##     prior(student_t(3, 0, 0.4), class = "sd")
    ##   )
    ##   mod_form <- bf(cover ~ fYear + (1 | Site) + (1 | Transect),
    ##     family = "Beta"
    ##   )
    ##   ## ----end
    ##   ## ---- brms_2
    ##   mod_brms_2 <- brm(mod_form,
    ##     data = benthos_fixed_locs_obs_2,
    ##     iter = 5000,
    ##     warmup = 1000,
    ##     chains = 3,
    ##     cores = 3,
    ##     prior = priors,
    ##     thin =  5,
    ##     control = list(adapt_delta = 0.99),
    ##     backend = "cmdstanr"
    ##   )
    ##   saveRDS(mod_brms_2,
    ##     file = paste0(data_path, "synthetic/mod_brms_2.rds")
    ##   ) 
    ##   ## ----end
    ##   mod_brms_2
    ## }),
    ## tar_target(emmeans_mod_brms_2, {
    ##   mod_brms_2 <- mod_brms_2_
    ##   newdata_2 <- site_replacements_newdata_2_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- brms_2_emmeans
    ##   brms_2_sum <- 
    ##     mod_brms_2 |> emmeans(~fYear, at = newdata_2, type = "response") |>
    ##     as.data.frame() |> 
    ##     rename(median = response, lower = lower.HPD, upper = upper.HPD) |>
    ##     mutate(Year = as.numeric(as.character(fYear))) |>
    ##     mutate(type = "brms")
    ##   saveRDS(brms_2_sum,
    ##     file = paste0(data_path, "synthetic/brms_2_sum.rds")
    ##   ) 
    ##   ## ----end
    ##   brms_2_sum 
    ## }),

    ## ## stan
    ## tar_target(mod_stan_2_, {
    ##   benthos_fixed_locs_obs_2 <- site_replacements_data_prep_2_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   site_extra_functions_
    ##   ## ---- stan_pre_2
    ##   benthos_fixed_locs_obs_2 <-
    ##     benthos_fixed_locs_obs_2 |>
    ##     mutate(
    ##       cYear = fYear,
    ##       grid_id = factor(Reef),
    ##       cSite = factor(interaction(Reef, Site)),
    ##       cReplicate = ifelse(is.na(Transect),
    ##         interaction(Site, Year),
    ##         interaction(cSite, Transect)),
    ##       Cover = cover,
    ##       area = 1,
    ##       sum = 1
    ##     )
    ##   saveRDS(benthos_fixed_locs_obs_2,
    ##     file = paste0(data_path, "synthetic/saveRDS(benthos_fixed_locs_obs_2_forstan.rds")
    ##   ) 
    ##   stan_data <- prepare_data_for_stan(benthos_fixed_locs_obs_2, yrs = NULL)
    ##   model_stan <- cmdstanr::cmdstan_model(stan_file = "model1.stan")
    ##   ## model_stan <- cmdstanr::cmdstan_model(stan_file = "mod1a.stan")
    ##   ## model_stan <- cmdstanr::cmdstan_model(stan_file = "mod2a.stan")
    ##   ## ----end
    ##   ## ---- stan_2
    ##   mod_stan_2 <- model_stan$sample(
    ##     data = stan_data,
    ##     seed = 123,
    ##     iter_sampling = 5000,
    ##     iter_warmup = 1000,
    ##     thin = 5,
    ##     chains = 3,
    ##     parallel_chains = 3,
    ##     adapt_delta = 0.99,
    ##     output_dir = paste0(data_path, "synthetic/"),
    ##   )
    ##   saveRDS(mod_stan_2,
    ##     file = paste0(data_path, "synthetic/mod_stan_2.rds")
    ##   ) 
    ##   ## ----end
    ##   mod_stan_2
    ## }),
    ## tar_target(pdp_mod_stan_2, {
    ##   mod_stan_2 <- mod_stan_2_
    ##   newdata_2 <- site_replacements_newdata_2_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- stan_2_pdp
    ##   stan_2_sum <-
    ##     mod_stan_2$draws(variables = "cellmeans") |>
    ##     posterior::as_draws_df() |>
    ##     posterior::summarise_draws(
    ##       median,
    ##       HDInterval::hdi,
    ##       ~ HDInterval::hdi(., credMass = c(0.9)),
    ##       rhat,
    ##       ess_bulk,
    ##       ess_tail
    ##     ) |>
    ##     rename(lower_90 = V4, upper_90 = V5) |>
    ##     bind_cols(newdata_2) |>
    ##     mutate(Year = as.numeric(as.character(fYear)))
    ##   saveRDS(stan_2_sum,
    ##     file = paste0(data_path, "synthetic/stan_2_sum.rds")
    ##   ) 
    ##   ## ----end
    ##   stan_2_sum 
    ## }),    

    ## ## gbm
    ## tar_target(mod_gbm_2_, {
    ##   benthos_fixed_locs_obs_2 <- site_replacements_data_prep_2_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- gbm_2
    ##   mod_gbm_2 <- gbm(cover ~ fYear,
    ##     data =  benthos_fixed_locs_obs_2,
    ##     distribution = "gaussian",
    ##     n.trees = 10000,
    ##     interaction.depth = 5,
    ##     shrinkage = 0.001,
    ##     bag.fraction = 0.5,
    ##     cv.folds = 5,
    ##     verbose = TRUE
    ##   )
    ##   saveRDS(mod_gbm_2,
    ##     file = paste0(data_path, "synthetic/mod_gbm_2.rds")
    ##   ) 
    ##   ## ----end
    ##   ## ---- gbm_post_2
    ##   n.trees <- gbm.perf(mod_gbm_2, method = "cv")
    ##   ## ----end
    ##   list(mod_gbm_2 = mod_gbm_2, n.trees = n.trees)
    ## }),
    ## tar_target(pdp_gbm_2_, {
    ##   mod_gbm_2 <- mod_gbm_2_$mod_gbm_2
    ##   n.trees <- mod_gbm_2_$n.trees
    ##   newdata_2 <- site_replacements_newdata_2_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- gbm_pdp_2
    ##   gbm_2_sum <- newdata_2 |>
    ##     mutate(median = predict(mod_gbm_2, newdata_2, n.trees = n.trees, type = "response")) |> 
    ##     mutate(Year = as.numeric(as.character(fYear))) |>
    ##     mutate(type = "gbm")
    ##   saveRDS(gbm_2_sum,
    ##     file = paste0(data_path, "synthetic/gbm_2_sum.rds")
    ##   ) 
    ##   ## ----end
    ## }),

    ## ## gbm 2
    ## tar_target(mod_gbm2_2_, {
    ##   benthos_fixed_locs_obs_2 <- site_replacements_data_prep_2_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- gbm2_2
    ##   mod_gbm2_2 <- gbm(cover ~ fYear + Latitude + Longitude + Reef + Site + Transect,
    ##     data =  benthos_fixed_locs_obs_2,
    ##     distribution = "gaussian",
    ##     n.trees = 10000,
    ##     interaction.depth = 5,
    ##     shrinkage = 0.001,
    ##     bag.fraction = 0.5,
    ##     cv.folds = 5,
    ##     verbose = TRUE
    ##   )
    ##   saveRDS(mod_gbm2_2,
    ##     file = paste0(data_path, "synthetic/mod_gbm2_2.rds")
    ##   ) 
    ##   ## ----end
    ##   ## ---- gbm2_post_2
    ##   n.trees <- gbm.perf(mod_gbm2_2, method = "cv")
    ##   ## ----end
    ##   list(mod_gbm2_2 = mod_gbm2_2, n.trees = n.trees)
    ## }),
    ## tar_target(pdp_gbm2_2_, {
    ##   mod_gbm2_2 <- mod_gbm2_2_$mod_gbm2_2
    ##   n.trees <- mod_gbm2_2_$n.trees
    ##   newdata2_2 <- site_replacements_newdata2_2_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- gbm2_pdp_2
    ##   gbm2_2_sum <- newdata2_2 |>
    ##     mutate(median = predict(mod_gbm2_2, newdata2_2, n.trees = n.trees, type = "response")) |>
    ##     group_by(fYear) |>
    ##     summarise(
    ##       median = median(median)
    ##     ) |>
    ##     ungroup() |> 
    ##     mutate(Year = as.numeric(as.character(fYear))) |>
    ##     mutate(type = "gbm2")
    ##   saveRDS(gbm2_2_sum,
    ##     file = paste0(data_path, "synthetic/gbm2_2_sum.rds")
    ##   ) 
    ##   ## ----end
    ## }),

    ## ## dbarts
    ## tar_target(mod_dbarts_2_, {
    ##   benthos_fixed_locs_obs_2 <- site_replacements_data_prep_2_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- dbarts_2
    ##   mod_dbarts_2 <- bart2(log(cover) ~ fYear,
    ##     data =   benthos_fixed_locs_obs_2,
    ##     keepTrees = TRUE
    ##   )
    ##   saveRDS(mod_dbarts_2,
    ##     file = paste0(data_path, "synthetic/mod_dbarts_2.rds")
    ##   ) 
    ##   ## ----end
    ##   mod_dbarts_2
    ## }),
    ## tar_target(pdp_dbarts_2_, {
    ##   mod_dbarts_2 <- mod_dbarts_2_
    ##   newdata_2 <- site_replacements_newdata_2_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- dbarts_pdp_2
    ##   preds <- predict(mod_dbarts_2, newdata_2, type = "ev") |>
    ##     exp() |> 
    ##     summarise_draws(median, HDInterval::hdi)
    ##   dbarts_2_sum <- newdata_2 |>
    ##     bind_cols(preds) |> 
    ##     mutate(Year = as.numeric(as.character(fYear))) |>
    ##     mutate(type = "dbarts")
    ##   saveRDS(dbarts_2_sum,
    ##     file = paste0(data_path, "synthetic/dbarts_2_sum.rds")
    ##   ) 
    ##   ## ----end
    ## }),
    
    ## ## dbarts 2
    ## tar_target(mod_dbarts2_2_, {
    ##   benthos_fixed_locs_obs_2 <- site_replacements_data_prep_2_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- dbarts2_2
    ##   mod_dbarts2_2 <- bart2(log(cover) ~ fYear + Latitude + Longitude + Reef + Site + Transect,
    ##     data =   benthos_fixed_locs_obs_2,
    ##     keepTrees = TRUE
    ##   )
    ##   saveRDS(mod_dbarts2_2,
    ##     file = paste0(data_path, "synthetic/mod_dbarts2_2.rds")
    ##   ) 
    ##   ## ----end
    ##   mod_dbarts2_2
    ## }),
    ## tar_target(pdp_dbarts2_2_, {
    ##   mod_dbarts2_2 <- mod_dbarts2_2_
    ##   newdata2_2 <- site_replacements_newdata2_2_
    ##   data_path <- site_replacement_global_parameters_$data_path
    ##   ## ---- dbarts2_pdp_2
    ##   preds <- predict(mod_dbarts2_2, newdata2_2, type = "ev") |>
    ##     exp()
    ##   ## summarise_draws(median, HDInterval::hdi)
    ##   dbarts2_2_sum <-
    ##     preds |>
    ##     as_tibble() |>
    ##     mutate(.draw = 1:n()) |>
    ##     pivot_longer(cols = -.draw, names_to = "Variable", values_to = "values") |>
    ##     left_join(
    ##       newdata2_2 |>
    ##         mutate(Variable = paste0("V", 1:n())),
    ##       by = c("Variable")
    ##     ) |>
    ##     group_by(.draw, fYear) |>
    ##     summarise(
    ##       median = median(values)
    ##     ) |>
    ##     ungroup() |>
    ##     group_by(fYear) |>
    ##     tidybayes::summarise_draws(median, HDInterval::hdi) |> 
    ##     mutate(Year = as.numeric(as.character(fYear))) |>
    ##     mutate(type = "dbarts2")
    ##   saveRDS(dbarts2_2_sum,
    ##     file = paste0(data_path, "synthetic/dbarts2_2_sum.rds")
    ##   ) 
    ##   ## ----end
    ## })
  )
  return(targets)
}
