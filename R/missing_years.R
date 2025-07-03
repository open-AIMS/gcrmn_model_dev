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
      ## ----end
      benthos_fixed_locs_obs_3 
    }),
    tar_target(missing_years_newdata_3_, {
      benthos_fixed_locs_obs_3 <- missing_years_data_prep_3_
      ## ---- newdata 3
      newdata_3 <-
        benthos_fixed_locs_obs_3 |>
        tidyr::expand(fYear)
      ## ----end
      newdata_3
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
    tar_target(emmeans_mod_brms_3_, {
      mod_brms_3 <- mod_brms_3_
      newdata_3 <- missing_years_newdata_3_
      data_path <- missing_years_global_parameters_$data_path
      ## ---- brms_3_emmeans
      brms_3_sum <- 
        mod_brms_3 |> emmeans(~fYear, at = newdata_3, type = "response") |>
        as.data.frame() |> 
        rename(median = response, lower = lower.HPD, upper = upper.HPD) |>
        mutate(Year = as.numeric(as.character(fYear))) |>
        mutate(type = "brms")
      saveRDS(brms_3_sum,
        file = paste0(data_path, "synthetic/brms_3_sum.rds")
      ) 
      ## ----end
      brms_3_sum 
    }),
    tar_target(emmeans_mod_brms_3_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      brms_3_sum <- emmeans_mod_brms_3_
      mod_simple_3 <- mod_simple_3_
      ## ---- brms_3_emmeans plot
      brms_3_sum <- readRDS(
        file = paste0(data_path, "synthetic/brms_3_sum.rds")
      )
      g1 <- 
        brms_3_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "brms")) +
        geom_line(data = mod_simple_3,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_3,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <- 
        brms_3_sum |>
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
          fig_path, "R_pdp_mod_brms_3.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
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
    tar_target(pdp_mod_stan_3_, {
      mod_stan_3 <- mod_stan_3_
      newdata_3 <- missing_years_newdata_3_
      data_path <- missing_years_global_parameters_$data_path
      ## ---- stan_3_pdp
      stan_3_sum <-
        mod_stan_3$draws(variables = "cellmeans") |>
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
        bind_cols(newdata_3) |>
        mutate(Year = as.numeric(as.character(fYear)))
      saveRDS(stan_3_sum,
        file = paste0(data_path, "synthetic/stan_3_sum.rds")
      ) 
      ## ----end
      stan_3_sum 
    }),    
    tar_target(pdp_mod_stan_3_plot_, {
      stan_3_sum <- pdp_mod_stan_3_
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_simple_3 <- mod_simple_3_
      ## ---- stan_3_pdp plot
      stan_3_sum <- readRDS(
        file = paste0(data_path, "synthetic/stan_3_sum.rds")
      )
      g1 <- stan_3_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "stan")) +
        geom_line(data = mod_simple_3,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_3,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <- stan_3_sum |>
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
          fig_path, "R_pdp_mod_stan_3.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
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
        cv.folds = 5,
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
    tar_target(pdp_gbm_3_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      gdm_3_sum <- pdp_gbm_3_
      mod_simple_3 <- mod_simple_3_
      ## ---- gbm_pdp_3
      gbm_3_sum <- readRDS(
        file = paste0(data_path, "synthetic/gbm_3_sum.rds")
      )
      g1 <-
        gbm_3_sum |>
        ggplot() +
        geom_line(aes(x = Year, y = median, color = "gbm")) +
        geom_line(data = mod_simple_3,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_3,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <-
        gbm_3_sum |>
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
          fig_path, "R_pdp_mod_gbm_3.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
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
        cv.folds = 5,
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
    tar_target(pdp_gbm_3b_, {
      benthos_fixed_locs_obs_3 <- missing_years_data_prep_3_
      mod_gbm_3b <- mod_gbm_3b_$mod_gbm_3b
      n.trees <- mod_gbm_3b_$n.trees
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      ## ---- gbm_pdp_3b
      newdata_3b <- benthos_fixed_locs_obs_3
      gbm_3b_sum <- newdata_3b |>
        mutate(median = predict(mod_gbm_3b, newdata_3b, n.trees = n.trees, type = "response")) |>
        mutate(Year = as.numeric(as.character(fYear))) |>
        ## group_by(Year, Site, Transect) |>
        ## summarise(median = median(median)) |>
        ## group_by(Year, Site, .add = FALSE) |>
        ## summarise(median = median(median)) |>
        group_by(Year, .add = FALSE) |> 
        summarise(median = median(median)) |> 
        mutate(type = "gbm")
      saveRDS(gbm_3b_sum,
        file = paste0(data_path, "synthetic/gbm_3b_sum.rds")
      ) 
      ## ----end
      gbm_3b_sum
    }),
    tar_target(pdp_gbm_3b_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      gbm_3b_sum <- pdp_gbm_3b_
      mod_simple_3 <- mod_simple_3_
      ## ---- gbm_pdp_3b plot
      gbm_3b_sum <- readRDS(
        file = paste0(data_path, "synthetic/gbm_3b_sum.rds")
      )
      g1 <-
        gbm_3b_sum |>
        ggplot() +
        geom_line(aes(x = Year, y = median, color = "gbm")) +
        geom_line(data = mod_simple_3,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_3,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <-
        gbm_3b_sum |>
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
          fig_path, "R_pdp_mod_gbm_3b.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(infl_gbm_3b_, {
      mod_gbm_3b <- mod_gbm_3b_$mod_gbm_3b
      n.trees <- mod_gbm_3b_$n.trees
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      ## ---- gbm_infl_3b
      ## gbm_3b_infl <- gbm::relative.influence(mod_gbm_3b, n.trees = n.trees, scale = TRUE, sort = TRUE)
      infl <- summary(mod_gbm_3b, n.trees =  n.trees, plot = FALSE)
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
          fig_path, "R_infl_mod_gbm_3b.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(pdp_gbm_3c_, {
      benthos_reefs_sf <- read_all_reefs_data_
      mod_gbm_3b <- mod_gbm_3b_$mod_gbm_3b
      n.trees <- mod_gbm_3b_$n.trees
      data_path <- missing_years_global_parameters_$data_path
      ## ---- gbm_pdp_3c
      newdata_3c <- benthos_reefs_sf |>
        mutate(
          Latitude = st_coordinates(geometry)[,2],
          Longitude = st_coordinates(geometry)[,1],
          fYear = as.factor(Year),
          )
      gbm_3c_sum <- newdata_3c |>
        mutate(median = predict(mod_gbm_3b, newdata_3c, n.trees = n.trees, type = "response")) |>
        mutate(Year = as.numeric(as.character(fYear))) |>
        ## group_by(Year, Site, Transect) |>
        ## summarise(median = median(median)) |>
        ## group_by(Year, Site, .add = FALSE) |>
        ## summarise(median = median(median)) |>
        group_by(Year, .add = FALSE) |> 
        summarise(median = median(median)) |> 
        mutate(type = "gbm")
      saveRDS(gbm_3c_sum,
        file = paste0(data_path, "synthetic/gbm_3c_sum.rds")
      ) 
      ## ----end
      gbm_3c_sum
    }),
    tar_target(pdp_gbm_3c_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      gbm_3c_sum <- pdp_gbm_3c_
      mod_simple_3 <- mod_simple_3_
      ## ---- gbm_pdp_3c plot
      gbm_3c_sum <- readRDS(
        file = paste0(data_path, "synthetic/gbm_3c_sum.rds")
      )
      g1 <-
        gbm_3c_sum |>
        ggplot() +
        geom_line(aes(x = Year, y = median, color = "gbm")) +
        geom_line(data = mod_simple_3,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_3,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <-
        gbm_3c_sum |>
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
          fig_path, "R_pdp_mod_gbm_3c.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(infl_gbm_3c_, {
      mod_gbm_3b <- mod_gbm_3b_$mod_gbm_3b
      n.trees <- mod_gbm_3b_$n.trees
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      ## ---- gbm_infl_3c
      ## gbm_3b_infl <- gbm::relative.influence(mod_gbm_3b, n.trees = n.trees, scale = TRUE, sort = TRUE)
      infl <- summary(mod_gbm_3b, n.trees =  n.trees, plot = FALSE)
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
          fig_path, "R_infl_mod_gbm_3c.png"
        ),
        g,
        width = 8, height = 6, dpi = 72
      )
      ## ----end
    }),
    
    ## dbarts --------------------------------------------------------
    tar_target(mod_dbarts_3_, {
      benthos_fixed_locs_obs_3 <- missing_years_data_prep_3_
      data_path <- missing_years_global_parameters_$data_path
      newdata_3 <- missing_years_newdata_3_
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
        exp() |> 
        summarise_draws(median, HDInterval::hdi)
      saveRDS(preds,
        file = paste0(data_path, "synthetic/mod_dbarts_3_preds.rds")
      ) 
      ## ----end
      list(
        mod_dbarts_3 = mod_dbarts_3,
        preds = preds
      )
    }),
    tar_target(dbarts_pdp_3, {
      data_path <- missing_years_global_parameters_$data_path
      newdata_3 <- missing_years_newdata_3_
      preds <- mod_dbarts_3_$preds
      ## ---- dbarts_pdp_3
      dbarts_3_sum <- newdata_3 |>
        bind_cols(preds) |> 
        mutate(Year = as.numeric(as.character(fYear))) |>
        mutate(type = "dbarts")
      saveRDS(dbarts_3_sum,
        file = paste0(data_path, "synthetic/dbarts_3_sum.rds")
      ) 
      # ----end
      dbarts_3_sum 
    }),
    tar_target(pdp_dbarts_3_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      dbarts_3_sum <- dbarts_pdp_3
      mod_simple_3 <- mod_simple_3_
      ## ---- dbarts_pdp_3 plot
      dbarts_3_sum <- readRDS(
        file = paste0(data_path, "synthetic/dbarts_3_sum.rds")
      )
      g1 <- 
        dbarts_3_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "dbarts")) +
        geom_line(data = mod_simple_3,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_3,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <- 
        dbarts_3_sum |>
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
          fig_path, "R_pdp_mod_dbarts_3.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    
    ## dbarts + covariates ---------------------------------------------
    tar_target(mod_dbarts_3b_, {
      benthos_fixed_locs_obs_3 <- missing_years_data_prep_3_
      benthos_reefs_sf <- read_all_reefs_data_
      data_path <- missing_years_global_parameters_$data_path
      newdata_3 <- missing_years_newdata_3_
      ## ---- dbarts_3b
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
      newdata_3b <- benthos_fixed_locs_obs_3
      preds <- predict(mod_dbarts_3b, newdata_3b, type = "ev") |>
        exp()
      saveRDS(preds,
        file = paste0(data_path, "synthetic/mod_dbarts_3b_preds.rds")
      ) 
      ## ----end
      ## ---- dbarts_pred_3c
      newdata_3c <- benthos_reefs_sf |>
        mutate(
          Latitude = st_coordinates(geometry)[, 2],
          Longitude = st_coordinates(geometry)[, 1],
          fYear = as.factor(Year),
        ) |>
        st_drop_geometry() |>
        group_by(Year, fYear, Reef) |>
        summarise(across(c(Latitude, Longitude, CYC, DHW, OTHER), mean))
      preds_3c <- predict(mod_dbarts_3b, newdata_3c, type = "ev") |>
        exp()
      saveRDS(preds_3c,
        file = paste0(data_path, "synthetic/mod_dbarts_3c_preds.rds")
      ) 
      ## ----end
      list(
        mod_dbarts_3b = mod_dbarts_3b,
        preds = preds,
        preds_3c = preds_3c
      )
    }),
    tar_target(pdp_dbarts_3b_, {
      preds <- mod_dbarts_3b_$preds
      data_path <- missing_years_global_parameters_$data_path
      benthos_fixed_locs_obs_3 <- missing_years_data_prep_3_
      ## ---- dbarts_pdp_3
      dbarts_3b_sum <- preds |> 
        t() |>
        as_tibble(.name_repair = ~ paste0("V", seq_along(.))) |>
        bind_cols(benthos_fixed_locs_obs_3) |>
        pivot_longer(
          cols = matches("^V[0-9]*"),
          names_to = ".draw", values_to = "fit"
        ) |>
        group_by(Year, .draw) |>
        summarise(fit = mean(fit)) |>
        ## dplyr::select(Year, fit, .draw) |> 
        group_by(Year, .add = FALSE) |> 
        summarise_draws(median, HDInterval::hdi)
      saveRDS(dbarts_3b_sum,
        file = paste0(data_path, "synthetic/dbarts_3b_sum.rds")
      ) 
      ## ----end
      dbarts_3b_sum
    }),
    tar_target(pdp_dbarts_3b_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      dbarts_3b_sum <- pdp_dbarts_3b_
      mod_simple_3 <- mod_simple_3_
      ## ---- dbarts_pdp_3 plot
      dbarts_3b_sum <- readRDS(
        file = paste0(data_path, "synthetic/dbarts_3b_sum.rds")
      )
      g1 <- 
        dbarts_3b_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "dbarts")) +
        geom_line(data = mod_simple_3,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_3,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <- 
        dbarts_3b_sum |>
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
          fig_path, "R_pdp_mod_dbarts_3b.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    tar_target(infl_dbarts_3b_plot_, {
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_dbarts_3b <- mod_dbarts_3b_$mod_dbarts_3b
      ## ---- dbarts_infl_3 plot
      rel_inf <- apply(mod_dbarts_3b$varcount, c(2, 3), sum) |>
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
          fig_path, "R_infl_mod_dbarts_3b.png"
        ),
        g,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
    }),
    
    tar_target(pdp_dbarts_3c_, {
      preds_3c <- mod_dbarts_3b_$preds_3c
      benthos_reefs_sf <- read_all_reefs_data_ |>
        mutate(
          Latitude = st_coordinates(geometry)[, 2],
          Longitude = st_coordinates(geometry)[, 1],
          fYear = as.factor(Year),
        ) |>
        st_drop_geometry() |>
        group_by(Year, fYear, Reef) |>
        summarise(across(c(Latitude, Longitude, CYC, DHW, OTHER), mean))
      data_path <- missing_years_global_parameters_$data_path
      ## ---- gbm_pdp_3c
      dbarts_3c_sum <- preds_3c |> 
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
      saveRDS(dbarts_3c_sum,
        file = paste0(data_path, "synthetic/dbarts_3c_sum.rds")
      ) 
      ## ----end
      dbarts_3c_sum
    }),
    tar_target(pdp_dbarts_3c_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      dbarts_3c_sum <- pdp_dbarts_3c_
      mod_simple_3 <- mod_simple_3_
      ## ---- dbarts_pdp_3 plot
      dbarts_3c_sum <- readRDS(
        file = paste0(data_path, "synthetic/dbarts_3c_sum.rds")
      )
      g1 <- 
        dbarts_3c_sum |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "dbarts")) +
        geom_line(data = mod_simple_3,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_3,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <- 
        dbarts_3c_sum |>
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
          fig_path, "R_pdp_mod_dbarts_3c.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
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
    tar_target(emmeans_mod_glmmTMB_4, {
      mod_glmmTMB_4 <- mod_glmmTMB_4_
      newdata_4 <- missing_years_newdata_4_
      data_path <- site_replacement_global_parameters_$data_path
      ## ---- glmmTMB_4_emmeans
      glmmTMB_4_sum <- 
        mod_glmmTMB_4 |> emmeans(~fYear, at = newdata_4, type = "response") |>
        as.data.frame() |> 
        rename(median = response, lower = asymp.LCL, upper = asymp.UCL) |>
        mutate(Year = as.numeric(as.character(fYear))) |>
        mutate(type = "glmmTMB")
      saveRDS(glmmTMB_4_sum,
        file = paste0(data_path, "synthetic/glmmTMB_4_sum.rds")
      ) 
      ## ----end
      glmmTMB_4_sum 
    }),
    tar_target(emmeans_mod_glmmTMB_4_plot_, {
      glmmTMB_4_sum <- emmeans_mod_glmmTMB_4
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      fig_path <- site_replacement_global_parameters_$fig_path
      mod_simple_4 <- mod_simple_4_
      ## ---- glmmTMB_4_emmeans plot
      all_years <- glmmTMB_4_sum |>
        reframe(Year = full_seq(Year, period = 1))
      glmmTMB_4_sum_sam_yrs <- glmmTMB_4_sum |>
        full_join(all_years, by = "Year") 
      gap_years_boundary <- glmmTMB_4_sum |>
        reframe(Year = range(setdiff(full_seq(Year, period = 1), Year)) + c(-1, 1))
      glmmTMB_4_sum_gap_yrs <- glmmTMB_4_sum |>
        right_join(gap_years_boundary, by = "Year") 

      g1 <- glmmTMB_4_sum_sam_yrs |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "glmmTMB")) +
        geom_ribbon(data = glmmTMB_4_sum_gap_yrs,
          aes(x = Year, ymin = lower, ymax = upper), alpha = 0.1) +
        geom_line(data = glmmTMB_4_sum_gap_yrs,
          aes(x = Year, y = median, color = "glmmTMB"), alpha = 0.4) +
        geom_line(data = mod_simple_4,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_4,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        ## geom_line(data = all_sampled_sum,
        ##   aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      g2 <- glmmTMB_4_sum_sam_yrs |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "glmmTMB")) +
        geom_ribbon(data = glmmTMB_4_sum_gap_yrs,
          aes(x = Year, ymin = lower, ymax = upper), alpha = 0.1) +
        geom_line(data = glmmTMB_4_sum_gap_yrs,
          aes(x = Year, y = median, color = "glmmTMB"), alpha = 0.4) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "true mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "true median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_glmmTMB_4.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
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
    tar_target(emmeans_mod_brms_4_, {
      mod_brms_4 <- mod_brms_4_
      newdata_4 <- missing_years_newdata_4_
      data_path <- missing_years_global_parameters_$data_path
      ## ---- brms_4_emmeans
      brms_4_sum <- 
        mod_brms_4 |> emmeans(~fYear, at = newdata_4, type = "response") |>
        as.data.frame() |> 
        rename(median = response, lower = lower.HPD, upper = upper.HPD) |>
        mutate(Year = as.numeric(as.character(fYear))) |>
        mutate(type = "brms")
      saveRDS(brms_4_sum,
        file = paste0(data_path, "synthetic/brms_4_sum.rds")
      ) 
      ## ----end
      brms_4_sum 
    }),
    tar_target(emmeans_mod_brms_4_plot_, {
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      brms_4_sum <- emmeans_mod_brms_4_
      mod_simple_4 <- mod_simple_4_
      ## ---- brms_4_emmeans plot
      brms_4_sum <- readRDS(
        file = paste0(data_path, "synthetic/brms_4_sum.rds")
      )
      all_years <- brms_4_sum |>
        reframe(Year = full_seq(Year, period = 1))
      brms_4_sum_sam_yrs <- brms_4_sum |>
        full_join(all_years, by = "Year") 
      gap_years_boundary <- brms_4_sum |>
        reframe(Year = range(setdiff(full_seq(Year, period = 1), Year)) + c(-1, 1))
      brms_4_sum_gap_yrs <- brms_4_sum |>
        right_join(gap_years_boundary, by = "Year") 

      g1 <- 
        brms_4_sum_sam_yrs |>
        full_join(all_years, by = "Year") |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "brms")) +
        geom_ribbon(data = brms_4_sum_gap_yrs,
          aes(x = Year, ymin = lower, ymax = upper), alpha = 0.1) +
        geom_line(data = brms_4_sum_gap_yrs,
          aes(x = Year, y = median, color = "brms"), alpha = 0.4) +
        geom_line(data = mod_simple_4,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_4,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <- 
        brms_4_sum_sam_yrs |>
        full_join(all_years, by = "Year") |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "brms")) +
        geom_ribbon(data = brms_4_sum_gap_yrs,
          aes(x = Year, ymin = lower, ymax = upper), alpha = 0.1) +
        geom_line(data = brms_4_sum_gap_yrs,
          aes(x = Year, y = median, color = "brms"), alpha = 0.4) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_brms_4.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
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

    ## stan -----------------------------------------------------------
    tar_target(mod_stan_4_, {
      source("model_functions.R")
      benthos_fixed_locs_obs_4 <- missing_years_data_prep_4_
      data_path <- missing_years_global_parameters_$data_path
      site_extra_functions_
      ## ---- stan_pre_4
      benthos_fixed_locs_obs_4 <-
        benthos_fixed_locs_obs_4 |>
        mutate(fYear = factor(fYear, levels = 1:12)) |> 
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
    tar_target(pdp_mod_stan_4_, {
      mod_stan_4 <- mod_stan_4_
      newdata_4 <- missing_years_newdata_4_
      data_path <- missing_years_global_parameters_$data_path
      newdata_4a <- newdata_4 |>
        reframe(fYear = factor(full_seq(as.numeric(as.character(fYear)), period = 1)))
      ## ---- stan_4_pdp
      stan_4_sum <-
        ## mod_stan_4$draws(variables = "cellmeans") |>
        mod_stan_4$draws(variables = "Years") |>
        posterior::as_draws_df() |>
        mutate(across(starts_with("Years"), plogis)) |> 
        posterior::summarise_draws(
          median,
          HDInterval::hdi,
          ~ HDInterval::hdi(., credMass = c(0.9)),
          rhat,
          ess_bulk,
          ess_tail
        ) |>
        rename(lower_90 = V4, upper_90 = V5) |>
        bind_cols(newdata_4a) |>
        mutate(Year = as.numeric(as.character(fYear)))
      saveRDS(stan_4_sum,
        file = paste0(data_path, "synthetic/stan_4_sum.rds")
      ) 
      ## ----end
      stan_4_sum 
    }),    
    tar_target(pdp_mod_stan_4_plot_, {
      stan_4_sum <- pdp_mod_stan_4_
      benthos_reefs_temporal_summary <- read_all_temporal_summary_
      all_sampled_sum <- sampled_simple_raw_means_
      data_path <- missing_years_global_parameters_$data_path
      fig_path <- missing_years_global_parameters_$fig_path
      mod_simple_4 <- mod_simple_4_
      newdata_4 <- missing_years_newdata_4_
      ## ---- stan_4_pdp plot
      stan_4_sum <- readRDS(
        file = paste0(data_path, "synthetic/stan_4_sum.rds")
      )
      all_years <- stan_4_sum |>
        reframe(Year = full_seq(Year, period = 1))
      stan_4_sum_sam_yrs <- newdata_4 |>
        left_join(stan_4_sum, by = "fYear") 
      gap_years <- stan_4_sum_sam_yrs |>
        reframe(Year = setdiff(full_seq(Year, period = 1), Year))
      gap_years <- gap_years |>
        reframe(Year = full_seq(range(Year) + c(-1, 1), period = 1))
      stan_4_sum_gap_yrs <- stan_4_sum |>
        right_join(gap_years, by = "Year") 
      stan_4_sum_sam_yrs <- stan_4_sum_sam_yrs |>
        full_join(all_years, by = "Year")
      g1 <- stan_4_sum_sam_yrs |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "stan")) +
        geom_ribbon(data = stan_4_sum_gap_yrs,
          aes(x = Year, ymin = lower, ymax = upper), alpha = 0.1) +
        geom_line(data = stan_4_sum_gap_yrs,
          aes(x = Year, y = median, color = "stan"), alpha = 0.4) +
        geom_line(data = mod_simple_4,
          aes(x = Year, y = Mean, colour = "simple data mean"), linetype = "dashed") +
        geom_line(data = mod_simple_4,
          aes(x = Year, y = Median, colour = "simple data median"), linetype = "dashed") +
        theme_bw()
      g2 <- stan_4_sum_sam_yrs |>
        ggplot() +
        geom_ribbon(aes(x = Year, ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(aes(x = Year, y = median, color = "stan")) +
        geom_ribbon(data = stan_4_sum_gap_yrs,
          aes(x = Year, ymin = lower, ymax = upper), alpha = 0.1) +
        geom_line(data = stan_4_sum_gap_yrs,
          aes(x = Year, y = median, color = "stan"), alpha = 0.4) +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Mean, colour = "all mean"), linetype = "dashed") +
        geom_line(data = benthos_reefs_temporal_summary,
          aes(x = Year, y = Median, colour = "all median"), linetype = "dashed") +
        geom_line(data = all_sampled_sum,
          aes(x = Year, y = response, colour = type), linetype = "dashed") +
        theme_bw()
      ggsave(
        filename = paste0(
          fig_path, "R_pdp_mod_stan_4.png"
        ),
        g1 + g2,
        width = 12, height = 6, dpi = 72
      )
      ## ----end
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
    })
    
  )
  return(targets)
}
