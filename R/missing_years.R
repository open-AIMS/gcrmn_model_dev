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
      ## ---- newdata 1
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
      ## ---- newdata 1
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
    })


    
  )
  return(targets)
}
