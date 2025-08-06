
## ---- make_brms_dharma_res_functions
make_brms_dharma_res <- function(brms_model, seed = 10, ...) {
                                        # equivalent to `simulateResiduals(lme4_model, use.u = FALSE)`
                                        # cores are set to 1 just to ensure reproducibility
    options(mc.cores = 1)
    on.exit(options(mc.cores = parallel::detectCores()))
    response <- brms::standata(brms_model)$Y
    ndraws <- nrow(as_draws_df(brms_model))
    manual_preds_brms <- matrix(0, ndraws, nrow(brms_model$data))
    random_terms <- insight::find_random(
                                 brms_model, split_nested = TRUE, flatten = TRUE
                             )
                                        # for this to have a similar output to `glmmTMB`'s default, we need to
                                        #   create new levels in the hierarchical variables, so then we can
                                        #   use `allow_new_levels = TRUE` and `sample_new_levels = "gaussian"` in
                                        #   `brms::posterior_epred`. This is equivalent to
                                        #   `simulateResiduals(lme4_model, use.u = FALSE)`. See details in
                                        #   `lme4:::simulate.merMod` and `glmmTMB:::simulate.glmmTMB`
                                        ## random_terms <- unlist(str_split(random_terms, ".\\+."))
  random_terms <- str_subset(random_terms, "\\+", negate = TRUE)
    new_data <- brms_model$data |>
        dplyr::mutate(across(
                   all_of(random_terms), \(x)paste0("NEW_", x) |> as.factor()
               ))
    set.seed(seed)
    brms_sims <- brms::posterior_predict(
                           brms_model, re_formula = NULL, newdata = new_data,
                           allow_new_levels = TRUE, sample_new_levels = "gaussian"
                       ) |>
        t()
    fitted_median_brms <- apply(brms_sims, 1, median)
    ## fitted_median_brms <- apply(
    ##     t(brms::posterior_epred(brms_model, ndraws = ndraws, re.form = NA)),
    ##     1,
    ##     mean)
    DHARMa::createDHARMa(
                simulatedResponse = brms_sims,
                observedResponse = response,
                fittedPredictedResponse = fitted_median_brms,
                ...
            )
}
## ----end


source("model_functions.R")
helper_functions <- function() {
  targets <- list(
    tar_target(DHARMa_glmmTMB_, {
      ## ---- DHARMa_glmmTMB_function
      DHARMa_glmmTMB <- function(mod, data_path, fig_path, filename = "mod_glmmTMB_12") {
        resids <- mod |> 
          simulateResiduals(n = 1000)
        saveRDS(resids,
          file = paste0(data_path, paste0("synthetic/", filename, "_dharma.rds"))
        )
        g <- wrap_elements(~testUniformity(resids)) +
          wrap_elements(~plotResiduals(resids)) +
          wrap_elements(~testDispersion(resids))
        ggsave(
          filename = paste0(
            fig_path, "R_dharma_", filename, ".png"
          ),
          g,
          width = 10, height = 4, dpi = 72
        )
      }
      ## ----end
      DHARMa_glmmTMB
    }),

    tar_target(pred_glmmTMB_, {
      ## ---- pred_glmmTMB_function
      pred_glmmTMB <- function(mod,
                               type = 1,
                               model_type = "",
                               newdata = newdata,
                               true_sum = true_sum,
                               data_path, fig_path,
                               filename = "mod_glmmTMB_12") {
        model_type <- ifelse(model_type == "covariates", "b", model_type)
        if (model_type == "") {
          
          pred_glmmTMB <- 
            mod |> emmeans(~fYear, at = newdata, type = "response") |>
            as.data.frame() |> 
            rename(median = response, lower = asymp.LCL, upper = asymp.UCL) |>
            mutate(Year = as.numeric(as.character(fYear))) |>
            mutate(type = "glmmTMB")
          saveRDS(pred_glmmTMB,
            file = paste0(data_path, paste0("synthetic/", filename, model_type, "_", type, ".rds"))
          ) 
        } else if (model_type == "b") {
          newdata <- newdata |>
            mutate(fYear = factor(Year, levels = unique(Year)),
              Site = NA, Transect = NA)
          if ("geometry" %in% names(newdata) & !"Latitude" %in% names(newdata)) {
            newdata <- newdata |>
              mutate(
                Longitude = st_coordinates(st_centroid(newdata))[, 1],
                Latitude = st_coordinates(st_centroid(newdata))[, 2],
                )
          }
          pred <- predict(mod, newdata = newdata,
            re.form =  ~ 0,
            allow.new.levels = TRUE
            ## se.fit = TRUE
          )

          pred_glmmTMB <- 
            newdata |>
            ungroup() |> 
            mutate(
              fit = plogis(pred)
            ) |>
            group_by(Year) |>
            summarise(median = median(fit, na.rm = TRUE),
              lower = HDInterval::hdi(fit)[1],
              upper = HDInterval::hdi(fit)[2]
              )
          saveRDS(pred_glmmTMB,
            file = paste0(data_path, paste0("synthetic/", filename, model_type, "_", type, ".rds"))
          ) 
          
        }

        g1 <-
          pred_glmmTMB |>
          ggplot() +
          geom_ribbon(aes(x = Year, ymin = lower, ymax = upper, fill = "glmmTMB", colour = "glmmTMB"),
            alpha = 0.8) +
          geom_line(aes(x = Year, y = median, color = "glmmTMB"),
            linewidth = 2) +
          geom_line(data = true_sum,
            aes(x = Year, y = Mean, colour = "simple data mean"),
            linetype = "dashed",
            linewidth = 1) +
          geom_line(data = true_sum,
            aes(x = Year, y = Median, colour = "simple data median"),
            linetype = "dashed",
            linewidth = 1) +
          scale_color_manual("",
            values = c("white", "green", "blue"),
            breaks = c("glmmTMB", "simple data mean", "simple data median")
            ) +
          scale_fill_manual("",
            values = c("orange", NA, NA),
            breaks = c("glmmTMB", "simple data mean", "simple data median")
            ) +
          scale_y_continuous("Coral cover (%)", labels = function(x) x * 100) +
          guides(
            fill = "none", #guide_legend(override.aes = list(color = "orange")),
            color = guide_legend(override.aes = list(fill = c("orange", NA, NA)))
          ) +
          theme_classic() 
        
        ggsave(
          filename = paste0(
            fig_path, "R_pdp_", type, "_", filename, model_type, ".png"
          ),
          g1,
          width = 6, height = 4, dpi = 72
        )
      }
      ## ----end
      pred_glmmTMB
    }),
    
    tar_target(mse_glmmTMB_, {
      calc_mse <- calc_mse_
      ## ---- mse_glmmTMB_function
      mse_glmmTMB <- function(mod, newdata = newdata, type, model_type = "covariates") {
        newdata <- newdata |>
          mutate(fYear = factor(Year, levels = unique(Year)),
            Site = NA, Transect = NA)
        if ("geometry" %in% names(newdata) & !"Latitude" %in% names(newdata)) {
          newdata <- newdata |>
            mutate(
              Longitude = st_coordinates(st_centroid(newdata))[, 1],
              Latitude = st_coordinates(st_centroid(newdata))[, 2],
              )
        }
        pred <- predict(mod, newdata = newdata,
          re.form =  ~ 0,
          allow.new.levels = TRUE
        )
        newdata <- newdata |>
          ungroup() |> 
          mutate(fit = plogis(pred))
        mse_glmmTMB <- calc_mse(newdata) |>
          mutate(model = "glmmTMB", type = type, model_type = model_type)
      }
      ## ----end
      mse_glmmTMB
    }),

    tar_target(pred_brms_, {
      ## ---- pred_brms_function
      pred_brms <- function(mod,
                               type = 1,
                               model_type = "",
                               newdata = newdata,
                               true_sum = true_sum,
                               data_path, fig_path,
                               filename = "mod_brms_12") {
        model_type <- ifelse(model_type == "covariates", "b", model_type)
        if (model_type == "") {
          
          pred_brms <- 
            mod |> emmeans(~fYear, at = newdata, type = "response") |>
            as.data.frame() |> 
            rename(median = response, lower = lower.HPD, upper = upper.HPD) |>
            mutate(Year = as.numeric(as.character(fYear))) |>
            mutate(type = "brms")
          saveRDS(pred_brms,
            file = paste0(data_path, paste0("synthetic/", filename, model_type, "_", type, ".rds"))
          ) 
        } else if (model_type == "b") {
          newdata <- newdata |>
            mutate(fYear = factor(Year, levels = unique(Year)))
          if ("geometry" %in% names(newdata) & !"Latitude" %in% names(newdata)) {
            newdata <- newdata |>
              mutate(
                Longitude = st_coordinates(st_centroid(newdata))[, 1],
                Latitude = st_coordinates(st_centroid(newdata))[, 2],
                )
          }
          newdata <- newdata |>
            dplyr::select(Reef, Longitude, Latitude, Year, fYear, CYC, DHW, OTHER)
          pred <- add_epred_draws(mod, newdata = newdata,
            re_formula =  ~ 0,
          )

          pred_brms <- 
            pred |> 
            ungroup() |> 
            group_by(Year, .draw) |>
            summarise(fit = mean(.epred)) |>
            ungroup() |>
            group_by(Year) |>
            summarise_draws(median, HDInterval::hdi)
            
          saveRDS(pred_brms,
            file = paste0(data_path, paste0("synthetic/", filename, model_type, "_", type, ".rds"))
          ) 
          
        }

        g1 <-
          pred_brms |>
          ggplot() +
          geom_ribbon(aes(x = Year, ymin = lower, ymax = upper, fill = "brms", colour = "brms"),
            alpha = 0.8) +
          geom_line(aes(x = Year, y = median, color = "brms"),
            linewidth = 2) +
          geom_line(data = true_sum,
            aes(x = Year, y = Mean, colour = "simple data mean"),
            linetype = "dashed",
            linewidth = 1) +
          geom_line(data = true_sum,
            aes(x = Year, y = Median, colour = "simple data median"),
            linetype = "dashed",
            linewidth = 1) +
          scale_color_manual("",
            values = c("white", "green", "blue"),
            breaks = c("brms", "simple data mean", "simple data median")
            ) +
          scale_fill_manual("",
            values = c("orange", NA, NA),
            breaks = c("brms", "simple data mean", "simple data median")
            ) +
          scale_y_continuous("Coral cover (%)", labels = function(x) x * 100) +
          guides(
            fill = "none", #guide_legend(override.aes = list(color = "orange")),
            color = guide_legend(override.aes = list(fill = c("orange", NA, NA)))
          ) +
          theme_classic() 
        
        ggsave(
          filename = paste0(
            fig_path, "R_pdp_", type, "_", filename, model_type, ".png"
          ),
          g1,
          width = 6, height = 4, dpi = 72
        )
      }
      ## ----end
      pred_brms
    }),

    tar_target(mse_brms_, {
      calc_mse <- calc_mse_
      ## ---- mse_brms_function
      mse_brms <- function(mod, newdata = newdata, type, model_type = "covariates") {
        newdata <- newdata |>
          mutate(fYear = factor(Year, levels = unique(Year)),
            Site = NA, Transect = NA)
        if ("geometry" %in% names(newdata) & !"Latitude" %in% names(newdata)) {
          newdata <- newdata |>
            mutate(
              Longitude = st_coordinates(st_centroid(newdata))[, 1],
              Latitude = st_coordinates(st_centroid(newdata))[, 2],
              )
        }

        pred <- add_epred_draws(mod, newdata = newdata,
          re_formula =  ~ 0,
          )
        newdata <- pred |>
          ungroup() |> 
          mutate(fit = .epred)
        mse_brms <- calc_mse(newdata) |>
          mutate(model = "brms", type = type, model_type = model_type)
      }
      ## ----end
      mse_brms
    }),

    tar_target(pred_stan_, {
      ## ---- pred_stan_function
      pred_stan <- function(mod,
                               type = 1,
                               model_type = "",
                               newdata = newdata,
                               true_sum = true_sum,
                               data_path, fig_path,
                               filename = "mod_stan_12") {
        model_type <- ifelse(model_type == "covariates", "b", model_type)
        if (model_type == "") {
          
          stan_sum <-
            mod$draws(variables = "cellmeans") |>
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
            bind_cols(newdata) |>
            mutate(Year = as.numeric(as.character(fYear)))


          pred_stan <- 
            stan_sum |> 
            mutate(type = "stan")
          saveRDS(pred_stan,
            file = paste0(data_path, paste0("synthetic/", filename, model_type, "_", type, ".rds"))
          ) 
        } else if (model_type == "b") {
          
        }

        g1 <-
          pred_stan |>
          ggplot() +
          geom_ribbon(aes(x = Year, ymin = lower, ymax = upper, fill = "stan", colour = "stan"),
            alpha = 0.8) +
          geom_line(aes(x = Year, y = median, color = "stan"),
            linewidth = 2) +
          geom_line(data = true_sum,
            aes(x = Year, y = Mean, colour = "simple data mean"),
            linetype = "dashed",
            linewidth = 1) +
          geom_line(data = true_sum,
            aes(x = Year, y = Median, colour = "simple data median"),
            linetype = "dashed",
            linewidth = 1) +
          scale_color_manual("",
            values = c("white", "green", "blue"),
            breaks = c("stan", "simple data mean", "simple data median")
            ) +
          scale_fill_manual("",
            values = c("orange", NA, NA),
            breaks = c("stan", "simple data mean", "simple data median")
            ) +
          scale_y_continuous("Coral cover (%)", labels = function(x) x * 100) +
          guides(
            fill = "none", #guide_legend(override.aes = list(color = "orange")),
            color = guide_legend(override.aes = list(fill = c("orange", NA, NA)))
          ) +
          theme_classic() 
        
        ggsave(
          filename = paste0(
            fig_path, "R_pdp_", type, "_", filename, model_type, ".png"
          ),
          g1,
          width = 6, height = 4, dpi = 72
        )
      }
      ## ----end
      pred_stan
    }),

    tar_target(mse_stan_, {
      calc_mse <- calc_mse_
      ## ---- mse_stan_function
      mse_stan <- function(mod, newdata = newdata, type, model_type = "covariates") {
        newdata <- newdata |>
          ungroup() |> 
          mutate(fYear = factor(Year, levels = unique(Year)),
            Site = NA, Transect = NA)
        if ("geometry" %in% names(newdata) & !"Latitude" %in% names(newdata)) {
          newdata <- newdata |>
            mutate(
              Longitude = st_coordinates(st_centroid(newdata))[, 1],
              Latitude = st_coordinates(st_centroid(newdata))[, 2],
              )
        }
        Xmat <- model.matrix(~ -1 + fYear, data = newdata)
        stan_sum <-
          mod$draws(variables = "beta") |>
          posterior::as_draws_df() |>
          dplyr::select(-.iteration, -.chain) |>
          group_by(.draw) |>
          nest() |>
          mutate(pred = map(.x = data,
            .f = ~ {
              fit <-
                .x |> 
                as.matrix() |>
                as.vector() %*% t(Xmat) |>
                as.vector() |>
                plogis()
              newdata |> mutate(fit = fit) 
            })) |>
          dplyr::select(-data) |>
          unnest(pred)
        
        newdata <- stan_sum |>
          ungroup() 
        mse_stan <- calc_mse(newdata) |>
          mutate(model = "stan", type = type, model_type = model_type)
      }
      ## ----end
      mse_stan
    }),

    tar_target(pred_gbm_, {
      ## ---- pred_gbm_function
      pred_gbm <- function(mod,
                           n.trees,
                               type = 1,
                               model_type = "",
                               newdata = newdata,
                               true_sum = true_sum,
                               data_path, fig_path,
                               filename = "mod_gbm_12") {
        model_type <- ifelse(model_type == "covariates", "b", model_type)
        if (model_type == "") {
          pred_gbm <- newdata |>
            mutate(median = predict(mod, newdata, n.trees = n.trees, type = "response")) |> 
            mutate(Year = as.numeric(as.character(fYear))) |>
            mutate(type = "gbm")

          saveRDS(pred_gbm,
            file = paste0(data_path, paste0("synthetic/", filename, model_type, "_", type, ".rds"))
          ) 
        } else if (model_type == "b") {
          newdata <- newdata |>
            mutate(fYear = factor(Year, levels = unique(Year))) |>
            ungroup()  
          pred_gbm <- newdata |>
            mutate(median = predict(mod, newdata, n.trees = n.trees, type = "response")) |> 
            ## mutate(Year = as.numeric(as.character(fYear))) |>
            group_by(Year) |>
            summarise(median = median(median)) |> 
            mutate(type = "gbm") |>
            ungroup()

          saveRDS(pred_gbm,
            file = paste0(data_path, paste0("synthetic/", filename, model_type, "_", type, ".rds"))
          ) 
        }

        g1 <-
          pred_gbm |>
          ggplot() +
          ## geom_ribbon(aes(x = Year, ymin = lower, ymax = upper, fill = "gbm", colour = "gbm"),
          ##   alpha = 0.8) +
          geom_line(aes(x = Year, y = median, color = "gbm"),
            linewidth = 2) +
          geom_line(data = true_sum,
            aes(x = Year, y = Mean, colour = "simple data mean"),
            linetype = "dashed",
            linewidth = 1) +
          geom_line(data = true_sum,
            aes(x = Year, y = Median, colour = "simple data median"),
            linetype = "dashed",
            linewidth = 1) +
          scale_color_manual("",
            values = c("orange", "green", "blue"),
            breaks = c("gbm", "simple data mean", "simple data median")
            ) +
          ## scale_fill_manual("",
          ##   values = c("orange", NA, NA),
          ##   breaks = c("gbm", "simple data mean", "simple data median")
          ##   ) +
          scale_y_continuous("Coral cover (%)", labels = function(x) x * 100) +
          ## guides(
          ##   fill = "none", #guide_legend(override.aes = list(color = "orange")),
          ##   color = guide_legend(override.aes = list(fill = c("orange", NA, NA)))
          ## ) +
          theme_classic() 
        
        ggsave(
          filename = paste0(
            fig_path, "R_pdp_", type, "_", filename, model_type, ".png"
          ),
          g1,
          width = 6, height = 4, dpi = 72
        )
      }
      ## ----end
      pred_gbm
    }),
    
    tar_target(mse_gbm_, {
      calc_mse <- calc_mse_
      ## ---- mse_gbm_function
      mse_gbm <- function(mod, newdata = newdata, n.trees, type, model_type = "covariates") {
        newdata <- newdata |>
          ungroup() |> 
          mutate(fYear = factor(Year, levels = unique(Year)),
            Site = NA, Transect = NA)
        if ("geometry" %in% names(newdata) & !"Latitude" %in% names(newdata)) {
          newdata <- newdata |>
            mutate(
              Longitude = st_coordinates(st_centroid(newdata))[, 1],
              Latitude = st_coordinates(st_centroid(newdata))[, 2],
              )
        }

        newdata <- newdata |>
          mutate(fit = predict(mod, newdata, n.trees = n.trees, type = "response")) |> 
          mutate(type = "gbm") |>
          ungroup()
        mse_gbm <- calc_mse(newdata) |>
          mutate(model = "gbm", type = type, model_type = model_type)
      }
      ## ----end
      mse_gbm
    }),

    tar_target(pred_dbarts_, {
      ## ---- pred_dbarts_function
      pred_dbarts <- function(preds,
                              type = 1,
                               model_type = "",
                               newdata = newdata,
                               true_sum = true_sum,
                               data_path, fig_path,
                               filename = "mod_dbarts_12") {
        model_type <- ifelse(model_type == "covariates", "b", model_type)
        if (model_type == "") {
          pred_dbarts <- newdata |>
            bind_cols(preds) |> 
            mutate(Year = as.numeric(as.character(fYear))) |>
            mutate(type = "dbarts")
          saveRDS(pred_dbarts,
            file = paste0(data_path, paste0("synthetic/", filename, model_type, "_", type, ".rds"))
          ) 
        } else if (model_type == "b") {
          full_preds <- 
            t(preds) |>
            as.data.frame() |>
            mutate(i = 1:n()) |>
            pivot_longer(cols = -i, names_to = ".draw") |>
            mutate(.draw = as.numeric(str_remove(.draw, "V"))) |>
            left_join(newdata |>
                        ungroup() |> 
                        mutate(i = 1:n()),
              by = "i") |>
            ungroup() 

          pred_dbarts <-
            full_preds |> 
            group_by(Year, .draw) |>
            summarise(median = median(value)) |>
            ungroup() |>
            group_by(Year) |>
            summarise_draws(median, HDInterval::hdi) |> 
            dplyr::select(-variable) |> 
            mutate(type = "dbarts")

          saveRDS(pred_dbarts,
            file = paste0(data_path, paste0("synthetic/", filename, model_type, "_", type, ".rds"))
          ) 
        }

        g1 <-
          pred_dbarts |>
          ggplot() +
          geom_ribbon(aes(x = Year, ymin = lower, ymax = upper, fill = "dbarts", colour = "dbarts"),
            alpha = 0.8) +
          geom_line(aes(x = Year, y = median, color = "dbarts"),
            linewidth = 2) +
          geom_line(data = true_sum,
            aes(x = Year, y = Mean, colour = "simple data mean"),
            linetype = "dashed",
            linewidth = 1) +
          geom_line(data = true_sum,
            aes(x = Year, y = Median, colour = "simple data median"),
            linetype = "dashed",
            linewidth = 1) +
          scale_color_manual("",
            values = c("white", "green", "blue"),
            breaks = c("dbarts", "simple data mean", "simple data median")
            ) +
          scale_fill_manual("",
            values = c("orange", NA, NA),
            breaks = c("dbarts", "simple data mean", "simple data median")
            ) +
          scale_y_continuous("Coral cover (%)", labels = function(x) x * 100) +
          guides(
            fill = "none", #guide_legend(override.aes = list(color = "orange")),
            color = guide_legend(override.aes = list(fill = c("orange", NA, NA)))
          ) +
          theme_classic() 
        
        ggsave(
          filename = paste0(
            fig_path, "R_pdp_", type, "_", filename, model_type, ".png"
          ),
          g1,
          width = 6, height = 4, dpi = 72
        )
      }
      ## ----end
      pred_dbarts
    }),

    tar_target(mse_dbarts_, {
      calc_mse <- calc_mse_
      ## ---- mse_dbarts_function
      mse_dbarts <- function(preds, newdata = newdata, type, model_type = "covariates") {

        newdata <- 
          t(preds) |>
          as.data.frame() |>
          mutate(i = 1:n()) |>
          pivot_longer(cols = -i, names_to = ".draw", values_to = "fit") |>
          mutate(.draw = as.numeric(str_remove(.draw, "V"))) |>
          left_join(newdata |>
                      ungroup() |> 
                      mutate(i = 1:n()),
            by = "i") |>
          ungroup() 
        mse_dbarts <- calc_mse(newdata) |>
          mutate(model = "dbarts", type = type, model_type = model_type)
      }
      ## ----end
      mse_dbarts
    }),

    tar_target(calc_mse_, {
      ## ---- calc_mse_function
      calc_mse <- function(dat) {
        dat |>
          mutate(
            se_r = (fit - HCC)^2,
            acc = exp(log(fit) - log(HCC)),
            acc = abs(acc - 1)
          ) |>
          summarise(
            mse_mean = mean(se_r),
            mse_median = median(se_r),
            acc_mean = mean(acc),
            acc_median = median(acc),
            lower = HDInterval::hdi(acc)[1],
            upper = HDInterval::hdi(acc)[2]
          )
      }
      ## ----end
      calc_mse
    }
    )

    
  )
  return(targets)
}


## ---- GCRMN_plotCellMeans_function
GCMRN_plotCellMeans <- function(cellMeans, pointsize=1, linesize=0.75,
                                ytitle='', title='', resp, clr='#0DC8F6',
                                baseSize=11, baseFam='Open Sans', guidelines=TRUE, ylims=c(NA,NA),
                                watermark=FALSE,
                                decimal.mark='.') {
  if (watermark) {
    library(grid)
    library(png)
    gcrmn_logo <- png::readPNG('../docs/resources/GCRMN_logo.png')
  } 
  cm = cellMeans %>%
    mutate(Year=as.numeric(Year)) %>%
    filter(Year<2020) %>%
    mutate(rleid = with(rle(Data), rep(seq_along(lengths), lengths)),
      group=as.integer(rleid)) #%>%
  #group_by(rleid)

  cm1 = cm %>% ungroup %>% #group_by(ECOREGION) %>%
    mutate(d=3) %>% uncount(d, .id='A') %>%
    mutate_at(vars(Year, value, .lower_0.8, .upper_0.8, .lower_0.95, .upper_0.95),
      function(x=.) ifelse(.$A==1,(x+lag(x))/2,
        ifelse(.$A==3, (x+lead(x))/2, x))) %>%
    group_by_at(group_vars(cm)) %>%
    filter(row_number()!=1, row_number()!=n())

  g1 <- cm1  %>% ungroup %>%
    ggplot()
  if(watermark) {
    g1 <- g1 +
      annotation_custom(rasterGrob(gcrmn_logo,  width=unit(1, 'npc'), height=unit(1,  'npc')), -Inf,  Inf, -Inf, Inf)
  }
  g1 <- g1 +
    geom_ribbon(aes(ymin=.lower_0.95, ymax=.upper_0.95, x=Year,
      fill=as.factor(Data), group=group),show.legend=FALSE,
      alpha=0.2) +
    geom_ribbon(aes(ymin=.lower_0.8, ymax=.upper_0.8, x=Year, fill=as.factor(Data), group=group),show.legend=FALSE, alpha=0.2) +
    geom_line(aes(y=value, x=Year, color=as.factor(Data), group=group), show.legend=FALSE, size=linesize) +
    geom_point(data=cm, aes(y=value, x=Year, color=as.factor(Data)), show.legend=FALSE, size=pointsize) +
    scale_fill_manual(breaks=c('0','1'), values=c('grey',clr)) +
    scale_color_manual(breaks=c('0','1'), values=c('grey',clr)) +
    ## scale_y_continuous(ytitle, labels=function(x) {x*100}, expand=c(0,0), breaks=scales::breaks_width(0.05), limits=ylims) +
    scale_y_continuous(ytitle,
      labels=scales::label_number(accuracy = 1, scale=100, decimal.mark=decimal.mark),
      expand=c(0,0), breaks=scales::breaks_width(0.05), limits=ylims) +
    scale_x_continuous(breaks=scales::breaks_width(5)) +
    theme_bw(baseSize, base_family = baseFam) +
    theme(axis.title.x = element_blank(),
      axis.title.y = element_text(size=rel(1.5), margin=margin(r=1, unit='lines'))
    ) +
    ggtitle(title)
  if (!guidelines) g1 <- g1 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
    panel.border=element_blank(),
    axis.line.x.bottom=element_line(lineend='square'), axis.line.y.left=element_line(lineend='square')) +
      scale_y_continuous(ytitle,
        labels=scales::label_number(accuracy = 1, scale=100, decimal.mark=decimal.mark), expand=c(0,0), limits=ylims) #+
  ## scale_y_continuous(ytitle, labels=function(x) {x*100}, expand=c(0,0), limits=ylims) #+
  if (resp %in% c('ca','cma')) {
    log_lims <- c(min(cm1$.lower_0.95), max(cm1$.upper_0.95))
    ## g1 = g1 + scale_y_continuous(ytitle,trans=scales::log2_trans(),  breaks=2^(c(0, log(c(0.1, 0.2, 0.5, 2, 5, 10, 20, 50, 100), 2))), labels=function(x) {x}, expand=c(0,0), limits=ylims) #g1 = g1 + scale_y_continuous(ytitle, labels=function(x) 2^x, expand=c(0,0))
    g1 = g1 + scale_y_continuous(ytitle,trans=scales::log2_trans(),  breaks=ceiling(2^pretty(log(log_lims,2))*2)/2, labels=scales::label_number(decimal.mark=decimal.mark, accuracy = 0.1, trim=TRUE), expand=c(0,0), limits=ylims) #g1 = g1 + scale_y_continuous(ytitle, labels=function(x) 2^x, expand=c(0,0))
    ## g1 = g1 + scale_y_continuous(ytitle,trans=scales::log2_trans(),  breaks=ceiling(2^pretty(log(log_lims,2))*2)/2, labels=function(x) {x}, expand=c(0,0), limits=ylims) #g1 = g1 + scale_y_continuous(ytitle, labels=function(x) 2^x, expand=c(0,0))
  }
  g1
}
## ----end


width_height_ratio <- function(object) {
  bb <- object |> st_bbox()
  ht <- bb[4] - bb[2]
  wd <- bb[3] - bb[1]
  wd / ht
}


## ---- prepare_data_for_stan_function
prepare_data_for_stan <- function(data) {
  data <- data |> st_drop_geometry()
  if (length(unique(data$cYear)) > 1) {
    X <- model.matrix(~ -1 + cYear, data = data |> droplevels())
  } else {
    X <- model.matrix(~ 1, data = data)
  }
  ## all_years <- 1978:2020
  all_years <- as.numeric(as.character(levels(data$cYear)))
  data_years <- which(all_years %in% (data |> pull(Year) |> unique() |> sort()))
  no_data_years <- which(!(all_years %in% (data |> pull(Year) |> unique() |> sort())))
  gap_years <- no_data_years[no_data_years > data_years[1] & no_data_years < data_years[length(data_years)]]
  prior_years <- no_data_years[no_data_years < data_years[1]]
  init_year <- data |>
    filter(Year == min(Year)) |>
    pull(Year) |>
    unique()
  init_year <- which(init_year == all_years)
  post_years <- which(all_years > init_year)
  init_cover <- binomial()$linkfun(data |> filter(Year == min(Year)) |>
                                     droplevels() |>
                                     summarise(avCover = mean(Cover)) |> as.numeric())
  non_init_year <- (1:length(all_years))[-init_year]
  between_years <- gap_years[gap_years < max(data_years)]
  after_years <- post_years[post_years > max(data_years)]
  N <- nrow(data)

  ## Get weights
  grid_wts <- data |>
    group_by(grid_id) |>
    summarise(area = unique(sum)) |>
    ungroup() |>
    mutate(wt = area / sum(area))
  stan_data <- list(
    N = N,
    Y = data$Cover,
    K = ncol(X),
    X = X,
    Z_1_1 = rep(0, N),
    Z_2_1 = rep(0, N),
    Z_3_1 = rep(0, N),
    J_1 = as.numeric(factor(data$grid_id)),    ## grid_id
    J_2 = as.numeric(factor(data$cSite)),      ## cSite
    J_3 = as.numeric(factor(data$cReplicate)), ## cReplicate
    N_1 = length(unique(factor(data$grid_id))),
    M_1 = 1,
    NC_1 =  0,
    wt_1 = grid_wts$wt,
    N_2 = length(unique(factor(data$cSite))),
    M_2 = 1,
    NC_2 =  0,
    N_3 = length(unique(factor(data$cReplicate))),
    M_3 = 1,
    NC_3 =  0,
    n_all_years = length(all_years),
    all_years = all_years,
    n_prior_years = length(prior_years),
    prior_years = prior_years,
    n_post_years = length(post_years),
    post_years = post_years,
    init_year = init_year,
    init_cover = init_cover,
    n_gap_years = length(gap_years),
    gap_years = gap_years,
    n_after_years = length(after_years),
    after_years = after_years,
    n_data_years = length(data_years),
    data_years = data_years
  )
  stan_data
}
## ----end

## ---- stan_draws_function
stan_draws <- function(fit) {
  nm <- str_replace(fit, "fit", "draws")
  fit <- readRDS(fit)
  draws <- fit$draws(variables = "Years", format = "df") |>
    mutate(across(starts_with("Year"), plogis))  
  saveRDS(draws, file = nm)
  nm
}
## ----end

## ---- calculate_cellmeans_function
calculate_cellmeans <- function(draws, data, stan_data) {
  nm <- str_replace(draws, "draws", "cellmeans")
  draws <- readRDS(draws)
  all_years <- stan_data$all_years
  cellmeans <-
    draws |> 
    summarise_draws(median, HDInterval::hdi, ~ HDInterval::hdi(., credMass = 0.8)) |>
    bind_cols(years = all_years) |>
    rename(
      value = median,
      .lower_0.8 = V4,
      .lower_0.95 = lower,
      .upper_0.8 = V5,
      .upper_0.95 = upper,
      Year = years
    ) |> 
    mutate(
      ECOREGION = unique(data$ECOREGION),
      GCRMN_region = unique(data$GCRMN_region),
      GCRMN_subregion = unique(data$GCRMN_subregion),
      Data = ifelse(Year %in% Year[stan_data$data_years], 1, 0)
    )
  saveRDS(cellmeans, file = nm)
  nm
}
## ----end
