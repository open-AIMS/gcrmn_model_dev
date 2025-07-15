
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
