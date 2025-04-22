prepare_data_for_stan <- function(data, yrs) {
  data <- data |> st_drop_geometry()
  if (length(unique(data$cYear)) > 1) {
    X <- model.matrix(~ -1 + cYear, data = data |> droplevels())
  } else {
    X <- model.matrix(~ 1, data = data)
  }
  ## all_years <- 1978:2020
  if (is.null(yrs)) {
    all_years <- as.numeric(as.character(levels(data$cYear)))
  } else {
    all_years <- yrs
  }
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

  ## Xmat <- model.matrix(~ factor(all_years))
  Xmat <- model.matrix(~ -1 + factor(all_years))

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
    Xmat = Xmat,
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
