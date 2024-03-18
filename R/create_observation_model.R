#' Create observation model
#'
#' @description This function constructs the observation model and defines
#'   likelihood over the observational data (e.g. cases, hospital
#'   admissions...). It first begets a timeseries of expected observations from
#'   a timeseries of new infections, the infection-to-observation delay
#'   distribution(s), and the proportion of infections to be observed.
#'   Specifically, the model uses the delay distribution(s) to convolve the
#'   infection timeseries into the expected observation timeseries, thus
#'   accounting for the probability of observation over time-since-infection. An
#'   additional multiplication by observation proportion is applied to the
#'   convolved timeseries, to adjust for effects such as case ascertainment.
#'   Optionally, a day-of-week effect may be included in the convolution process
#'   to introduce weekly periodicity. The expected observation timeseries is
#'   treated as the mean of a negative binomial distribution from which the data
#'   is observed, thus allowing us to define likelihood over the data, and link
#'   the data to the unknown infection timeseries.
#'
#' @param infection_timeseries greta array of infection timeseries
#' @param proportion_observed long form data, for all dates of infection
#'   timeseries, of expected proportions of infections observed in count data
#' @param dow_model optional module of greta arrays defining day-of-week model
#' @param data_id optional name label identifying data type for greta arrays
#' @param convolution_distribution distribution for convolution, be it delay
#'   distribution for neg binom count data or decay curve for sero/wastewater
#'   data
#' @param obs_data observational data, of any type, in long data format indexed
#'   by date and jurisdiction
#' @param sample_size total sample size in binomial survey data, or
#'   normalisation denominator in the wastewater concentration data. This input
#'   data needs be in long data format and match size with obs_data
#' @param state_pop state level denominator for converting observed infections
#'   into proportion binom data. e.g. total jurisdiction pop for calculating
#'   infection prevalence for seropositivity sample data. This argument needs to
#'   be a vector matching to jurisdiction names in the infection matrix
#' @param likelihood choice for likelihood over data
#' @param noise optional cauchy noise
#'
#' @importFrom greta %*% as_data negative_binomial normal sweep zeros
#'
#' @return greta arrays of observation model
#' @export
create_observation_model <- function (infection_timeseries,
                                      convolution_distribution,
                                      proportion_observed = NULL,
                                      obs_data,
                                      sample_size = NULL,
                                      state_pop = NULL,
                                      dow_model = NULL,
                                      likelihood = c('negbinom',
                                                     'binom',
                                                     'logNormal'),
                                      noise = FALSE,
                                      data_id = NULL) {

  # add if statements to check that infection_days is long enough to cover
  # the right period


  # need to avoid match.arg? can change later
  likelihood <- match.arg(likelihood)

  if (likelihood == 'negbinom') {
    message("negative binomial likelihood over data!")
  }

  if (likelihood == 'binom') {
    message("binomial likelihood over data!")
  }


  if (likelihood == 'logNormal') {
    message("log-Normal likelihood over data!")
  }

  infection_days <- as.Date(rownames(infection_timeseries))
  jurisdictions <- colnames(infection_timeseries)
  n_jurisdictions <- ncol(infection_timeseries)
  n_dates <- nrow(infection_timeseries)

  # construct null prop for wastewater and sero data
  if (is.null(proportion_observed)) {
    proportion_observed <- GPreff::expand_constant_value(
      dates = infection_days,
      jurisdictions = jurisdictions,
      constant_val = 1,
      col_name = 'prop')
  }

  # for wastewater data, need to account for site level effect, so need to
  # modify the data_to_matrix function to accommodate the extra column

  if (likelihood == 'logNormal') {
    obs_mat <- data_with_site_to_matrix(obs_data, 'concentration') # dummy fun
  } else {
    obs_mat <- as_matrix(obs_data, 'count')
  }
  dist_mat <- as_matrix(convolution_distribution, 'delay_fxn')
  prop_mat <- as_matrix(proportion_observed, 'prop')

  # for hospitalisation, multiply the car by the chr to create ihr
  # this is temp code for covid live demo
  if (!is.null(proportion_observed$value)) {
    prop_mat <- prop_mat * proportion_observed$value
  }

  ## add a check for correct dow arrays
  if (!is.null(dow_model)) {
    dow_correction <- implement_day_of_week(infection_days, dow_model)
    prop_mat <- prop_mat * dow_correction
  }

  # add check that ncol(infection_timeseries and below yield same. number of juris)
  # n_jurisdictions <- length(unique(convolution_distribution$jurisdiction))
  #ncol(infection_timeseries)
  # n_dates <- nrow(infection_timeseries)


  convolution_matrices <- lapply(
    unique(delay_distribution$jurisdiction),
    function(x) {
      get_convolution_matrix(convolution_distribution,
                             x,
                             n_dates)
    })

  # compute expected cases of the same length
  expected_obs_list <- lapply(
    1:n_jurisdictions,
    function(x) {
      convolution_matrices[[x]] %*% infection_timeseries[, x] * prop_mat[, x]
    })

  expected_obs <- do.call(
    cbind,
    expected_obs_list
  )

  if (noise == TRUE) {
    # Cauchy noise in obs model
    # add cauchy noise
    noise_sigma <- 0.001
    noise_raw <- cauchy(0, 1, dim = c(nrow(expected_obs),
                                      ncol(expected_obs))
    )
    noise <- noise_raw * noise_sigma
    expected_obs <- expected_obs * exp(noise)
  }

  if (likelihood %in% c('binom')) {
    # divide by denominator to get the correct rates
    expected_obs <- sweep(expected_obs,
                          2,
                          state_pop,
                          FUN = "/")
  }


  data_idx <- infection_days %in% rownames(obs_mat)
  expected_obs_idxed <- expected_obs[data_idx, ]

  # tally the number of days in the observed data
  n_days_idxed <- nrow(obs_mat)

  # construct validity matrix
  valid_mat <- obs_mat
  valid_mat[is.na(obs_mat)] <- FALSE
  valid_mat[!is.na(obs_mat)] <- TRUE
  valid_idx <- as.logical(as.numeric(valid_mat))
  # use it to get the valid array of data
  obs_mat_array <- greta::as_data(
    as.numeric(obs_mat)[valid_idx])

  if (likelihood == 'negbinom') {
    # negative binomial parameters
    sqrt_inv_size <- greta::normal(0, 0.5,
                                   truncation = c(0, Inf),
                                   dim = n_jurisdictions)
    sqrt_inv_size <- greta::sweep(greta::zeros(n_dates,
                                               n_jurisdictions),
                                  2, sqrt_inv_size,
                                  FUN = "+")

    size <- 1 / sqrt(sqrt_inv_size)
    prob <- 1 / (1 + expected_obs / size)

    # truncate size and prob to observed date range
    size_idxed <- size[data_idx, ]
    prob_idxed <- prob[data_idx, ]

    greta::distribution(obs_mat_array) <- greta::negative_binomial(
      size_idxed[valid_idx],
      prob_idxed[valid_idx])

    greta_arrays <- list(
      size,
      prob,
      convolution_matrices
    )

    names(greta_arrays) <- c(
      paste0(data_id, '_negbinom_size'),
      paste0(data_id, '_negbinom_prob'),
      paste0(data_id, '_convolution_matrices')
    )
  }

  if (likelihood == 'binom') {

    # get the binomial sample size from survey in matrix format
    # colname doesn't need to be hard coded
    sample_sizes_mat <- data_to_matrix(sample_size, 'total')

    greta::distribution(obs_mat_array) <- greta::binomial(
      size = sample_sizes_mat[valid_idx],
      prob = expected_obs_idxed[valid_idx]
    )

    greta_arrays <- list(
      expected_obs, # suspect we want to output the full length version?
      convolution_matrices
    )

    names(greta_arrays) <- c(
      paste0(data_id, '_binom_prob'),
      paste0(data_id, '_convolution_matrices')
    )
  }


  if (likelihood == 'logNormal') {

    # time constant nuisance params
    logN_sd <- greta::normal(0,1,truncation = c(.Machine$double.eps,Inf))

    # get the normalisation denominator in matrix format
    sample_sizes_mat <- data_to_matrix(sample_size, 'PMMoV')

    greta::distribution(obs_mat_array) <- greta::lognormal(
      meanlog = expected_obs_idxed[valid_idx]/sample_sizes_mat[valid_idx],
      sdlog = logN_sd # dummy var, time constant
    )

    greta_arrays <- list(
      expected_obs,
      convolution_matrices
    )

    names(greta_arrays) <- c(
      paste0(data_id, '_logNormal_mean'),
      paste0(data_id, '_convolution_matrices')
    )
  }

  return(greta_arrays)
}
