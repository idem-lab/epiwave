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
#'
#' @param infection_timeseries
#' @param delay_distribution
#' @param proportion_observed
#' @param count_data
#' @param dow_model
#' @param data_id optional name label identifying data type for greta arrays
#'
#' @importFrom greta %*%
#'
#' @return
#' @export
create_observation_model <- function (infection_timeseries,
                                      delay_distribution,
                                      proportion_observed,
                                      count_data,
                                      dow_model = NULL,
                                      data_id = NULL) {

  # add if statements to check that infection_days is long enough to cover
  # the right period

  case_mat <- data_to_matrix(count_data, 'count')
  # delay_mat <- data_to_matrix(delay_distribution)
  prop_mat <- data_to_matrix(proportion_observed, 'prop')

  infection_days <- as.Date(rownames(prop_mat))

  ## check for correct dow arrays
  if (!is.null(dow_model)) {
    dow_correction <- implement_day_of_week(infection_days, dow_model)
    prop_mat <- prop_mat * dow_correction
  }

  n_jurisdictions <- ncol(infection_timeseries)
  n_dates <- nrow(infection_timeseries)

  convolution_matrices <- lapply(1:n_jurisdictions, function(x)
    get_convolution_matrix(delay_distribution[, x],
                           n_dates))

  # compute expected cases of the same length
  # note not all of these dates would have been observed
  expected_cases <- do.call(
    cbind,
    lapply(1:n_jurisdictions, function(x) {
      convolution_matrices[[x]] %*% infection_timeseries[, x] *
        prop_mat[, x]
    }))

  data_idx <- infection_days %in%
    rownames(case_mat)
  expected_cases_idx <- expected_cases[data_idx, ]

  n_days <- nrow(case_mat)

  # negative binomial parameters - need to change from mean and variance
  # specification to size and prob
  sqrt_inv_size <- greta::normal(0, 0.5,
                                 truncation = c(0, Inf),
                                 dim = n_jurisdictions)
  sqrt_inv_size <- greta::sweep(greta::zeros(n_days,
                                             n_jurisdictions),
                                2, sqrt_inv_size,
                                FUN = "+")

  size <- 1 / sqrt(sqrt_inv_size)
  prob <- 1 / (1 + expected_cases_idx / size)

  valid_mat <- case_mat
  valid_mat[is.na(case_mat)] <- FALSE
  valid_mat[!is.na(case_mat)] <- TRUE
  valid_idx <- as.logical(as.numeric(valid_mat))

  case_mat_array <- greta::as_data(
    as.numeric(case_mat)[valid_idx])

  greta::distribution(case_mat_array) <- greta::negative_binomial(
    size[valid_idx],
    prob[valid_idx])

  greta_arrays <- list(
    size,
    prob,
    convolution_matrices
  )

  names(greta_arrays) <- c(
    paste0(data_id, '_size'),
    paste0(data_id, '_prob'),
    paste0(data_id, '_convolution_matrices')
  )

  return(greta_arrays)
}
