#' Create observation model
#'
#' @description This function constructs the observation model and defines
#'  likelihood over the observational data (e.g. cases, hospital
#'  admissions...). It first begets a timeseries of expected observations from
#'  a timeseries of new infections, the infection-to-observation delay
#'  distribution(s), and the proportion of infections to be observed.
#'  Specifically, the model uses the delay distribution(s) to convolve the
#'  infection timeseries into the expected observation timeseries, thus
#'  accounting for the probability of observation over time-since-infection. An
#'  additional multiplication by observation proportion is applied to the
#'  convolved timeseries, to adjust for effects such as case ascertainment.
#'  Optionally, a day-of-week effect may be included in the convolution process
#'  to introduce weekly periodicity. The expected observation timeseries is
#'  treated as the mean of a negative binomial distribution from which the data
#'  is observed, thus allowing us to define likelihood over the data, and link
#'  the data to the unknown infection timeseries.
#'
#' @param infection_timeseries greta array of infection timeseries
#' @param delay_distribution distribution of delays, covering length of
#'  infection timeseries, including incubation period, if applicable
#' @param proportion_observed long form data, for all dates of infection
#'  timeseries, of expected proportions of infections observed in count data
#' @param count_data long form data of counts of notifications
#' @param dow_model optional module of greta arrays defining day-of-week model
#' @param data_id optional name label identifying data type for greta arrays
#'
#' @importFrom greta %*% as_data negative_binomial normal sweep zeros
#'
#' @return greta arrays of observation model
#' @export
create_observation_model <- function (infection_timeseries,
                                      delay_distribution,
                                      proportion_observed,
                                      count_data,
                                      dow_model = NULL,
                                      data_id = NULL) {

  # add if statements to check that infection_days is long enough to cover
  # the right period

  case_mat <- data_to_matrix(count_data, 'value')
  delay_mat <- data_to_matrix(delay_distribution, 'value')
  prop_mat <- data_to_matrix(proportion_observed, 'value')

  infection_days <- as.Date(rownames(prop_mat))

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

  n_jurisdictions <- ncol(infection_timeseries)
  n_dates <- nrow(infection_timeseries)

  convolution_matrices <- lapply(
    1:n_jurisdictions,
    function(x) {
      get_convolution_matrix(delay_mat[, x],
                             n_dates)
    })

  # compute expected cases of the same length
  expected_cases_list <- lapply(
    1:n_jurisdictions,
    function(x) {
      convolution_matrices[[x]] %*% infection_timeseries[, x] * prop_mat[, x]
    })

  expected_cases <- do.call(
    cbind,
    expected_cases_list
  )

  data_idx <- infection_days %in% rownames(case_mat)
  expected_cases_idx <- expected_cases[data_idx, ]

  n_days <- nrow(case_mat)

  # negative binomial parameters
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
