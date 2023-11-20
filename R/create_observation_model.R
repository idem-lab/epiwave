#' Create observation model
#'
#' @param infection_timeseries
#' @param infection_delay_distribution
#' @param timevarying_proportion
#' @param count_data
#' @param data_id optional name label identifying data type for greta arrays
#'
#' @return
#' @export
create_observation_model <- function (infection_timeseries,
                                      infection_delay_distribution,
                                      timevarying_proportion, # long format dataframes/matrices
                                      count_data, # specific data type
                                      data_id) {

  case_mat <- data_to_matrix(count_data, 'count')

  valid_mat <- count_data
  valid_mat[is.na(count_data)] <- FALSE
  valid_mat[!is.na(count_data)] <- TRUE

  delay_mat <- data_to_matrix(infection_delay_distribution, 'delay')
  prop_mat <- data_to_matrix(timevarying_proportion, 'prop')

  n_jurisdictions <- ncol(infection_timeseries)

  ext_data_idx <- infection_timeseries_dates %in%
    rownames(infection_delay_distribution)
  infection_timeseries_ext_data <- infection_timeseries[ext_data_idx, ]

  n_ext_dates <- nrow(infection_timeseries_ext_data)

  convolution_matrices <- lapply(1:n_jurisdictions, function(x)
    get_convolution_matrix(infection_delay_distribution[, x],
                           n_ext_dates))

  # compute expected cases of the same length
  # note not all of these dates would have been observed
  expected_cases <- do.call(
    cbind,
    lapply(1:n_jurisdictions, function(x) {
      convolution_matrices[[x]] %*% infection_timeseries_ext_data[, x] *
        timevarying_proportion[, x]
    }))

  # here index expected cases to data that is brought in.
  # # get indices for subset of infection days that end up in observed/forecast data
  # observable_infection_idx <- which(full_infection_dates %in% observable_infection_dates)
  #
  # # subset infection timeseries to these indices
  # # note this is the number of infection days that get observed in the observed
  # # data, it is not the same length as the actual number of observed dates in
  # # notification series, which is below
  # infection_observable <- infections_timeseries[observable_infection_idx,]
  #
  # # observed days in the obs data itself
  # # basically this backs out extra left and right days
  # obs_data_idx <- which(observable_infection_dates %in% as.Date(rownames(count_data)))
  # # build convolution matrix
  # n_days_infection_observable <- length(observable_infection_idx)



  # use validity matrix - default is all valid but can use this to nullify
  # expected mean cases for dates when reporting is stopped
  # extend the valid mat to include forecast days by repeating the last row
  # this ensures for states where data collection has stopped, the forecast is also 0

  # if (!is.matrix(valid_mat)) {
  #
  #     valid_mat <- count_data
  #     valid_mat[] <- TRUE
  #     valid_idx <- which(as.logical(valid_mat))
  # } else {
  #     valid_idx <- which(valid_mat,arr.ind = FALSE)
  # }

  data_idx <- rownames(infection_delay_distribution) %in%
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
