#' copied from the main function here
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
#' @param data_id name of observation model data in list
#' @param observation_model_data list of observation model data lists
#' @param infection_days infection dates that cover more than the data dates
#' @param infection_model greta arrays that define the infection model
#'
#' @importFrom greta %*% as_data negative_binomial normal sweep zeros
#'
#' @return greta arrays of observation model
#' @export
create_observation_model <- function (data_id = 'sero',
                                      observation_model_data,
                                      infection_days,
                                      infection_model
) {

  observations <- observation_model_data[[data_id]]

  delays <- observations$delays # technically conversions, but keeping name consistent

  # this needs to be only the days we actually conduct sero survey
  case_mat <- observations$case_mat

  # if data is sero, needs to have a size mat
  size_mat <- observations$size_mat

  prop_mat <- observations$prop_mat # assuming this is 1 throughout

  n_dates <- length(infection_days)

  convolution_matrices <- lapply(
    unique(delays$jurisdiction),
    function(x) {
      get_convolution_matrix(delays,
                             x,
                             n_dates)
    })

  n_jurisdictions <- length(unique(delays$jurisdiction))

  expected_cases_list <- lapply(
    1:n_jurisdictions,
    function(x) {
      (convolution_matrices[[x]]) %*% infection_model[, x] * prop_mat[, x]
    })

  expected_cases <- do.call(
    cbind,
    expected_cases_list
  )

  data_idx <- infection_days %in% as.Date(rownames(case_mat))
  expected_cases_idx <- expected_cases[data_idx, ]

  # remember to MAKE this
  # a vector of population by jurisdiction (static)
  total_pop <- observation_model_data$total_pop

  prob <- sweep(expected_cases_idx,2,total_pop,"/")

  n_days <- nrow(case_mat)


  valid_mat <- case_mat
  valid_mat[is.na(case_mat)] <- FALSE
  valid_mat[!is.na(case_mat)] <- TRUE
  valid_idx <- as.logical(as.numeric(valid_mat))

  case_mat_array <- greta::as_data(
    as.numeric(case_mat)[valid_idx])

  size_mat_array <- greta::as_data(
    as.numeric(size_mat)[valid_idx])

  greta::distribution(case_mat_array) <- greta::binomial(
    size_mat_array,
    prob[valid_idx])

  greta_arrays <- list(
    #size,
    prob,
    convolution_matrices
  ) # diff output depending on data_id type

  names(greta_arrays) <- c(
    #paste0(data_id, '_size'),
    paste0(data_id, '_prob'),
    paste0(data_id, '_convolution_matrices')
  )

  return(greta_arrays)

}
