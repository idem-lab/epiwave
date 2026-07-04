#' Create small seroprevalence observation model
#'
#' @description This function constructs the observation model and defines
#'  likelihood over seroprevalence survey data. It mirrors
#'  `create_observation_model()`: a timeseries of expected observations is
#'  computed from the infection timeseries, the infection-to-seroconversion
#'  delay distribution(s), and the proportion of infections to be observed.
#'  Unlike case/hospitalisation notification (a one-time event best modelled
#'  as a `discrete_pmf`), seroconversion is typically persistent -- an
#'  infected person may test positive for many consecutive days -- so
#'  `delay_from_infection` for a sero stream is usually a `discrete_weights`/
#'  `discrete_weights_series` object (e.g. probability of testing seropositive
#'  by day since infection), rather than a normalised `discrete_pmf`. The
#'  expected observation timeseries, divided by jurisdiction population,
#'  is treated as the probability of a binomial distribution over the survey
#'  sample size, from which the observed seropositive count is drawn.
#'
#' @param data_id name of observation model data in list
#' @param observation_model_data list of observation model data lists, as
#'  produced by `stack_jurisdictions()`
#' @param infection_model greta arrays that define the infection model
#'
#' @importFrom greta %*% as_data binomial sweep
#'
#' @return greta arrays of observation model
#' @export
create_small_sero_model <- function (data_id = 'sero',
                                      observation_model_data,
                                      infection_model
) {

  observations <- observation_model_data[[data_id]]

  case_mat <- observations$case_mat
  size_mat <- observations$size_mat
  prop_mat <- observations$prop_mat
  convolution_matrices <- observations$convolution_matrices
  total_pop <- observations$total_pop

  n_jurisdictions <- length(convolution_matrices)

  expected_cases_list <- lapply(
    1:n_jurisdictions,
    function(x) {
      (convolution_matrices[[x]]) %*% infection_model[, x] * prop_mat[, x]
    })

  expected_cases <- do.call(
    cbind,
    expected_cases_list
  )

  prob <- sweep(expected_cases, 2, total_pop, "/")

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
    prob,
    convolution_matrices
  )

  names(greta_arrays) <- c(
    paste0(data_id, '_prob'),
    paste0(data_id, '_convolution_matrices')
  )

  return(greta_arrays)

}
