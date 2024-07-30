#' Prepare observation data
#'
#' @description Prepare list of data objects needed for observation model.
#'
#' @param timeseries_data long form data of counts of notifications
#' @param delays distribution of delays, covering length of
#'  infection timeseries, including incubation period, if applicable
#' @param proportion_observed long form data, for all dates of infection
#'  timeseries, of expected proportions of infections observed in count data
#' @param type count or prevalence
#' @param dow_model optional module of greta arrays defining day-of-week model
#'
#' @importFrom greta %*% as_data negative_binomial normal sweep zeros
#'
#' @return greta arrays of observation model
#' @export
prepare_observation_data <- function (timeseries_data,
                                      delays,
                                      proportion_observed,
                                      type = c('count', 'prevalence'),
                                      dow_model = NULL) {

  # see how we are defining the likelihood
  type <- match.arg(type)

  # check the data is a valid timeseries object (or can be coerced to one) and
  # that it contains data consistent with the observation model type
  # timeseries_data <- check_valid_timeseries_data(timeseries_data,
  #                                                type = type)

  # add if statements to check that infection_days is long enough to cover
  # the right period

  case_mat <- as_matrix(timeseries_data)
  prop_mat <- as_matrix(proportion_observed)

  obs_infection_days <- as.Date(rownames(prop_mat))

  ## add a check for correct dow arrays
  if (!is.null(dow_model)) {
    dow_correction <- implement_day_of_week(obs_infection_days, dow_model)
    prop_mat <- prop_mat * dow_correction
  }

  out <- list(timeseries_data = timeseries_data,
              delays = delays,
              case_mat = case_mat,
              prop_mat = prop_mat)
  out

}
