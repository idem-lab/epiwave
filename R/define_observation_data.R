#' Define observation data
#'
#' @description Bundle one jurisdiction's data for a single observation
#'  stream (e.g. cases). Call this once per jurisdiction per stream; combine
#'  multiple jurisdictions later via `stack_jurisdictions()`.
#'
#' @param timeseries_data timeseries data for data of interest, for one
#'  jurisdiction
#' @param delay_from_infection delay data
#' @param proportion_infections proportion data
#' @param dow_model logical indicating whether to apply a DOW
#'
#' @return list of observation data for one data type
#' @export
define_observation_data <- function (timeseries_data,
                                     delay_from_infection,
                                     proportion_infections,
                                     dow_model = FALSE) {

  out <- list(timeseries_data = timeseries_data,
              delay_from_infection = delay_from_infection,
              proportion_infections = proportion_infections,
              dow_model = dow_model)
  out

}

#' Define seroprevalence observation data
#'
#' @description Bundle one jurisdiction's seroprevalence survey data. Call
#'  this once per jurisdiction; combine multiple jurisdictions later via
#'  `stack_jurisdictions()`.
#'
#' @param timeseries_data seroprevalence survey timeseries data, for one
#'  jurisdiction
#' @param total_pop population size of this jurisdiction (a single numeric
#'  value)
#' @param size_vec sample size of the seroprevalence survey, for one
#'  jurisdiction
#' @param delay_from_infection delay data
#' @param proportion_infections proportion data
#' @param dow_model logical indicating whether to apply a DOW
#'
#' @return list of observation data for one data type
define_sero_data <- function (timeseries_data,
                              total_pop,
                              size_vec,
                              delay_from_infection,
                              proportion_infections,
                              dow_model = FALSE) {

  out <- list(timeseries_data = timeseries_data,
              total_pop = total_pop,
              size_vec = size_vec,
              delay_from_infection = delay_from_infection,
              proportion_infections = proportion_infections,
              dow_model = dow_model)
  out

}
