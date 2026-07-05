#' Define observation data
#'
#' @description Bundle one jurisdiction's data for a single observation
#'  stream (e.g. cases). Call this once per jurisdiction per stream; combine
#'  multiple jurisdictions later via `stack_jurisdictions()`.
#'
#' @param timeseries_data timeseries data for data of interest, for one
#'  jurisdiction. A plain data.frame/tibble with `date` and `value` columns
#'  is fine (gaps in date coverage are fine too) -- it doesn't need to be
#'  pre-classed, see `as_epiwave_timeseries()`
#' @param delay_from_infection a `discrete_pmf`/`discrete_weights` object
#'  (replicated across dates), or an already time-varying
#'  `discrete_pmf_series`/`discrete_weights_series` object
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
#'  jurisdiction. A plain data.frame/tibble with `date` and `value` columns
#'  is fine, see `as_epiwave_timeseries()`
#' @param total_pop population size of this jurisdiction (a single numeric
#'  value)
#' @param size_vec sample size of the seroprevalence survey, for one
#'  jurisdiction. A plain data.frame/tibble with `date` and `value` columns
#'  is fine, see `as_epiwave_timeseries()`
#' @param delay_from_infection typically a `discrete_weights`/
#'  `discrete_weights_series` object (e.g. probability of testing
#'  seropositive by day since infection, which need not sum to 1 -- unlike
#'  case/hospitalisation notification, seroconversion is usually persistent
#'  rather than a one-time event). A `discrete_pmf`/`discrete_pmf_series` is
#'  also accepted if a normalised delay is more appropriate.
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
