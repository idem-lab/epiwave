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
#' @param delay_from_infection a `discrete_pmf` object (replicated across
#'  dates), or an already time-varying `discrete_pmf_series` object
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
