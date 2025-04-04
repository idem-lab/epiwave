#' Define observation data
#'
#' @param timeseries_data timeseries data for data of interest
#' @param delay_from_infection delay data
#' @param proportion_infections proportion data
#' @param dow_model logical indicating whether to apply at DOW
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
