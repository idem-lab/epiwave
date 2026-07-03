#' Prepare observation data
#'
#' @description Prepare the data objects needed for one jurisdiction's
#'  observation model for a single data stream (e.g. cases). Multiple
#'  jurisdictions are combined later, by `stack_jurisdictions()`.
#'
#' @param observation_data data for one jurisdiction, one stream, as
#'  returned by `define_observation_data()`/`define_sero_data()`
#' @param target_infection_dates full date sequence of infection timeseries,
#'  shared across all jurisdictions in a fit
#'
#' @return a list describing this jurisdiction's stream: `convolution_matrix`,
#'  `case_vec`, `prop_vec`, `dow_model`, `inits_values`, `observable_idx`, and
#'  (for seroprevalence streams) `total_pop`/`size_vec`
#' @export
prepare_observation_data <- function (observation_data,
                                      target_infection_dates) {

  delays <- observation_data$delay_from_infection
  if (!('epiwave_massfun_timeseries' %in% class(delays))) {
    delays <- create_epiwave_massfun_timeseries(
      dates = target_infection_dates,
      value = delays)
  } else if (!identical(as.Date(delays$date), as.Date(target_infection_dates))) {
    stop('`delay_from_infection` dates must match `target_infection_dates`')
  }

  prop <- observation_data$proportion_infections
  if (!('epiwave_timeseries' %in% class(prop))) {
    prop <- create_epiwave_timeseries(
      dates = target_infection_dates,
      value = prop)
  } else if (!identical(as.Date(prop$date), as.Date(target_infection_dates))) {
    stop('`proportion_infections` dates must match `target_infection_dates`')
  }

  case_vec <- as_matrix(observation_data$timeseries_data, target_infection_dates)
  prop_vec <- as_matrix(prop, target_infection_dates)

  pmf_series <- epiwave.params::new_discrete_series(
    values = delays$value,
    index = delays$date
  )
  convolution_matrix <- new_convolution_matrix(pmf_series)

  inits <- inits_by_jurisdiction(
    case_vec,
    delays,
    prop_vec,
    target_infection_dates)

  inits_values <- rep(NA_real_, length(target_infection_dates))
  observable_idx <- rep(FALSE, length(target_infection_dates))
  inits_values[inits$observable_idx] <- inits$inits_values
  observable_idx[inits$observable_idx] <- TRUE

  out <- list(convolution_matrix = convolution_matrix,
              case_vec = case_vec,
              prop_vec = prop_vec,
              dow_model = observation_data$dow_model,
              inits_values = inits_values,
              observable_idx = observable_idx)

  if ('total_pop' %in% names(observation_data)) {
    out$total_pop <- observation_data$total_pop
  }
  if ('size_vec' %in% names(observation_data)) {
    out$size_vec <- as_matrix(observation_data$size_vec, target_infection_dates)
  }

  out
}
