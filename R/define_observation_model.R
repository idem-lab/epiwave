#' Define observation model
#'
#' @description
#' Bundle one jurisdiction's observation streams (e.g. cases,
#' hospitalisations, sero) together. Call this once per jurisdiction;
#' combine multiple jurisdictions later via `stack_jurisdictions()`.
#'
#' @param target_infection_dates sequence of infection dates
#' @param ... observation data sets for this jurisdiction (as returned by
#'  `define_observation_data()`/`define_sero_data()`), named by stream
#'
#' @return list describing one jurisdiction's observation model
#' @export
#'
define_observation_model <- function (target_infection_dates = NULL, ...) {

  observation_list <- list(...)

  prepared_observation_model_data <- lapply(observation_list,
                                            prepare_observation_data,
                                            target_infection_dates)

  incidence_observable <- Reduce(
    `|`,
    lapply(prepared_observation_model_data, function(x) x$observable_idx))

  inits_stream_mat <- do.call(
    cbind,
    lapply(prepared_observation_model_data, function(x) x$inits_values))
  incidence_observable_inits <- rowMeans(inits_stream_mat, na.rm = TRUE)

  list(observation_model_data = prepared_observation_model_data,
       target_infection_dates = target_infection_dates,
       incidence_observable_inits = incidence_observable_inits,
       incidence_observable = incidence_observable)

}
