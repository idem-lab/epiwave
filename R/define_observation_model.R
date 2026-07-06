#' Define observation model
#'
#' @description
#' Bundle one jurisdiction's observation streams (e.g. cases,
#' hospitalisations) together. Call this once per jurisdiction. A
#' single-jurisdiction fit can pass this straight to `fit_waves()`; to
#' combine multiple jurisdictions, pass several of these to
#' `stack_jurisdictions()` first.
#'
#' @param target_infection_dates sequence of infection dates
#' @param ... observation data sets for this jurisdiction (as returned by
#'  `define_observation_data()`), named by stream
#'
#' @return list describing one jurisdiction's observation model, with class
#'  `epiwave_observation_model`
#' @export
#'
define_observation_model <- function (target_infection_dates = NULL, ...) {

  observation_list <- list(...)

  prepared_observation_model_data <- lapply(observation_list,
                                            prepare_observation_data,
                                            target_infection_dates)

  out <- list(observation_model_data = prepared_observation_model_data,
             target_infection_dates = target_infection_dates)

  class(out) <- c("epiwave_observation_model", class(out))
  out

}
