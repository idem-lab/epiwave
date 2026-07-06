#' Stack per-jurisdiction observation models together
#'
#' @description Combine several per-jurisdiction observation model bundles
#'  (each produced by `define_observation_model()`) into the
#'  `date x jurisdiction` matrices used by the model-fitting functions. This
#'  is where DOW correction is actually applied (it needs to know
#'  `n_jurisdictions`, which earlier, per-jurisdiction steps don't).
#'
#'  Only needed when combining more than one jurisdiction -- a
#'  single-jurisdiction fit can pass `define_observation_model()`'s output
#'  straight to `fit_waves()` without calling this at all.
#'
#'  Every jurisdiction must supply the same set of observation streams, and
#'  `dow_model` must agree across jurisdictions within a given stream; both
#'  are intentional limitations for now rather than permanent design
#'  decisions.
#'
#' @param ... per-jurisdiction observation model bundles (output of
#'  `define_observation_model()`), named by jurisdiction, e.g.
#'  `stack_jurisdictions(VIC = observation_model_vic, NSW = observation_model_nsw)`
#'
#' @return list of stacked observation model data, in the shape expected by
#'  `fit_waves()`, with class `epiwave_stacked_observations`
#' @export
stack_jurisdictions <- function (...) {
  stack_jurisdictions_list(list(...))
}

#' @noRd
stack_jurisdictions_list <- function (observations_by_jurisdiction) {

  jurisdictions <- names(observations_by_jurisdiction)

  if (length(observations_by_jurisdiction) == 1 &&
      (is.null(jurisdictions) || identical(jurisdictions, ""))) {
    # a single jurisdiction has no real identity to preserve -- fit_waves()
    # routes a raw define_observation_model() object through here unnamed
    jurisdictions <- "1"
  } else if (is.null(jurisdictions) || any(jurisdictions == "")) {
    stop('when combining more than one jurisdiction, all arguments to ',
         'stack_jurisdictions() must be named, named by jurisdiction')
  }

  target_infection_dates <- observations_by_jurisdiction[[1]]$target_infection_dates
  dates_match <- vapply(
    observations_by_jurisdiction,
    function(x) identical(as.Date(x$target_infection_dates),
                          as.Date(target_infection_dates)),
    logical(1))
  if (!all(dates_match)) {
    stop('All jurisdictions must share the same target_infection_dates')
  }

  stream_lists <- lapply(observations_by_jurisdiction,
                         function(x) names(x$observation_model_data))
  stream_ids <- stream_lists[[1]]
  streams_match <- vapply(
    stream_lists,
    function(s) identical(sort(s), sort(stream_ids)),
    logical(1))
  if (!all(streams_match)) {
    stop('Every jurisdiction must supply the same set of observation streams')
  }

  observation_model_data <- lapply(
    stream_ids,
    stack_stream,
    observations_by_jurisdiction,
    jurisdictions,
    target_infection_dates)
  names(observation_model_data) <- stream_ids

  out <- list(
    observation_model_data = observation_model_data,
    target_infection_dates = target_infection_dates,
    target_jurisdictions = jurisdictions
  )
  class(out) <- c("epiwave_stacked_observations", class(out))
  out
}

#' Stack one observation stream across jurisdictions
#'
#' @param stream_id name of the stream to stack (e.g. 'cases')
#' @param observations_by_jurisdiction named list of per-jurisdiction
#'  observation model bundles
#' @param jurisdictions jurisdiction labels, in stacking order
#' @param target_infection_dates full date sequence shared by all
#'  jurisdictions
#'
#' @return list describing this stream's stacked observation model data
#' @noRd
stack_stream <- function (stream_id,
                          observations_by_jurisdiction,
                          jurisdictions,
                          target_infection_dates) {

  per_jurisdiction <- lapply(
    observations_by_jurisdiction,
    function(x) x$observation_model_data[[stream_id]])

  case_mat <- do.call(cbind, lapply(per_jurisdiction, `[[`, "case_vec"))
  prop_mat <- do.call(cbind, lapply(per_jurisdiction, `[[`, "prop_vec"))
  colnames(case_mat) <- colnames(prop_mat) <- jurisdictions
  rownames(case_mat) <- rownames(prop_mat) <- as.character(target_infection_dates)

  # kept alongside the (possibly DOW-corrected) prop_mat below specifically
  # for compute_flat_prior_inits(): GAM-based inits must be computed from the
  # raw proportion, not the DOW-corrected one (DOW correction wasn't even
  # known/applicable at the point inits used to be computed, pre-refactor;
  # using the corrected version here would make inits depend on a
  # prior-predictive draw of the DOW effect instead of being deterministic)
  prop_mat_raw <- prop_mat

  dow_flags <- vapply(per_jurisdiction, `[[`, logical(1), "dow_model")
  if (length(unique(dow_flags)) > 1) {
    stop(sprintf(
      '`dow_model` must be the same for all jurisdictions within the "%s" stream',
      stream_id))
  }
  if (isTRUE(dow_flags[1])) {
    dow_priors <- create_dow_priors(length(jurisdictions))
    dow_correction <- implement_day_of_week(target_infection_dates, dow_priors)
    prop_mat <- prop_mat * dow_correction
  }

  convolution_matrices <- lapply(per_jurisdiction, `[[`, "convolution_matrix")
  delays <- lapply(per_jurisdiction, `[[`, "delays")

  list(
    convolution_matrices = convolution_matrices,
    delays = delays,
    case_mat = case_mat,
    prop_mat = prop_mat,
    prop_mat_raw = prop_mat_raw
  )
}
