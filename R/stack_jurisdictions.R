#' Stack per-jurisdiction observation models together
#'
#' @description Combine a named list of per-jurisdiction observation model
#'  bundles (each produced by `define_observation_model()`) into the
#'  `date x jurisdiction` matrices used by the model-fitting functions. A
#'  single-jurisdiction fit is just a length-1 list. This is where DOW
#'  correction is actually applied (it needs to know `n_jurisdictions`, which
#'  earlier, per-jurisdiction steps don't).
#'
#'  Every jurisdiction must supply the same set of observation streams, and
#'  `dow_model` must agree across jurisdictions within a given stream; both
#'  are intentional limitations for now rather than permanent design
#'  decisions.
#'
#' @param observations_by_jurisdiction a named list of per-jurisdiction
#'  observation model bundles (output of `define_observation_model()`),
#'  named by jurisdiction
#'
#' @return list of stacked observation model data, in the shape expected by
#'  `fit_waves()`
#' @export
stack_jurisdictions <- function (observations_by_jurisdiction) {

  jurisdictions <- names(observations_by_jurisdiction)
  if (is.null(jurisdictions) || any(jurisdictions == "")) {
    stop('`observations_by_jurisdiction` must be a fully named list, ',
         'named by jurisdiction')
  }

  target_infection_dates <- observations_by_jurisdiction[[1]]$target_infection_dates
  dates_match <- vapply(
    observations_by_jurisdiction,
    function(x) identical(as.Date(x$target_infection_dates), as.Date(target_infection_dates)),
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

  incidence_observable <- do.call(
    cbind,
    lapply(observations_by_jurisdiction, `[[`, "incidence_observable"))
  incidence_observable_inits <- do.call(
    cbind,
    lapply(observations_by_jurisdiction, `[[`, "incidence_observable_inits"))
  colnames(incidence_observable) <- jurisdictions
  colnames(incidence_observable_inits) <- jurisdictions

  observation_model_data <- lapply(
    stream_ids,
    stack_stream,
    observations_by_jurisdiction,
    jurisdictions,
    target_infection_dates)
  names(observation_model_data) <- stream_ids

  list(
    observation_model_data = observation_model_data,
    target_infection_dates = target_infection_dates,
    target_jurisdictions = jurisdictions,
    incidence_observable = incidence_observable,
    incidence_observable_inits = incidence_observable_inits
  )
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

  sero_flags <- vapply(per_jurisdiction, function(x) !is.null(x$total_pop), logical(1))
  if (length(unique(sero_flags)) > 1) {
    stop(sprintf(
      'Jurisdictions disagree on whether "%s" is a seroprevalence stream ',
      '(define_sero_data() vs define_observation_data())',
      stream_id))
  }

  out <- list(
    convolution_matrices = convolution_matrices,
    case_mat = case_mat,
    prop_mat = prop_mat
  )

  if (isTRUE(sero_flags[1])) {
    out$total_pop <- vapply(per_jurisdiction, `[[`, numeric(1), "total_pop")
    names(out$total_pop) <- jurisdictions
    out$size_mat <- do.call(cbind, lapply(per_jurisdiction, `[[`, "size_vec"))
    colnames(out$size_mat) <- jurisdictions
    rownames(out$size_mat) <- as.character(target_infection_dates)
  }

  out
}
