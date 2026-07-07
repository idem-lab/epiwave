#' Stack per-jurisdiction observation models together
#'
#' @description Combine several per-jurisdiction observation model bundles
#'  (each produced by `define_observation_model()`) into the
#'  `date x jurisdiction` matrices used by the model-fitting functions.
#'
#'  This is where the shared `target_infection_dates` axis is actually
#'  decided -- it's emergent, not something the caller has to compute: it's
#'  the union of every jurisdiction's implied date range (each jurisdiction's
#'  own range is itself the union of its streams' implied ranges -- see
#'  `prepare_observation_data()`), optionally extended forward by
#'  `forecast_horizon`. Every stream's data stays in its own native,
#'  date-keyed form right up until this point (real `Date`s, not
#'  positions) -- this function is also where each stream is finally
#'  aligned to the shared axis (by matching dates, not by assuming
#'  pre-alignment) and where convolution matrices get built, since both
#'  need to know the axis, which isn't decided any earlier than this.
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
#' @param forecast_horizon optional number of days to extend the emergent
#'  date axis beyond the latest date implied by the data (e.g. for
#'  nowcasting/forecasting). Defaults to `0` (no extension) -- how far
#'  forward to estimate is a modelling choice, not something derivable from
#'  the data.
#'
#' @return list of stacked observation model data, in the shape expected by
#'  `fit_waves()`, with class `epiwave_stacked_observations`
#' @export
stack_jurisdictions <- function (..., forecast_horizon = 0) {
  stack_jurisdictions_list(list(...), forecast_horizon = forecast_horizon)
}

#' @noRd
stack_jurisdictions_list <- function (observations_by_jurisdiction, forecast_horizon = 0) {

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

  # the emergent date axis: union of every jurisdiction's implied range,
  # optionally extended forward for a nowcast/forecast horizon
  starts <- do.call(c, lapply(observations_by_jurisdiction, function(x) x$implied_range[1]))
  ends <- do.call(c, lapply(observations_by_jurisdiction, function(x) x$implied_range[2]))
  target_infection_dates <- seq(min(starts), max(ends) + forecast_horizon, by = "day")

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

#' Align and stack one observation stream across jurisdictions
#'
#' @description Each jurisdiction's stream arrives still in its native,
#'  date-keyed form (see `prepare_observation_data()`) -- this is where it's
#'  finally aligned to the shared `target_infection_dates` axis (by
#'  matching real dates, not by positional assumption), where the
#'  convolution matrix is built (its size depends on the axis, so it can't
#'  be built any earlier), and where jurisdictions are combined into
#'  `date x jurisdiction` matrices.
#'
#' @param stream_id name of the stream to stack (e.g. 'cases')
#' @param observations_by_jurisdiction named list of per-jurisdiction
#'  observation model bundles
#' @param jurisdictions jurisdiction labels, in stacking order
#' @param target_infection_dates the emergent, shared date sequence
#'
#' @return list describing this stream's stacked observation model data
#' @noRd
stack_stream <- function (stream_id,
                          observations_by_jurisdiction,
                          jurisdictions,
                          target_infection_dates) {

  per_jurisdiction_raw <- lapply(
    observations_by_jurisdiction,
    function(x) x$observation_model_data[[stream_id]])

  per_jurisdiction <- lapply(per_jurisdiction_raw, function (stream) {
    align_stream_to_axis(stream, target_infection_dates)
  })

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

#' Align one jurisdiction's stream onto the shared date axis
#'
#' @param stream the native, date-keyed stream object (output of
#'  `prepare_observation_data()`)
#' @param target_infection_dates the emergent, shared date sequence
#'
#' @return list with `case_vec`, `prop_vec`, `delays` (aligned
#'  `discrete_pmf_series`), `convolution_matrix`, `dow_model`
#' @noRd
align_stream_to_axis <- function (stream, target_infection_dates) {

  delays <- stream$delays
  if (inherits(delays, 'discrete_pmf_series')) {
    delays <- delays[as.Date(target_infection_dates)]
    if (!identical(as.Date(delays$index), as.Date(target_infection_dates))) {
      stop('`delay_from_infection` dates must cover the derived ',
           'target_infection_dates axis')
    }
  } else {
    delays <- epiwave.params::new_discrete_series(
      values = delays,
      index = target_infection_dates)
  }

  prop <- stream$proportion_infections
  if (inherits(prop, 'greta_proportion')) {
    # a greta-backed proportion (e.g. IHR = CAR x CHR) can only be
    # dimensioned once the axis is known -- resolve it now
    car <- prop$car
    chr_prior <- prop$chr_prior
    dim(chr_prior) <- length(target_infection_dates)
    prop_vec <- car * chr_prior
  } else {
    if (inherits(prop, 'epiwave_timeseries')) {
      if (!identical(as.Date(prop$date), as.Date(target_infection_dates))) {
        stop('`proportion_infections` dates must cover the derived ',
             'target_infection_dates axis')
      }
    }
    prop_vec <- as_matrix(prop, target_infection_dates)
  }

  case_vec <- as_matrix(stream$timeseries_data, target_infection_dates)
  convolution_matrix <- new_convolution_matrix(delays)

  list(
    case_vec = case_vec,
    prop_vec = prop_vec,
    delays = delays,
    convolution_matrix = convolution_matrix,
    dow_model = stream$dow_model
  )
}
