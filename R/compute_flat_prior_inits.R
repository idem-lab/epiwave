#' Compute initial values for the flat_prior infection model
#'
#' @description GAM-smoothed initial values and the observable-date window
#'  are only used by `infection_model_type = 'flat_prior'` -- GP-based
#'  infection models (`gp_infections`/`gp_growth_rate`/
#'  `gp_growth_rate_deriv`) never reference `observable_idx`/`inits_values`.
#'  This is deliberately not computed during data preparation
#'  (`prepare_observation_data()`/`define_observation_model()`/
#'  `stack_jurisdictions()`), since doing so unconditionally would run an
#'  `mgcv::gam()` fit per stream per jurisdiction for no purpose whenever a
#'  GP model is used instead. `fit_waves()` calls this lazily, only inside
#'  its `flat_prior` branch, on the fully stacked multi-jurisdiction object.
#'
#' @param observations a stacked observations object, as produced by
#'  `stack_jurisdictions()`
#'
#' @return list with `incidence_observable` (`date x jurisdiction` logical
#'  matrix, `TRUE` where at least one stream has nearby information to infer
#'  an infection) and `incidence_observable_inits` (`date x jurisdiction`
#'  numeric matrix, mean GAM-smoothed initial value across streams)
#' @noRd
compute_flat_prior_inits <- function (observations) {

  target_infection_dates <- observations$target_infection_dates
  jurisdictions <- observations$target_jurisdictions
  n_dates <- length(target_infection_dates)
  stream_ids <- names(observations$observation_model_data)

  per_jurisdiction <- lapply(seq_along(jurisdictions), function (j) {

    stream_inits <- lapply(stream_ids, function (stream_id) {
      stream <- observations$observation_model_data[[stream_id]]
      # prop_mat_raw, not prop_mat: inits must use the pre-DOW-correction
      # proportion (see stack_stream()) -- using the corrected one would make
      # this depend on a prior-predictive draw of the DOW effect instead of
      # being a deterministic approximation
      inits_by_jurisdiction(
        stream$case_mat[, j],
        stream$delays[[j]],
        stream$prop_mat_raw[, j],
        target_infection_dates)
    })

    observable_idx <- rep(FALSE, n_dates)
    inits_stream_mat <- matrix(NA_real_, n_dates, length(stream_ids))

    for (s in seq_along(stream_inits)) {
      idx <- stream_inits[[s]]$observable_idx
      observable_idx[idx] <- TRUE
      inits_stream_mat[idx, s] <- stream_inits[[s]]$inits_values
    }

    list(
      observable_idx = observable_idx,
      inits_values = rowMeans(inits_stream_mat, na.rm = TRUE)
    )
  })

  incidence_observable <- do.call(
    cbind, lapply(per_jurisdiction, `[[`, "observable_idx"))
  incidence_observable_inits <- do.call(
    cbind, lapply(per_jurisdiction, `[[`, "inits_values"))
  colnames(incidence_observable) <- jurisdictions
  colnames(incidence_observable_inits) <- jurisdictions

  list(
    incidence_observable = incidence_observable,
    incidence_observable_inits = incidence_observable_inits
  )
}
