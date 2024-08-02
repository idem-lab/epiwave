compute_reff <- function (wavylistything, gi) {

  ## Reff
  reff_model_objects <- estimate_reff(
    infection_timeseries = wavylistything$infection_model$infection_timeseries,
    generation_interval_mass_fxns = gi)

  reff_out <- epiwave.pipelines::generate_long_estimates(
    reff_model_objects$reff,
    wavylistything$fit,
    wavylistything$infection_days,
    wavylistything$jurisdictions)

  reff_out
}
