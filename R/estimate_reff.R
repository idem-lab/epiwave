#' Estimate R_effective
#'
#' @description Calculate a timeseries of effective reproduction number, given a
#'  timeseries of infections, and the generation interval distribution. This
#'  function first convolves infections into infectiousness using the
#'  generation interval distribution. This weights each infectious individual
#'  by their relative infectivity profile for each day post infection. Then,
#'  infection timeseries is divided by the infectiousness timeseries to
#'  calculate instantaneous reproduction number. This interpretation of the
#'  instantaneous reproduction number is well established and consistent with
#'  classical discrete renewal models. Note that here the generation interval
#'  distribution acts as a proxy for infectivity profile distribution since the
#'  latter cannot be learned from typical epidemiological data. This is
#'  justified because across population, the distribution of the timing of
#'  secondary infections relative to the timing of primary infection (i.e. the
#'  generation interval) should correlate with when the primary individual is
#'  most infectious. Also note that because time of infection is rarely
#'  documented and generation interval is therefore hard to measure, in
#'  practice the serial interval is often used as a further proxy.
#'
#' @param infection_timeseries greta array of infection timeseries
#' @param generation_interval_mass_fxns mass functions describing the
#'  generation interval
#'
#' @importFrom greta %*%
#'
#' @return greta arrays for infectiousness and $R_eff$
#' @export
estimate_reff <- function (infection_timeseries,
                           generation_interval_mass_fxns) {

  convolution_matrix <- get_convolution_matrix(
      generation_interval_mass_fxns,
      nrow(infection_timeseries))
  infectiousness <- convolution_matrix %*% infection_timeseries

  reff <- infection_timeseries / (infectiousness + .Machine$double.eps)

  greta_arrays <- module(
    infectiousness,
    reff)

  return(greta_arrays)
}
