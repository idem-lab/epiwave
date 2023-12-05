#' Estimate R_effective
#'
#' @description Calculate a timeseries of effective reproduction number, given a
#'   timeseries of infections, and the generation interval distribution. This
#'   function first convolves infections into infectiousness using the
#'   generation interval distribution. This weights each infectious individual
#'   by their relative infectivity profile for each day post infection. Then,
#'   infection timeseries is divided by the infectiousness timeseries to
#'   calculate instantaneous reproduction number. This is a well established
#'   approach consistent with classical discrete renewal models, however note
#'   that the generation interval distribution acts as a proxy for infectivity
#'   profile distribution here, and should the serial interval distribution be
#'   available, it should be used instead.
#'
#'
#' @param infection_timeseries
#' @param generation_interval_mass_fxns
#'
#' @importFrom greta %*%
#'
#' @return greta arrays for infectiousness and $R_eff$
#' @export
estimate_reff <- function (infection_timeseries,
                           generation_interval_mass_fxns) {

  convolution_matrix <- get_convolution_matrix(
      generation_interval_mass_fxns, nrow(infection_timeseries))
  infectiousness <- convolution_matrix %*% infection_timeseries

  reff <- infection_timeseries / (infectiousness + .Machine$double.eps)

  greta_arrays <- module(
    infectiousness,
    reff)

  return(greta_arrays)

}
