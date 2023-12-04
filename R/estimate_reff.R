#' Estimate R_effective
#'
#' @description A short description...
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
