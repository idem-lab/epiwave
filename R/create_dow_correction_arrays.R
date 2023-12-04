#' Create day of week correction arrays
#'
#' @description A short description...
#'
#'
#' @param n_jurisdictions number of jurisdictions, defaults to 1
#' @param data_id optional name label identifying data type for greta arrays
#'
#' @return named greta arrays for day of week corrections
#' @export
create_dow_correction_arrays <- function (n_jurisdictions = 1,
                                          data_id = '') {

  # prior for dweek correction
  dow_alpha <- greta::normal(1, 1,
                             truncation = c(0, Inf),
                             dim = c(1, 7))

  dow_dist <- greta::dirichlet(dow_alpha,
                               n_realisations = n_jurisdictions)

  greta_arrays <- module(
    dow_alpha,
    dow_dist)

  return(greta_arrays)
}
