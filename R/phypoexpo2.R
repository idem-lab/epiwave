#' distribution function for 2-component hypoexponential distribution, with the
#' two rates as its parameters. A hypoexponential distribution is the sum of two
#' exponentially distributed RV with different rates. A special case when the
#' two rates are the same is the Erlang distribution.
#'
#' @param q vector of quantiles
#' @param rate_1 exponential rate for component 1
#' @param rate_2 exponential rate for component 2
#'
#' @return vector of probabilities
#' @export
#'
#' @examples
#'
phypoexpo2 <- function(q, rate_1, rate_2) {
  denom <- rate_2 - rate_1
  1 - (exp(-rate_1 * q) * rate_2 / denom) + (exp(-rate_2 * q) * rate_1 / denom)
}
