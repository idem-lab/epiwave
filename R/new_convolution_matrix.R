#' Construct a forward convolution matrix from a discrete PMF or weights
#'
#' @description A common operation in discrete-time epidemic models is forward
#' convolution through a delay distribution. This can be represented as a
#' linear operator in matrix form. This function constructs a forward
#' convolution matrix from a `discrete_pmf`/`discrete_pmf_series` object, or
#' from a `discrete_weights`/`discrete_weights_series` object (e.g. an
#' unnormalised persistence/detectability curve, such as the probability of
#' testing seropositive by day since infection). The resulting matrix can be
#' multiplied by a time series vector to convolve it forward through a delay
#' distribution.
#'
#' When `pmf` is a single `discrete_pmf`/`discrete_weights`, the same object
#' is applied at every timepoint and `n` must be supplied. When `pmf` is a
#' `discrete_pmf_series`/`discrete_weights_series`, `n` is derived from the
#' series index and does not need to be supplied.
#'
#' @srrstats {PD3.5} Discrete summation is used following Cori et al. (2013,
#'   AJE) and Gostic et al. (2020, PLOS Comp Bio), who demonstrate its
#'   suitability for discrete-time epidemiological models. The input PMF or
#'   weight function is finite by construction, guaranteeing a finite sum
#'   (PD3.5a).
#'
#' @param pmf a `discrete_pmf`, `discrete_weights`, `discrete_pmf_series`, or
#'   `discrete_weights_series` object. A single `discrete_pmf`/
#'   `discrete_weights` is applied uniformly across all timepoints. A
#'   `discrete_pmf_series`/`discrete_weights_series` enables time-varying
#'   delays, with `n` derived from the series index.
#' @param n a positive integer giving the number of timepoints, determining
#'   the dimensions of the output matrix. Required when `pmf` is a single
#'   `discrete_pmf`/`discrete_weights` object; inferred from the timeseries
#'   index when `pmf` is a series object.
#'
#' @return an `n` by `n` numeric matrix. To apply the convolution, multiply
#'   the matrix by a time series vector of length `n`, e.g.
#'   `new_convolution_matrix(pmf, n) %*% time_series`.
#'
#' @examples
#' pmf <- epiwave.params::as_discrete_pmf(
#'        distributional::dist_gamma(shape = 3, rate = 0.5))
#' new_convolution_matrix(pmf, n = 10)
#'
#' @export
new_convolution_matrix <- function(pmf, n = NULL) {

  if (!inherits(pmf, c("discrete_pmf", "discrete_pmf_series",
                       "discrete_weights", "discrete_weights_series"))) {
    cli::cli_abort(
      "`pmf` must be a {.cls discrete_pmf}, {.cls discrete_pmf_series}, {.cls discrete_weights}, or {.cls discrete_weights_series} object."
    )
  }

  n <- resolve_n(pmf, n)

  day_diff <- matrix(NA, n, n)
  day_diff <- row(day_diff) - col(day_diff)

  evaluate(pmf, day_diff)
}

#' @noRd
resolve_n <- function(pmf, n) {
  if (inherits(pmf, c("discrete_pmf_series", "discrete_weights_series"))) {
    if (!is.null(n)) {
      cli::cli_warn(
        "`n` is ignored when `pmf` is a series object; derived from `pmf$index`."
      )
    }
    # length(pmf$index) is guaranteed valid by new_discrete_series construction
    return(length(pmf$index))
  }
  validate_n(n)
}

#' @noRd
validate_n <- function(n) {
  if (is.null(n)) {
    cli::cli_abort(
      "`n` must be supplied when `pmf` is a single {.cls discrete_pmf} or {.cls discrete_weights} object."
    )
  }
  if (!is.numeric(n) || length(n) != 1) {
    cli::cli_abort("`n` must be a single numeric value.")
  }
  if (is.infinite(n) || is.nan(n)) {
    cli::cli_abort("`n` must be a finite value.")
  }
  if (n < 1 || n != as.integer(n)) {
    cli::cli_abort("`n` must be a positive integer.")
  }
  invisible(n)
}
