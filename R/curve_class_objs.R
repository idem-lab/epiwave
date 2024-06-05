#' Create distribution object from parametric distribution values
#'
#' @param dist distributional package object
#' @param min_delay optional specification for minimum delay
#' @param max_delay optional specification for maxiumum delay
#'
#' @return lowerGPreff_distribution object
#' @export
parametric_dist_to_curve <- function (dist,
                                      min_delay = NULL,
                                      max_delay = NULL) {

  # this quantile is the S3 method for quantile from distributional pkg
  if (is.null(min_delay)) {
    min_delay <- floor(quantile(dist, 0))
  }
  if (is.null(max_delay)) {
    max_delay <- ceiling(quantile(dist, 1))
  }

  # when distributional::cdf is applied to a sequence (of x) it returns a list
  cdf_fun <- function(x) distributional::cdf(dist, x)[[1]]

  out <- create_lowerGPreff_massfun(
    min_delay, max_delay,
    cdf_fun, normalise = FALSE)
  out

}


#' CHANGE ME Create distribution object from parametric distribution values
#'
#' @description For sero data the value col will be probability of testing
#'  positive. For wastewater data the value col will be amount of viral
#'  material shed by individual compared to peak amount.
#'
#' @param day_diff_col days since onset of infection
#' @param value_col value on that many days after infection onset
#'
#' @return lowerGPreff_distribution object
#' @export
data_to_curve <- function (day_diff_col,
                           value_col) {

  curve_df <- data.frame(delays = day_diff_col,
                         mass = value_col)

  class(curve_df) <- c("lowerGPreff_curve_massfun",
                            "lowerGPreff_massfun",
                            class(curve_df))

  curve_df
}

