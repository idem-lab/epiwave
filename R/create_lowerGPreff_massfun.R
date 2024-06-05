#' Create lowerGPreff_massfun object
#'
#' @param min_delay miniumum delay
#' @param max_delay maxiumum delay
#' @param cdf_fun function for creating the cdf that defines the mass
#' @param normalise whether the mass should be normalised
#'
#' @return lowerGPreff_massfun object
#' @export
create_lowerGPreff_massfun <- function (min_delay,
                                        max_delay,
                                        cdf_fun,
                                        normalise = c(TRUE, FALSE)) {

  delay_massfun <- data.frame(
    delays = seq(min_delay, max_delay)
  ) %>%
    dplyr::mutate(
      upper = cdf_fun(delays + 1),
      lower = cdf_fun(delays),
      mass = upper - lower
    )

  if (normalise) {
    delay_massfun <- delay_massfun %>%
      dplyr::mutate(
        correction = sum(mass),
        mass = mass / correction
      )
  }

  delay_massfun <- delay_massfun %>%
    dplyr::select(
      delays, mass
    )

  class(delay_massfun) <- c("lowerGPreff_massfun", class(delay_massfun))

  if (normalise) {
    class(delay_massfun) <- c("lowerGPreff_distribution_massfun", class(delay_massfun))
  }
  if (!normalise) {
    class(delay_massfun) <- c("lowerGPreff_curve_massfun", class(delay_massfun))
  }

  delay_massfun

}

#' Create distribution object from data
#'
#' @param data delay data
#' @param min_delay optional specification for minimum delay
#' @param max_delay optional specification for maxiumum delay
#'
#' @return lowerGPreff_distribution object
#' @export
data_to_distribution <- function (data,
                                  min_delay = NULL,
                                  max_delay = NULL) {

  day_diff <- data$notif_date - data$sym_date
  cdf_fun <- stats::ecdf(day_diff)

  if (is.null(min_delay)) {
    min_delay <- floor(quantile(cdf_fun, 0.00))
  }
  if (is.null(max_delay)) {
    max_delay <- ceiling(quantile(cdf_fun, 0.99))
  }

  out <- create_lowerGPreff_massfun(
    min_delay, max_delay,
    cdf_fun, normalise = TRUE)
  out

}

#' Create distribution object from parametric distribution values
#'
#' @param dist distributional package object
#' @param min_delay optional specification for minimum delay
#' @param max_delay optional specification for maxiumum delay
#'
#' @importFrom distributional cdf
#'
#' @return lowerGPreff_distribution object
#' @export
parametric_dist_to_distribution <- function (dist,
                                             min_delay = NULL,
                                             max_delay = NULL) {

  # this quantile is the S3 method for quantile from distributional pkg
  if (is.null(min_delay)) {
    min_delay <- floor(quantile(dist, 0.00))
  }
  if (is.null(max_delay)) {
    max_delay <- ceiling(quantile(dist, 0.99))
  }

  # when distributional::cdf is applied to a sequence (of x) it returns a list
  cdf_fun <- function(x) distributional::cdf(dist, x)[[1]]

  out <- create_lowerGPreff_massfun(
    min_delay, max_delay,
    cdf_fun, normalise = TRUE)
  out

}

#' Default plot function for lowerGPreff_massfun object
#'
#' @param x x value
#' @param y y value
#' @param ... additional arguments
#'
#' @return plot
#'
#' @export
plot.lowerGPreff_massfun <- function (x, y, ...) {
  barplot(x$mass, width = 1, names.arg = x$delays,
          xlab = "delay (days)",
          ylab = "probability")
}
