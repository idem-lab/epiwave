#' Construct delays
#'
#' @description Given a date range to compute delays over, and an empirical cdf,
#'  compute the probability of delay having that length. Outputs either the
#'  probability or the cumulative density, as either a data frame lined up
#'  against days-of-day or as a step function in the style of R native ecdf
#'  function. Optionally, two cdfs can be used as input, in which case the
#'  probability and cdf of the additive combined delay are calculated.
#'
#' @param ecdf_1 empirical cdf
#' @param ecdf_2 optional second empirical cdf to combine
#' @param delay_range number of days to create distribution for
#' @param output output style, choice of "probability" or "cumulative density"
#' @param stepfun_output logical whether or not to output as step function
#'
#' @importFrom tidyr expand_grid
#' @importFrom dplyr filter group_by summarise pull
#' @importFrom tibble tibble
#'
#' @return delay distribution in range of forms
#' @export
construct_delays <- function (ecdf_1,
                              ecdf_2 = NULL,
                              delay_range = c(-3, 28),
                              output = c("probability", "cumulative density"),
                              stepfun_output = FALSE) {

  if (is.list(ecdf_1)) {
    ecdf_1 <- ecdf_1[[1]]
  }

  # days of delay
  days <- seq(delay_range[1], delay_range[2])

  if (is.null(ecdf_2)) {

    # get discretised probability
    p <- ecdf_1(days + 1) - ecdf_1(days)

  } else {

    # get discretised probabilities for both cdfs
    p1 <- ecdf_1(days + 1) - ecdf_1(days)
    p2 <- ecdf_2(days + 1) - ecdf_2(days)

    p1 <- approxfun(days,p1, rule = 2)
    p2 <- approxfun(days,p2, rule = 2)

    # compute combined p
    p <- tidyr::expand_grid(
      x = days,
      z = days) %>%
      dplyr::filter(x - z >= 0) %>%
      dplyr::group_by(x) %>%
      dplyr::summarise(p = sum(p1(z ) * p2(x - z ))) %>%
      dplyr::pull(p)

  }

  # remove negative delay prob since notification cannot precede infection assuming that some
  # preventive measure has taken place once the "would-be" infectee is notified
  p[days < 0] <- 0

  # remove extremely long but unlikely delays, ususally a result of
  # parametric specification of the input cdf --- need to check if this is
  # right
  p[days > delay_range[2]] <- 0

  # normalise remaining probs
  p <- p / sum(p)

  if (output == "probability" & stepfun_output) {
    delay <- approxfun(days, p,
                       rule = 2,
                       method = "constant",
                       yleft = 0,
                       yright = 0,
                       f = 0)
  }

  if (output == "probability" & !stepfun_output) {
    delay <- tibble::tibble(days, p)
  }

  if (output != "probability" & stepfun_output) {
    delay <- make_ecdf(p, days)
  }

  if (output != "probability" & !stepfun_output) {
    delay <- tibble::tibble(days, cumsum(p))
  }

  return(delay)
}

#' Make empirical cumulative distribution function
#'
#' @description Convert a vector of cumulative probabilities into an ecdf object.
#'
#' @param y probability
#' @param x vector of one or more elements from which to choose
#'
#' @importFrom stats ecdf
#'
#' @return ecdf function
make_ecdf <- function (y, x) {

  sims <- sample(x,
                 100,
                 prob = y,
                 replace = TRUE)

  ecdf_null <- stats::ecdf(sims)
  envir <- environment(ecdf_null)

  # rebuild an ecdf object, the slow way
  method <- 2L
  yleft <- 0
  yright <- 1
  f <- envir$f
  n <- envir$nobs

  rval <- function (v) {
    stats:::.approxfun(x, y, v, method, yleft, yright, f)
  }

  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- attr(ecdf_null, "call")
  rval
}
