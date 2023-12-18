#' Create day of week correction priors
#'
#' @description This function proposes prior distributions for the weights of
#'  day-of-week effect. Weekly periodicity is often observed in the case
#'  timeseries, and is likely typically caused by some genuine weekly cyclical
#'  pattern of infections (due to changes in contact pattern), and some changes
#'  in test seeking behaviour resulting in more test seeking on some days of
#'  the week. The latter is likely the dominant effect, so this function
#'  proposes DOW periodicity as an effect within the observation model, and
#'  constructs DOW weights such that infections on different days of the week
#'  would have longer or shorter than usual observation delay, thus producing a
#'  DOW oscillation in the convolved observation timeseries. The DOW weight
#'  prior is proposed as a 7-variate Dirichlet distribution, but alternative
#'  proposals may be implemented in the future.
#'
#' @param n_jurisdictions number of jurisdictions, defaults to 1
#'
#' @importFrom greta normal dirichlet
#'
#' @return named greta arrays for day of week corrections
#' @export
create_dow_priors <- function (n_jurisdictions = 1) {

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
