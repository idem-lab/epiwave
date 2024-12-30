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

  if (n_jurisdictions == 1) {

    # nation-wide alpha
    dow_alpha <- greta::normal(1, 1,
                               truncation = c(0, Inf),
                               dim = c(1, 6))
  }

  if (n_jurisdictions > 1) {

    # nation-wide alpha
    dow_alpha_raw <- greta::normal(-1.8, # corresponds to roughly 1/7
                                   0.25,
                                   truncation = c(0, Inf),
                                   dim = c(6, 1))

    # state effect
    dow_alpha_state_raw <- greta::normal(0, # corresponds to roughly 1/7
                                         1,
                                         dim = c(n_jurisdictions, 6))

    dow_alpha_sd <- greta::normal(0, 1, truncation = c(0, Inf))

    # combine
    dow_alpha_state_scaled <- dow_alpha_state_raw*dow_alpha_sd
    dow_alpha <- greta::sweep(dow_alpha_state_scaled,
                              2, dow_alpha_raw, "+")
  }

  dow_dist <- greta::imultilogit(dow_alpha)

  dow_effect <- t(dow_dist * 7)

  greta_arrays <- module(
    dow_dist)

  return(greta_arrays)
}
