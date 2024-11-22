#' Define initials values by jurisdiction
#'
#' @param n_juris_ID numbered jurisdiction index
#' @param infection_inits_data matrix form of data to inform infection initials
#' @param delays epiwave timeseries delay object
#' @param obs_prop single numeric value of proportion
#' @param target_infection_dates infection date sequence
#' @param smooth logical whether to apply smoothing function of data to inform
#'    initials
#'
#' @return initial values as matrix
#' @export
inits_by_jurisdiction <- function (n_juris_ID,
                                   infection_inits_data,
                                   delays,
                                   obs_prop,
                                   target_infection_dates,
                                   smooth = FALSE) {

  case_dates <- as.Date(rownames(cases))
  case_idx <- which(target_infection_dates %in% case_dates)

  cases_by_juris <- cases[, n_juris_ID]
  infection_approx <- cases_by_juris / obs_prop

  avg_delay <- mean(unlist(lapply(delays, function (x)
    round(
      sum(x$delay * x$mass)
    ))))

  inits_idx <- case_idx - avg_delay

  if (smooth) {
    smooth_approx <- mgcv::gam(
      val ~ s(idx) + s(idx, bs = "cc"),
      data = data.frame(val = infection_approx, idx = seq_along(infection_approx)),
      family = mgcv::nb(link = "log")
    )

    smooth_pred <- mgcv::predict.gam(
      smooth_approx,
      newdata = data.frame(idx = seq_along(infection_approx)),
      type = "response"
    )

    infection_approx[] <- smooth_pred
  }

  inits_values <- rep(0, length(target_infection_dates))
  inits_values[inits_idx] <- infection_approx
  inits_values[is.na(inits_values)] <- 0
  inits_values <- pmax(inits_values, .Machine$double.eps)

  return(inits_values)
}
