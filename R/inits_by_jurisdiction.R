#' Define initial values for the infection timeseries
#'
#' @description Called internally by `compute_flat_prior_inits()`, once per
#'  stream per jurisdiction.
#'
#' @param obs_data numeric vector of observed data, aligned to
#'   `target_infection_dates` (`NA` for dates without observations)
#' @param delays a `discrete_pmf_series` delay object, aligned to
#'   `target_infection_dates`
#' @param obs_prop numeric vector of proportions, aligned to
#'   `target_infection_dates`, or a greta array of the same
#' @param target_infection_dates infection date sequence
#'
#' @importFrom mgcv gam predict.gam
#'
#' @return initial values as a list
#' @noRd
inits_by_jurisdiction <- function (obs_data,
                                   delays,
                                   obs_prop,
                                   target_infection_dates) {

  # currently creates inits to cover observable idx.
  if (inherits(obs_prop, 'greta_array')) {
    obs_prop_sim <- greta::calculate(obs_prop, nsim = 100)
    obs_prop <- as.numeric(apply(obs_prop_sim$obs_prop, 2:3, mean))
  }

  observed <- !is.na(obs_data)
  cases_by_juris <- obs_data[observed]
  inits_prop_by_juris <- obs_prop[observed]
  infection_approx <- cases_by_juris / inits_prop_by_juris

  # dates observable for this data type
  case_dates <- as.Date(target_infection_dates)[observed]

  # NOTE here we are using forward delay distribution as if it's a backward delay for purpose of imputation
  # since we are only taking the mean and only specifying inits, it doesn't affect model
  # if we wanted to improve the inits, we'd have to reverse the delay dist.

  # if delays are time-varying it is still only using the average delay
  # to shift observation data for calculation of inits
  delays_juris <- delays[case_dates]

  expected_delay_vals <- unlist(lapply(delays_juris$values, function (x)
    round(mean(x))
  ))

  max_delay_vals <- unlist(lapply(delays_juris$values, function (x)
    max(x$step)
  ))
  df <- data.frame(max_delay_vals, case_dates) # consider max_delay_vals + 1

  all_dates_all_delays <- lapply(seq_len(nrow(df)), function(x) df$case_dates[x] - 0:df$max_delay_vals[x])
  observable_dates <- sort(unique(do.call(c, all_dates_all_delays)))

  # inputing infection dates by average delay
  inferred_infection_dates <- case_dates - expected_delay_vals
  inferred_infection_idx <- match(inferred_infection_dates, observable_dates)

  smooth_approx <- mgcv::gam(
    val ~ s(idx),
    data = data.frame(val = infection_approx, idx = inferred_infection_idx),
    family = mgcv::nb(link = "log")
  )

  smooth_pred <- mgcv::predict.gam(
    smooth_approx,
    newdata = data.frame(idx = seq_along(observable_dates)),
    type = "response"
  )

  inits_values <- smooth_pred + 1

  observable_idx <- match(observable_dates, target_infection_dates)

  # drop observable dates that fall outside target_infection_dates (e.g. a
  # delay-shifted date walking off the front edge of the window)
  keep <- !is.na(observable_idx)
  observable_idx <- observable_idx[keep]
  inits_values <- inits_values[keep]

  out <- list(inits_values = inits_values,
              observable_idx = observable_idx)

  return(out)
}
