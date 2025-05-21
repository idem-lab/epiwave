#' Define initials values by jurisdiction
#'
#' @param n_juris_ID numbered jurisdiction index
#' @param obs_data matrix form of data to inform infection initials
#' @param delays epiwave timeseries delay object
#' @param obs_prop single numeric value of proportion
#' @param target_infection_dates infection date sequence
#'
#' @importFrom mgcv gam predict.gam
#' @importFrom dplyr filter
#'
#' @return initial values as matrix
#' @export
inits_by_jurisdiction <- function (n_juris_ID,
                                   obs_data,
                                   delays,
                                   obs_prop,
                                   target_infection_dates) {

  # currently creates inits to cover observable idx.
  if (inherits(obs_prop, 'greta_array')) {
    obs_prop_sim <- greta::calculate(obs_prop, nsim = 100)
    obs_prop <- apply(obs_prop_sim$obs_prop, 2:3, mean)
  }

  cases_by_juris <- obs_data[, n_juris_ID]
  inits_prop_by_juris <- obs_prop[n_juris_ID]
  infection_approx <- cases_by_juris / inits_prop_by_juris

  # dates observable for this data type
  case_dates <- as.Date(rownames(obs_data))

  # NOTE here we are using forward delay distribution as if it's a backward delay for purpose of imputation
  # since we are only taking the mean and only specifying inits, it doesn't affect model
  # if we wanted to improve the inits, we'd have to reverse the delay dist.

  # if delays are time-varying it is still only using the average delay
  # to shift observation data for calculation of inits
  delays_juris <- delays |>
    dplyr::filter(jurisdiction == delays$jurisdiction[n_juris_ID],
                  date %in% case_dates)

  expected_delay_vals <- unlist(lapply(delays_juris$value, function (x)
    round(
      sum(x$delay * x$mass)
    )))

  max_delay_vals <- unlist(lapply(delays_juris$value, function (x)
    max(x$delay)
  ))
  df <- data.frame(max_delay_vals, case_dates) # consider max_delay_vals + 1

  all_dates_all_delays <- lapply(seq_len(nrow(df)), function(x) df$case_dates[x] - 0:df$max_delay_vals[x])
  observable_dates <- sort(unique(do.call(c, all_dates_all_delays)))

  # inputing infection dates by average delay
  inferred_infection_dates <- case_dates - expected_delay_vals
  inferred_infection_idx <- match(inferred_infection_dates, observable_dates)

  # dates_not_in <- any(!(observable_dates %in% target_infection_dates))
  # if(dates_not_in) stop("target_infection_dates does not cover the entire range of observable dates.")

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

  inits_values <- smooth_pred + 1#infection_approx #
  # len_iv <- length(inits_values)
  # inits_values[(len_iv + 1) : (len_iv + avg_delay)] <- 0
  #
  # # inits_values <- smooth_pred# ???
  # # inits_values <- rep(0, length(observable_idx))
  # # inits_values[] <- smooth_pred
  # # inits_values[is.na(inits_values)] <- 0
  # inits_values <- pmax(inits_values, .Machine$double.eps)
  #
  #
  observable_idx <- match(observable_dates, target_infection_dates)


  out <- list(inits_values = inits_values,
              observable_idx = observable_idx)


  return(out)
}
