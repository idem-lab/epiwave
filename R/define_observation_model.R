#' Define observation model
#'
#' @description
#' Placeholder function for now to create list of the observation datasets.
#'
#' @param x observation data set
#' @param ... optional additional observation data sets
#'
#' @return list of datasets
#' @export
#'
define_observation_model <- function (target_infection_dates = NULL,
                                      target_jurisdictions,
                                      data_inform_inits = 'cases',
                                      x, ...) {

  observation_list <- list(...)

  inits_data_exist <- data_inform_inits %in% names(observation_list)
  if(!inits_data_exist) {
    stop('data_inform_inits must match name of a dataset in observations')
  }

  prepared_observation_model_data <- lapply(observation_list,
                                            prepare_observation_data,
                                            target_infection_dates,
                                            target_jurisdictions)

  # NEW COMMENTS
  # maybe here can allow for indexing of infection dates vs. case dates for informing inits
  # move creation of inits here

  # inits
  inits_by_jurisdiction <- function (n_juris_ID,
                                     cases,
                                     # can change name later
                                     delays,
                                     obs_prop = 0.75,
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

  inits_data_mat <- prepared_observation_model_data[[data_inform_inits]]$case_mat
  inits_delays <- prepared_observation_model_data[[data_inform_inits]]$delays
  n_jurisdictions <- length(target_jurisdictions)
  inits_list <- lapply(1:n_jurisdictions, inits_by_jurisdiction,
                       inits_data_mat, inits_delays$value,
                       obs_prop = 0.75, target_infection_dates)

  inits_values_mat <- as.matrix(do.call(cbind, inits_list))

  observations <- list(observation_model_data = prepared_observation_model_data,
                       target_infection_dates = target_infection_dates,
                       target_jurisdictions = target_jurisdictions,
                       inits_values_mat = inits_values_mat)

# check these are valid objects
  # observation_list <- lapply(observation_list,
  #                            check_valid_observation_object)

  # sanitise dates across all observations in the list

  # check they have the same temporal resolution if needed?

  # define the variable and fixed proportions to ensure identifiability
  ## this will be important.
  ## if prop_mat is a greta array not actual values, then this should be relative to each other

}
