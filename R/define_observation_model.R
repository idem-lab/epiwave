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

  # NEW COMMENTS
  # maybe here can allow for indexing of infection dates vs. case dates for informing inits
  # move creation of inits here

  # inits
  inits_by_jurisdiction <- function (n_juris_ID, cases, delays) {

    cases_by_juris <- cases[, n_juris_ID]
    incidence_approx_unshifted <- cases_by_juris / 0.75
    avg_delay <- mean(unlist(
      lapply(delays, function (x) round(sum(x$delay * x$mass)))
    ))
    shift_index <- pmin(avg_delay + seq_along(cases_by_juris), length(cases_by_juris))
    incidence_approx <- incidence_approx_unshifted[shift_index]

    return(incidence_approx)
  }


  inits_data_mat <- observation_list[[data_inform_inits]]$case_mat
  inits_delays <- observation_list[[data_inform_inits]]$delays
  n_jurisdictions <- length(target_jurisdictions)
  inits_list <- lapply(1:n_jurisdictions, inits_by_jurisdiction,
                       inits_data_mat, inits_delays$value)
  inits_df <- do.call(cbind, inits_list)

  inits_dates <- as.Date(rownames(inits_data_mat))
  inits_idx <- target_infection_dates %in% inits_dates



  observations <- list(observation_model_data = observation_list,
                       target_infection_dates = target_infection_dates,
                       target_jurisdictions = target_jurisdictions,
                       inits_df = inits_df)

# check these are valid objects
  # observation_list <- lapply(observation_list,
  #                            check_valid_observation_object)

  # sanitise dates across all observations in the list

  # check they have the same temporal resolution if needed?

  # define the variable and fixed proportions to ensure identifiability
  ## this will be important.
  ## if prop_mat is a greta array not actual values, then this should be relative to each other

}
