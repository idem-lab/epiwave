#' Define observation model
#'
#' @description
#' Placeholder function for now to create list of the observation datasets.
#'
#' @param target_infection_dates sequence of infection dates
#' @param target_jurisdictions jurisdictions
#' @param data_inform_inits which dataset to inform inits with
#' @param x observation data set
#' @param ... optional additional observation data sets
#'
#' @return list of datasets
#' @export
#'
define_observation_model <- function (target_infection_dates = NULL,
                                      target_jurisdictions,
                                      data_inform_inits,
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
#
#   # NEW COMMENTS
#   # maybe here can allow for indexing of infection dates vs. case dates for informing inits
#   # move creation of inits here
#
#   inits_data_mat <- prepared_observation_model_data[[data_inform_inits]]$case_mat
#   inits_delays <- prepared_observation_model_data[[data_inform_inits]]$delays
#   inits_prop <- prepared_observation_model_data[[data_inform_inits]]$inits_prop_vec
#   n_jurisdictions <- length(target_jurisdictions)
#   inits_list <- lapply(1:n_jurisdictions,
#                        inits_by_jurisdiction,
#                        inits_data_mat,
#                        inits_delays$value,
#                        obs_prop = inits_prop,
#                        target_infection_dates)
#
#   inits_values_mat <- as.matrix(do.call(cbind, inits_list))
  # inits_values_list <- lapply(inits_list,
  #                             function(x) x$inits_values)
  # inits_values_mat <- as.matrix(do.call(cbind, inits_values_list))

  # inits_idx_list <- lapply(inits_list,
  #                          function(x) x$inits_idx)
  # inits_idx_mat <- as.matrix(do.call(cbind, inits_idx_list))

  observations <- list(observation_model_data = prepared_observation_model_data,
                       target_infection_dates = target_infection_dates,
                       target_jurisdictions = target_jurisdictions#,
                       #inits_values_mat = inits_values_mat
                       )#,
                       # inits_idx_mat = inits_idx_mat)

# check these are valid objects
  # observation_list <- lapply(observation_list,
  #                            check_valid_observation_object)

  # sanitise dates across all observations in the list

  # check they have the same temporal resolution if needed?

  # define the variable and fixed proportions to ensure identifiability
  ## this will be important.
  ## if prop_mat is a greta array not actual values, then this should be relative to each other

}
