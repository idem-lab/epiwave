#' Define observation model
#'
#' @description
#' Placeholder function for now to create list of the observation datasets.
#'
#' @param target_infection_dates sequence of infection dates
#' @param target_jurisdictions jurisdictions
#' @param x observation data set
#' @param ... optional additional observation data sets
#'
#' @importFrom abind abind
#'
#' @return list of datasets
#' @export
#'
define_observation_model <- function (target_infection_dates = NULL,
                                      target_jurisdictions,
                                      x, ...) {

  observation_list <- list(...)

  prepared_observation_model_data <- lapply(observation_list,
                                            prepare_observation_data,
                                            target_infection_dates,
                                            target_jurisdictions)
<<<<<<< HEAD
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
=======

  # CONSIDER could also use `any` in place of `pmax` if they are 0,1
  # HERE MAY NEED TO HAVE A USE FOR INITS FLAG
  incidence_observable <- do.call(
    pmax, lapply(prepared_observation_model_data,
                 function(x) x$observable_idx_mat))
  incidence_observable <- array(as.logical(incidence_observable),
                                dim = dim(incidence_observable))

  list_inits_mats <- lapply(prepared_observation_model_data,
                            function(x) x$inits_values_mat)
  abind_args <- c(list_inits_mats, along = 3)
  inits_array <- do.call(abind::abind, abind_args)
  incidence_observable_inits <- apply(inits_array, 1:2, mean, na.rm = TRUE)

  observations <- list(observation_model_data = prepared_observation_model_data,
                       target_infection_dates = target_infection_dates,
                       target_jurisdictions = target_jurisdictions,
                       incidence_observable_inits = incidence_observable_inits, # mean init across data types per juris
                       incidence_observable = incidence_observable) # 1 if any juris has data that day
>>>>>>> origin/main

# check these are valid objects
  # observation_list <- lapply(observation_list,
  #                            check_valid_observation_object)

  # sanitise dates across all observations in the list

  # check they have the same temporal resolution if needed?

  # define the variable and fixed proportions to ensure identifiability
  ## this will be important.
  ## if prop_mat is a greta array not actual values, then this should be relative to each other

}
