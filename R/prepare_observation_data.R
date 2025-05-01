#' Prepare observation data
#'
#' @description Prepare list of data objects needed for observation model.
#'
#' @param observation_data long form data of counts of notifications
#' @param target_infection_dates full date sequence of infection timeseries
#' @param target_jurisdictions jurisdictions
#'
#' @importFrom greta %*% as_data negative_binomial normal sweep zeros
#' @importFrom epiwave.params create_epiwave_massfun_timeseries create_epiwave_timeseries
#'
#' @return greta arrays of observation model
#' @export
prepare_observation_data <- function (observation_data,
                                      target_infection_dates,
                                      target_jurisdictions) {


  # check the data is a valid timeseries object (or can be coerced to one) and
  # that it contains data consistent with the observation model type

  # add if statements to check that infection_days is long enough to cover
  # the right period

  # can check for juris and dates for all input data

  delays <- observation_data$delay_from_infection
  if (!('epiwave_massfun_timeseries' %in% class(delays))) {
    delays <- epiwave.params::create_epiwave_massfun_timeseries(
      dates = target_infection_dates,
      jurisdictions = target_jurisdictions,
      value = delays)
  }

  prop <- observation_data$proportion_infections
  if (!('epiwave_timeseries' %in% class(prop))) {
    prop <- epiwave.params::create_epiwave_timeseries(
      dates = target_infection_dates,
      jurisdictions = target_jurisdictions,
      value = prop)
  }

  case_mat <- as_matrix(observation_data$timeseries_data)

  prop_mat <- as_matrix(prop)

  # check jurisdictions
  jurisdictions_cases_mismatch <- !setequal(colnames(case_mat), target_jurisdictions)
  # jurisdictions_prop_mismatch <- !setequal(colnames(prop_mat), target_jurisdictions)
  if (jurisdictions_cases_mismatch){ # | jurisdictions_prop_mismatch) {
    stop('Jurisdictions in data do not match specified jurisdictions')
  }

  # reorder columns if necessary
  case_mat <- as.matrix(as.data.frame(case_mat[, target_jurisdictions]))
  # prop_mat <- as.matrix(as.data.frame(prop_mat[, target_jurisdictions]))

  # check dates
  # prop_dates <- as.Date(rownames(prop_mat))
  #
  # if (!identical(prop_dates, target_infection_dates)) {
  #   stop('dates supplied in proportion_infections should match target_infection_dates')
  # }


  ## we want to check the case dates are within the infection sequence, and also
  ## that with delays it covers the right period
  # if (case_dates ....)

    # check the data is a valid timeseries object (or can be coerced to one) and
  # that it contains data consistent with the observation model type

  # add if statements to check that infection_days is long enough to cover
  # the right period





  n_jurisdictions <- length(target_jurisdictions)
  inits_list <- lapply(1:n_jurisdictions,
                       inits_by_jurisdiction,
                       case_mat, # obs_mat
                       delays,
                       obs_prop = prop_mat,
                       target_infection_dates)

  inits_values_mat <- matrix(NA,
                             length(target_infection_dates),
                             n_jurisdictions)
  observable_idx_mat <- matrix(FALSE,
                               length(target_infection_dates),
                               n_jurisdictions)

  for (juris in 1:n_jurisdictions) {

    idx <- inits_list[[juris]]$observable_idx
    inits <- inits_list[[juris]]$inits_values
    inits_values_mat[idx, juris] <- inits
    observable_idx_mat[idx, juris] <- TRUE

  }

  if (observation_data$dow_model) {
    dow_priors <- create_dow_priors(length(target_jurisdictions))
    dow_correction <- implement_day_of_week(target_infection_dates, dow_priors)
    prop_mat <- prop_mat * dow_correction
  }

  out <- list(timeseries_data = observation_data$timeseries_data,
              delays = delays,
              case_mat = case_mat,
              prop_mat = prop_mat,
              inform_inits = observation_data$inform_inits,
              inits_values_mat = inits_values_mat,
              observable_idx_mat = observable_idx_mat)

  out

}
