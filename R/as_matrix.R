#' Turn long form data into a date-aligned vector
#'
#' @description Take long form input data and align it to a shared master
#'  date sequence (`target_infection_dates`), returning a numeric vector with
#'  one entry per date. This function assumes that users supply long form
#'  count data where dates with 0 cases are explicit; dates in
#'  `target_infection_dates` without a corresponding row in `data` are
#'  returned as `NA`.
#'
#' @param data long data form
#' @param target_infection_dates full date sequence to align the data to
#' @param ... extra args
#'
#' @return a numeric vector of length `length(target_infection_dates)`,
#'  named by date
#'
#' @export
as_matrix <- function(data, target_infection_dates, ...) {
  UseMethod("as_matrix")
}

#' @export
as_matrix.numeric <- function (data, target_infection_dates, ...) {

  n <- length(target_infection_dates)
  if (!(length(data) %in% c(1, n))) {
    stop('`data` must have length 1 or length(target_infection_dates)')
  }

  out <- rep(data, length.out = n)
  names(out) <- as.character(target_infection_dates)
  out
}

#' @export
as_matrix.epiwave_timeseries <- function (data, target_infection_dates, ...) {

  data_dates <- as.Date(data$date)
  target_dates <- as.Date(target_infection_dates)

  idx <- match(data_dates, target_dates)
  if (anyNA(idx)) {
    warning('Some dates in `data` fall outside `target_infection_dates` ',
            'and will be dropped.')
  }

  out <- rep(NA_real_, length(target_dates))
  out[idx[!is.na(idx)]] <- data$value[!is.na(idx)]
  names(out) <- as.character(target_dates)
  out
}
