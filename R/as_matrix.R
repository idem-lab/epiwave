#' Turn wide data to matrix
#'
#' @description Take long form input data and turn into wide form matrix with
#'  continuous sequence of dates. This function assumes that users supply long
#'  form count data where dates with 0 cases are explicit.
#'
#' @param long_data long data form
#' @param ... extra args
#'
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#'
#' @return data in wide matrix form, with continuous seq of dates
#'
#' @export
as_matrix <- function(long_data, ...) {
  UseMethod("as_matrix")
}

#' @export
as_matrix.numeric <- function (value, ...) {
  value
}

#' @export
as_matrix.epiwave_timeseries <- function (long_data, ...) {

  keep_df <- as.data.frame(long_data[c('date', 'jurisdiction', 'value')])
  wide_data <- keep_df %>%
    tidyr::pivot_wider(id_cols = date,
                       names_from = jurisdiction,
                       values_from = value,
                       values_fill = NA)
  wide_all_dates <- tibble::tibble(date = fill_date_gaps(wide_data)) |>
    dplyr::left_join(wide_data) |>
    tibble::column_to_rownames(var = 'date') %>%
    as.matrix()

  wide_all_dates
}

#' @export
as_matrix.epiwave_greta_timeseries <- function (prop, ...) {

  dates <- fill_date_gaps(prop$timeseries$date)

  jurisdictions <- unique(prop$timeseries$jurisdiction)
  ihr <- prop$ihr

  dim(ihr) <- c(length(unique(dates)), length(unique(jurisdictions)))

  ihr
}

#' Fill date gaps
#'
#' @description Fill in date gaps in data, if they exist, so that there is a
#'  continuous sequence of dates in matrix that is fed to the model.
#'
#' @param df Wide dataframe
#'
#' @importFrom methods is
#'
#' @keywords internal
#' @return Wide dataframe with rows filled in so it has continuous seq of dates

fill_date_gaps <- function (df) {

  if (!methods::is(df$date, 'Date')) {
    df$date <- as.Date(df$date)
  }
  dates <- seq(min(df$date),
                           max(df$date),
                           by = "days")

  dates
}

