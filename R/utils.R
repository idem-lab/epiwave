#' Calculate a convolution matrix
#'
#' @description Get a matrix to use for forward convolution, given the mass function(s) (integrating to 1 for
#' convolution) and the number of timepoints to convolve. This function can take either a single
#' mass function, or a list of mass functions with length = n. The latter parametrisation is
#' necessary when the delay mass function is time varying. Note that we need to think about if we
#' should make explicit if the user is inputting a single mass function or a list of mass functions
#'
#' @param mass_functions
#' @param n dimensions of output convolution matrix
#'
#' @return
get_convolution_matrix <- function (mass_functions, n) {

  if (!length(mass_functions) %in% c(1, n)) {
    stop('number of mass functions do not match the number of timepoints!')
  }

  # get a matrix of time differences between pairs of days
  day_diff <- matrix(NA, n, n)
  day_diff <- row(day_diff) - col(day_diff)

  if (length(mass_functions) == 1) {

    # apply the single mass function
    message("using a single mass function for all time points!")
    con_mat <- matrix(mass_functions(day_diff), n, n)
  }

  if (length(mass_functions) == n) {

    # apply the matching mass functions to these delays
    con_mat <- matrix(as.numeric(
      mapply(function(f, x) do.call(f, list(x)),
             mass_functions,
             day_diff)), n, n)
  }

  return(con_mat)
}

#' Create module of greta arrays
#'
#' @description Sets up a named list of greta arrays.
#'
#' @keywords internal
#' @return Named list
module <- function (..., sort = TRUE) {

  dots <- list(...)
  names <- names(dots)
  cl <- match.call()
  nm <- as.character(as.list(cl)[-1])
  if (is.null(names)) {
    names(dots) <- nm
  }
  else {
    blank_names <- names == ""
    names[blank_names] <- nm[blank_names]
    names(dots) <- names
  }
  if (sort) {
    dots <- dots[order(names(dots))]
  }
  dots

}
#' Turn wide data to matrix
#'
#' @param long_data long data form
#' @param extra_col for column name
#'
#' @keywords internal
#' @return matrix form
data_to_matrix <- function (long_data, extra_col) {

  wide_data <- long_data |>
    tidyr::pivot_wider(id_cols = date,
                       names_from = jurisdiction,
                       values_from = !!extra_col,
                       values_fill = NA) |> # assuming rows with 0 cases are explicit
    fill_date_gaps() |>
    tibble::column_to_rownames(var = 'date') |>
    as.matrix()

  wide_data
}

#' Title
#'
#' @param df
#'
#' @keywords internal
#' @return
fill_date_gaps <- function (df) {

  if (!is(df$date, 'Date')) df$date <- as.Date(df$date)
  date_sequence <- data.frame(
    date = seq(min(df$date),
               max(df$date),
               by = "days"))
  df_with_all_dates <- merge(df, date_sequence,
                             by = 'date',
                             all.x = FALSE, all.y = TRUE)
  df_with_all_dates
}

