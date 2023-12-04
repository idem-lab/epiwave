#' Turn wide data to matrix
#'
#' @param long_data long data form
#' @param extra_col for column name
#'
#' @return matrix form
#'
#' @export
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
#' @return
#'
#' @export
fill_date_gaps <- function (df) {

  if (class(df$date) != 'Date') df$date <- as.Date(df$date)
  date_sequence <- data.frame(
    date = seq(min(df$date),
               max(df$date),
               by = "days"))
  df_with_all_dates <- merge(df, date_sequence,
                             by = 'date',
                             all.x = FALSE, all.y = TRUE)
  df_with_all_dates
}
