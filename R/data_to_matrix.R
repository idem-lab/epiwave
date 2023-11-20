#' Turn wide data to matrix
#'
#' @param long_data long data form
#' @param extra_col for column name
#'
#' @return matrix form
data_to_matrix <- function (long_data, extra_col) {

  ## make jurisdiction an argument option
  ## and test with only one

  wide_data <- long_data |>
    dplyr::group_by(date, jurisdictions) |>
    dplyr::summarise(extra_col = dplyr::n()) |>
    tidyr::pivot_wider(id_cols = date,
                       names_from = jurisdictions,
                       values_from = !!extra_col,
                       values_fn = sum,
                       values_fill = 0) |>

  # if date gap exists, fill

    fill_date_gaps() |>
    tibble::column_to_rownames(var = date) |>
    data.frame()

  wide_data[is.na(wide_data)] <- 0

  wide_data
}

#' Title
#'
#' @param df
#'
#' @return
fill_date_gaps <- function (df) {

  date_sequence <- data.frame(
    date = seq(min(df$date),
               max(df$date),
               by = "days"))
  df_with_all_dates <- merge(df, date_sequence,
                             by = 'date',
                             all.x = FALSE, all.y = TRUE)
  df_with_all_dates
}
