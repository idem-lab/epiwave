#' Visualise observation data coverage against the derived date axis
#'
#' @description Since `target_infection_dates` is emergent rather than
#'  user-supplied (see `stack_jurisdictions()`), it's worth being able to
#'  see what actually got derived and how each stream/jurisdiction's data
#'  aligns to it before committing to a (potentially slow) `fit_waves()`
#'  run -- this is the visual form of the alignment checks used throughout
#'  the package's own tests.
#'
#' @param observations either an already-stacked `epiwave_stacked_observations`
#'  object (output of `stack_jurisdictions()`), or a single jurisdiction's raw
#'  `epiwave_observation_model` (output of `define_observation_model()`),
#'  which is stacked on its own to derive its axis first
#'
#' @importFrom dplyr bind_rows
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_tile facet_wrap scale_fill_manual
#'  scale_x_date labs
#' @importFrom cowplot theme_cowplot panel_border
#'
#' @return a ggplot object: a date x jurisdiction coverage tile per stream
#' @export
plot_observation_coverage <- function (observations) {

  if (!inherits(observations, 'epiwave_stacked_observations')) {
    observations <- stack_jurisdictions_list(list(observations))
  }

  target_infection_dates <- observations$target_infection_dates

  coverage_df <- dplyr::bind_rows(lapply(
    names(observations$observation_model_data),
    function (stream_id) {
      case_mat <- observations$observation_model_data[[stream_id]]$case_mat
      # go via column names (jurisdiction labels), not positional/vector
      # order, to avoid silently mismatching dates to jurisdictions
      observed_df <- as.data.frame(!is.na(case_mat), check.names = FALSE)
      observed_df$date <- target_infection_dates
      observed_long <- tidyr::pivot_longer(
        observed_df,
        cols = -date,
        names_to = "jurisdiction",
        values_to = "observed")
      observed_long$stream <- stream_id
      observed_long
    }))

  ggplot2::ggplot(
    coverage_df,
    ggplot2::aes(x = date, y = jurisdiction, fill = observed)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~stream, ncol = 1) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "steelblue", `FALSE` = "grey90")) +
    cowplot::theme_cowplot() +
    cowplot::panel_border(remove = TRUE) +
    ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
    ggplot2::labs(x = NULL, y = NULL, fill = "observed")
}

#' Print a summary of a stacked observation model
#'
#' @description Reports the derived `target_infection_dates` axis and, per
#'  stream and jurisdiction, the observed date range and coverage against
#'  that axis -- a quick text-only check of what `stack_jurisdictions()`
#'  actually derived, for when a full `plot_observation_coverage()` isn't
#'  needed.
#'
#' @param x an `epiwave_stacked_observations` object
#' @param ... unused, present for consistency with the `print()` generic
#'
#' @return `x`, invisibly
#' @export
print.epiwave_stacked_observations <- function (x, ...) {

  target_infection_dates <- x$target_infection_dates

  cat("<epiwave_stacked_observations>\n")
  cat(sprintf(
    "target_infection_dates: %s to %s (%d days)\n",
    format(min(target_infection_dates)),
    format(max(target_infection_dates)),
    length(target_infection_dates)))
  cat("jurisdictions:", paste(x$target_jurisdictions, collapse = ", "), "\n")

  for (stream_id in names(x$observation_model_data)) {
    cat("\n", stream_id, ":\n", sep = "")
    case_mat <- x$observation_model_data[[stream_id]]$case_mat
    for (jurisdiction in x$target_jurisdictions) {
      observed <- !is.na(case_mat[, jurisdiction])
      if (any(observed)) {
        observed_dates <- target_infection_dates[observed]
        cat(sprintf(
          "  %s: %s to %s (%d of %d days observed)\n",
          jurisdiction,
          format(min(observed_dates)),
          format(max(observed_dates)),
          sum(observed),
          length(target_infection_dates)))
      } else {
        cat(sprintf("  %s: no observed data\n", jurisdiction))
      }
    }
  }

  invisible(x)
}

#' Print a summary of a single jurisdiction's observation model
#'
#' @description Derives this jurisdiction's own axis (as `fit_waves()` would
#'  for a single-jurisdiction fit) and reports it the same way as
#'  `print.epiwave_stacked_observations()`.
#'
#' @param x an `epiwave_observation_model` object
#' @param ... unused, present for consistency with the `print()` generic
#'
#' @return `x`, invisibly
#' @export
print.epiwave_observation_model <- function (x, ...) {
  print(stack_jurisdictions_list(list(x)))
  invisible(x)
}
