#' Evaluate probability mass for a matrix of delay values
#'
#' Looks up the probability mass values corresponding to a matrix of delay
#' values from a `discrete_pmf` object, returning a numeric matrix of the
#' same dimensions. Delay values not present in the lookup return 0.
#'
#' @param pmf a `discrete_pmf` object
#' @param day_diff a matrix of integer delay values
#' @param ... unused
#'
#' @return a numeric matrix of probability mass values with the same
#'   dimensions as `day_diff`
#'
#' @noRd
evaluate <- function(pmf, day_diff, ...) {
  UseMethod("evaluate")
}

#' @noRd
evaluate.discrete_pmf <- function(pmf, day_diff, ...) {
  pmf <- as.data.frame(pmf)
  day_diff[] <- pmf$prob[
    match(unlist(day_diff), pmf$step)
  ]
  day_diff[is.na(day_diff)] <- 0
  day_diff
}

#' Evaluate probability mass column-wise for a time-varying PMF series
#'
#' @param pmf a `discrete_pmf_series` object
#' @param day_diff a matrix of integer delay values
#' @param ... unused
#'
#' @noRd
evaluate.discrete_pmf_series <- function(pmf, day_diff, ...) {
  con_list <- lapply(
    seq_len(ncol(day_diff)),
    function(col_idx) {
      evaluate(pmf$values[[col_idx]], day_diff[, col_idx])
    }
  )
  do.call(cbind, con_list)
}
