#' Evaluate the mass for given delay
#'
#' @param lookup epiwave_distribution object with mass functions
#' @param day_diff matrix of
#' @param ... extra args
#'
#' @return matrix of delay days with mass vals filled in
#'
#' @export
evaluate <- function(lookup,
                     day_diff, ...) {
  UseMethod("evaluate")
}

#' @export
evaluate.epiwave_distribution_massfun <- function (lookup,
                                                   day_diff,
                                                   ...) {

  lookup <- as.data.frame(lookup)
  day_diff[] <- lookup$mass[
    match(unlist(day_diff),
          lookup$delay)]
  day_diff[is.na(day_diff)] <- 0

  class(day_diff) <- 'numeric' # check and correct
  day_diff

}

#' @export
evaluate.epiwave_massfun_timeseries <- function (lookup,
                                                 day_diff,
                                                 ...) {

  con_list <- lapply(
    1:ncol(day_diff),
    function (x) {
      distribution <- lookup[x, 'value'][[1]]
      evaluate(distribution, day_diff[, x])
    }
  )
  con_mat <- do.call(cbind, con_list)
  con_mat

}
