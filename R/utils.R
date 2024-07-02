#' Calculate a convolution matrix
#'
#' @description Get a matrix to use for forward convolution, given the mass
#'  function(s) (integrating to 1 for convolution) and the number of timepoints
#'  to convolve. This function can take either a single mass function, or a
#'  vector of mass functions with length = n. The latter parameterisation is
#'  necessary when the delay mass function is time varying. Note that we need
#'  to think about if we should make explicit if the user is inputting a single
#'  mass function or a list of mass functions.
#'
#' @param delay_distribution single function or vector of functions that describe
#'  probability mass across delays
#' @param jurisdiction default no jurisdiction specified (eg for gi or one location)
#' @param n dimensions of output convolution matrix
#'
#' @return matrix for forward convolution
get_convolution_matrix <- function (delay_distribution,
                                    jurisdiction = NULL,
                                    n) {

  # if (!length(mass_functions) %in% c(1, n)) {
  #   stop('number of mass functions do not match the number of timepoints!')
  # }

  if(is.null(jurisdiction)) subset_distribution <- delay_distribution
  if(!is.null(jurisdiction)) {
    subset_distribution <- delay_distribution[
      delay_distribution$jurisdiction == jurisdiction, ]
  }

  # get a matrix of time differences between pairs of days
  day_diff <- matrix(NA, n, n)
  day_diff <- row(day_diff) - col(day_diff)

  con_mat <- evaluate(subset_distribution, day_diff)

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
