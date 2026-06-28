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
