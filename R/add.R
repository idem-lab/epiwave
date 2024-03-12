#' Default combine function for lowerGPreff_massfun objects
#'
#' @param delay_massfun1 first massfun object
#' @param delay_massfun2 second massfun object
#' @param ... extra args
#'
#' @return combined lowerGPreff massfuns
#'
#' @export
add <- function(delay_massfun1, delay_massfun2, ...) {
  UseMethod("add")
}

#' @export
add.lowerGPreff_massfun <- function (delay_massfun1,
                                     delay_massfun2, ...) {

  # add check that both inputs are same class (ie both distribution or
  # both curve)

  names(delay_massfun1) <- c('delay1', 'massfun1')
  names(delay_massfun2) <- c('delay2', 'massfun2')

  p <- tidyr::expand_grid(
    delay1 = delay_massfun1$delay1, # ie. incubation period
    delay2 = delay_massfun2$delay2) %>% # ie. sym to notif
    dplyr::left_join(delay_massfun1) %>%
    dplyr::left_join(delay_massfun2) %>%
    dplyr::mutate(delay = delay1 + delay2) %>%
    dplyr::group_by(delay) %>%
    dplyr::summarise(mass =
                       sum(massfun1 * massfun2))

  # is there a way for p to inherit class of delay_massfun1 or 2?
  class(p) <- c(class(delay_massfun1), class(p))
  p

}
