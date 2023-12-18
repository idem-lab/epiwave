#' Create infection hospitalisation rate prior
#'
#' @description Create a list of the same structure as case ascertainment rate.
#'  Generate an infection hospitalisation rate as a greta array.
#'
#' @param car case ascertainment rate
#'
#' @importFrom greta uniform
#'
#' @return infection fatality rate
#' @export
create_ihr_prior <- function (car) {
  ihr <- as.list(car)
  ihr_car_ratio <- greta::uniform(0, 1)
  ihr$ratio <- ihr_car_ratio
  ihr
}
