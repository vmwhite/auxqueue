#alpha_i function to calculate A matrix
#' Title
#'
#' @param p
#' @param i
#'
#' @return
#' @export
#'
#' @examples
calc_alpha_i <- function(p, i) {
  q <- 1-p
  alpha <- p*q^i
  return(alpha)
}
