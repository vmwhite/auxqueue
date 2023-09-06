#b_ij function to calculate A matrix
#' Title
#'
#' @param i := first subscript
#' @param j := second subscript
#' @param mu_p := primary server service rate
#' @param p := percentage of customers that are of type primary
#' @param mu_aux := auxiliary server service rate
#'
#' @return matrix b_ij
#' @export
#'
#' @examples
#' i < - 5
#' j <- 2
#' mu_p <- .5
#' p <- .8
#' mu_aux <- .2
#' calc_b_ij<- function(i, j, mu_p, p, mu_aux)
calc_b_ij<- function(i, j, mu_p, p, mu_aux) {
  q <- 1-p
  b <- (-i * mu_p * q) - (j * mu_aux)
  return(b)
}
