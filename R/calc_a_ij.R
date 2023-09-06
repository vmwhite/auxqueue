#' a_ij intermediate function to aid in calculating A matrix
#'
#' @param i := first subscript
#' @param j := second subscript
#' @param mu_p := primary server service rate
#' @param mu_aux := auxiliary server service rate
#' @param p := percentage of customers that are of type primary
#'
#' @return a_ij constant
#' @export
#'
#' @examples
#' i <- 1
#' j <- 2
#' mu_p <-.4
#' mu_aux <-.5
#' p <-.9
#' calc_a_ij(i, j, mu_p, mu_aux, p)
calc_a_ij <- function(i, j, mu_p, mu_aux, p) {
  q <- 1-p
  a <- (-i * mu_p * q) - (j * mu_aux * p)
  return(a)
}
