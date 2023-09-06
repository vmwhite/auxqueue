#' mu_mn intermediate function to calculate A matrix
#'
#' @param m :=first index
#' @param n := second index
#' @param mu_p := primary server service rate
#' @param mu_aux := auxiliary server service rate
#'
#' @return value of mu_mn
#' @export
#'
#' @examples
#' m < - 5
#' n <- 2
#' mu_p <- .5
#' mu_aux <- .2
calc_mu_mn<- function(m, n, mu_p, mu_aux) {
  mu <- (m* mu_p)+ (n * mu_aux)
  return(mu)
}
