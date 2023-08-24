#' mu_mn function to calculate A matrix
#'
#' @param m
#' @param n
#' @param mu_p
#' @param mu_aux
#'
#' @return
#' @export
#'
#' @examples
calc_mu_mn<- function(m, n, mu_p, mu_aux) {
  mu <- (m* mu_p)+ (n * mu_aux)
  return(mu)
}
