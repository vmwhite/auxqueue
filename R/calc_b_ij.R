#b_ij function to calculate A matrix
#' Title
#'
#' @param i
#' @param j
#' @param mu_p
#' @param p
#' @param mu_aux
#'
#' @return
#' @export
#'
#' @examples
calc_b_ij<- function(i, j, mu_p, p, mu_aux) {
  q <- 1-p
  b <- (-i * mu_p * q) - (j * mu_aux)
  return(b)
}
