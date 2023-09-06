#' alpha_i intermediate function to aid in calculating A_m matrix
#'
#' @param p := percentage of customers that are of type primary
#' @param i := index
#'
#' @return alpha
#' @export
#'
#' @examples
#' p <- .8
#' i <- 3
#' calc_alpha_i(p, i)
calc_alpha_i <- function(p, i) {
  q <- 1-p
  alpha <- p*q^i
  return(alpha)
}
