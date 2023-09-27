#' matrix intermediary function
#'
#' @param A := Matrix to raise to the n power
#' @param n := power to raise A matrix to
#'
#' @return completed matrix operation
#' @export
#'
#' @examples
#' K <- 13
#' r <- 10
#' s <- 5
#' lambda <- .2
#' p <-.9
#' lambda_aux <- (p-1) * lambda
#' lambda_p <- p * lambda
#' mu_p <-.4
#' mu_aux <-.5
#' A <- Calc_Am(K,s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p)
#' n <- 3
#' matrix_power(A, n)
matrix_power <- function(A, n) {
  count = 0
  while(count <= n){
    if (count == 0){
      A_power <- diag(nrow(A))
    }else if (count == 1){
      A_power <- A
    }else{
      A_power <- A_power %*% A
    }
    count <- count + 1
    }
  return(A_power)
}
