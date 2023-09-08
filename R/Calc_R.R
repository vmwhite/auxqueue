###### solve for R
#' Title
#'
#' @param A := A_m matrix that represents non-boundary behavior
#' @param K := truncation parameter
#' @param s := number of auxiliary servers
#'
#' @return R matrix
#' @export
#'
#' @examples
#' K <- 13
#' s <- 5
#' r <- 10
#' lambda <- .2
#' p <-.9
#' lambda_aux <- (p-1) * lambda
#' lambda_p <- p * lambda
#' mu_p <-.4
#' mu_aux <-.5
#' A <- Calc_Am (K,s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p)
#' Calc_R(A, K,s)
Calc_R <- function(A, K,s){
  matrix_size <- (K) + 1
  # A go from 0 to K-s+2
  A_m_lim <- K - s + 2 + 1 #add one since R indexes at 1
  #I is the identity matrix
  I <- diag(1,nrow=matrix_size)
  if (A_m_lim < 2 ){
    Z <- diag(0,nrow=matrix_size)
    temp <-  Z - A[1,,]
    temp <-  (I - temp)
  }else{
    temp <- (I - A[2,,])
  }
  temp <- solve(temp)
  R_Nminusone <- matrix(0, nrow=matrix_size, ncol=matrix_size)
  K_val = K-s+2 +1
  #time_limit <- 10
  #start_time <- Sys.time()
  #while (difftime(Sys.time(), start_time, units = "secs") < time_limit ) {
  for (N in 1:(A_m_lim)){
      R_N <-matrix(0, nrow=matrix_size, ncol=matrix_size)
      R_N <- R_N + A[1,,]
      if (A_m_lim > 3){
        for (m in 3:(A_m_lim)){
          R_N <- R_N + (matrix_power(R_Nminusone,m-1)%*% A[m,,])
        }
        R_N <- R_N %*% temp
      } else if (A_m_lim == 3){
        R_N <- R_N + (matrix_power(R_Nminusone,m-1)%*% A[m,,])
        R_N <- R_N %*% temp
      }
      R_Nminusone <- R_N
  }
  R<- R_N
  return(R)
}
