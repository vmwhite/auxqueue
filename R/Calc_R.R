###### solve for R
#' Title
#'
#' @param A := A_m matrix that represents non-boundary behavior
#' @param s
#' @param tolerance
#' @param K := truncation parameter
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
calc_R <- function(A, K,s, tolerance = 0.0001){
  matrix_size <- (K) + 1
  # A go from 0 to K-s+2
  A_m_lim <- K - s + 2 + 1 #add one since R indexes at 1
  #I is the identity matrix
  I <- diag(1,nrow=matrix_size)
  temp <- (I - A[2,,])
  temp <- solve(temp)

  R_Nminusone <- matrix(0, nrow=matrix_size, ncol=matrix_size)
  #initial R_0 GUESS
  R_N <- A[1,,] %*%  temp

  #update Guess
  while ( max(abs(diff(R_N - R_Nminusone))) > tolerance){
    tol <- max(abs(diff(R_N - R_Nminusone)))
    # last guess
    R_Nminusone <- R_N
    #R^0 * A_0
    R_N <- A[1,,]
    # sum m = 2:A_m_lim R^k A_k
    m = 2+1 # since R indexes at 1
    while (m <= A_m_lim){
      R_N <- R_N + (matrix_power(R_Nminusone,m-1) %*% A[m,,])
      m <- m +1
    }
    #multiply sum by determinate of (I - A_1)
    R_N <- R_N %*% temp
    R <- R_N
    # test if resulting R matrix is inversable
    t <- try(solve(I - R_N))
    if("try-error" %in% class(t)){
      R_N <-  R_Nminusone
      R <- R_Nminusone
      R <- FALSE
      print(paste0("Reducing size by 2" ))
      #print(paste0("Using prior R approximation instead. current difference in approximation iterations:", tol))
      break
  }
}
  return(R)
}
