#' Solve for Truncation parameter K ######
#' @param s := number of auxiliary servers
#' @param r := number of primary servers
#' @param lambda := arrival rate of all customers
#' @param lambda_aux := arrival rate of customers of type auxiliary
#' @param lambda_p := arrival rate of customers of type primary
#' @param mu_p := primary server rate
#' @param mu_aux := auxiliary server rate
#' @param p := percentage of customers that are of type primary
#' @param N := temporary starting truncation parameter
#'
#' @return final K truncation parameter with key metric stability
#' @export
#'
#' @examples
#' r <- 10
#' s <- 5
#' lambda <- .2
#' p <-.9
#' lambda_aux <- (p-1) * lambda
#' lambda_p <- p * lambda
#' mu_p <-.4
#' mu_aux <-.5
#' K <- Solve_K(s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p, N = s+1)
Solve_K <- function(s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p, N = s+1) {
N <- max(s+1, 11) # K > s
L_p_N_val <- 1
L_a_N_val <- 1
L_P_N_1 <- .0001
L_A_N_1 <-.0001
while ( L_a_N_val > 0.02 || L_p_N_val> 0.02  ){
  N <- N + 1
  A <- Calc_Am(N,s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p)
  B <- Calc_Bmn(N,s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p)
  R <- Calc_R(A, N, s)
  X <- Calc_X(N,s,r,A,B,R)

  # calc L_P
  L_P_N <- L_P_N_1
  L_P_N_1 <- 0
  for (i in 1:nrow(X)){
    for (j in 1:(N+1)){
      L_P_N_1 <- L_P_N_1 + (i -1)*X[i,j]
    }
  }
  L_p_N_val <- abs(L_P_N_1  - L_P_N) / L_P_N

  #Calc L_A
  L_A_N <- L_A_N_1
  L_A_N_1 <- 0
  for (i in 1:nrow(X)){
    for (j in 1:(N+1)){
      L_A_N_1 <- L_A_N_1 + (j -1)*X[i,j]
    }
  }
  L_a_N_val <- abs(L_A_N_1  - L_A_N) / L_A_N

  print(paste0("N = ", N, ", L_p_diff =", L_p_N_val, ", L_a_diff = ", L_a_N_val))
}
K <- N
return(X)
}
