#### Solve for Truncation parameter K ######
#' Title
#'
#' @param s
#' @param r
#' @param lambda
#' @param lambda_aux
#' @param lambda_p
#' @param mu_p
#' @param mu_aux
#' @param p
#' @param N
#'
#' @return
#' @export
#'
#' @examples
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
return(c(K, X))
}
