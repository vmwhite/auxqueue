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
Solve_K <- function(s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p, N = s+1, Stability) {
if (Stability == TRUE){
  N <- max(s+1, 11, N) # K > s
  L_p_N_val <- 1
  L_a_N_val <- 1
  L_P_N_1 <- .0001
  L_A_N_1 <-.0001

  while ( L_a_N_val > 0.02 || L_p_N_val> 0.02  ){
    N <- N + 1
    A <- Calc_Am(N,s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p)
    R <- Calc_R(A, N, s)
    reduced_by = 0
    count = 0
    while (R == FALSE){
      count =count +1
      reduced_by <- count^2
      N <- N - s
      s <- round(s/reduced_by)
      N <- N + s
      r <- round(r /reduced_by)
      lambda  <- lambda /reduced_by
      lambda_aux <- lambda *(1-p)
      lambda_p  <- lambda*p
      A <- Calc_Am(N,s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p)
      R <- Calc_R(A, N, s)
    }
    B <- Calc_Bmn(N,s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p)
    X <- Calc_X(N,s,r,A,B,R)

    # calc L_P
    L_P_N <- L_P_N_1
    L_P_N_1 <- 0
    for (i in 0:(nrow(X)-1)){
      for (j in 0:N){
        L_P_N_1 <- L_P_N_1 + i*X[i+1,j+1] #since R indexes at 1
      }
    }
    L_p_N_val <- abs(L_P_N_1  - L_P_N) / L_P_N

    #Calc L_A
    L_A_N <- L_A_N_1
    L_A_N_1 <- 0
    for (i in (0:(nrow(X)-1))){
      for (j in 0:N){
        L_A_N_1 <- L_A_N_1 + j*X[i+1,j+1] #since R indexes at 1
      }
    }
    L_a_N_val <- abs(L_A_N_1  - L_A_N) / L_A_N

    print(paste0("N = ", N, ", L_p_diff =", L_p_N_val, ", L_a_diff = ", L_a_N_val))
  }
  K <- N
}else{
  A <- Calc_Am(N,s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p)
  R <- Calc_R(A, N, s)
  reduced_by = 0
  count = 0
  while (length(R) == 1){
    count =count +1
    reduced_by <- count^2
    N <- N -s
    s <- round(s/reduced_by)
    N <- N +s
    r <- round(r /reduced_by)
    lambda  <- lambda /reduced_by
    lambda_aux <- lambda *(1-p)
    lambda_p  <- lambda*p
    A <- Calc_Am(N,s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p)
    R <- Calc_R(A, N, s)
  }
  B <- Calc_Bmn(N,s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p)
  X <- Calc_X(N,s,r,A,B,R)
}
  newList <- list("X" = X, "reduced_by" = reduced_by)
return(newList)
}
