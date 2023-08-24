#' Title
#'
#' @param c_p
#' @param c_aux
#' @param mu_p
#' @param mu_aux
#' @param lambda
#' @param p
#' @param K
#'
#' @return
#' @export
#'
#' @examples
auxqueue <- function(c_p,c_aux,mu_p, mu_aux, lambda, p, K=13){
  #### redefine parameters
  q <- 1-p
  lambda_aux <- q*lambda
  lambda_p <-p*lambda
  s<- c_aux
  r <- c_p
  Checks(lambda, mu_p, mu_aux, c_p, c_aux,p)
  #### Solve for Truncation parameter K ######
  X_i <- Solve_K(s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p, N = s+1)
  K <- ncol(X_i) - 1
  ### solve for A matrix ####
  #a_ij <- calc_a_ij()
  #A <- Calc_Am(K,s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p)

  ### solve for R
  #R <- Calc_R(A, K,s)

  ###solve for B
#  B <- Calc_Bmn(K,s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p)

  ###solve for X
#  X_i <- Calc_X(K,s,r, A,B,R)

  ### calculate results
  results <- Calc_results(p,lambda,lambda_aux,r,s,mu_p,mu_aux,K,X_i)
return(results)
  }
