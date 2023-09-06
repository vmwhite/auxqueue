#' Calculate metrics of A QUEUEING SYSTEM WITH AUXILIARY SERVERS
#'
#' @param lambda := arrival rate of all customers
#' @param lambda_aux := arrival rate of customers of type auxiliary
#' @param r := number of primary servers
#' @param s := number of auxiliary servers
#' @param mu_p := primary server rate
#' @param mu_aux := auxiliary server rate
#' @param K := truncation parameter
#' @param X_i := Array of transition probabilities
#'
#' @return returns queue key performance metrics in a data frame
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
#' R <- Calc_R(A, K,s)
#' B <-  Calc_Bmn(K,s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p)
#' X_i <- Calc_X(K,s,r, A,B,R)
#' Calc_results(p,lambda,lambda_aux,r,s,mu_p,mu_aux,K,X_i)

Calc_results <-function(p,lambda,lambda_aux,r,s,mu_p,mu_aux,K,X_i){
  ### L^p_N = mean number of individuals in the primary queue when the truncation is N #skip
  ### L^a_N = mean number of individuals in the auxiliary queue when the truncation is N #skip
  ### Traffic intensity of type A customer
  rho_A <- lambda_aux /(s*mu_aux)
  ### Traffic intensity of all customers
  rho <- lambda / (r*mu_p)
  ### the probability of a Type A customer waiting while a primary server is free
  alpha <- 0
  for(i in 1:((r-s))){
    for(j in (s+2):ncol(X_i)){
      alpha <- alpha + X_i[i,j]
    }
  }

  ### steady state number customers waiting in the primary queue
  L_P_q <- 0
  for (j in 1:(s-1+1)){
    for (i in (r-j+1):nrow(X_i)){
      L_P_q <- L_P_q + (i - r + j - 2)*X_i[i,j]
    }
  }
  for (j in (s+1):(K+1)){
    for (i in (r-s+1+1):nrow(X_i)){
      L_P_q <- L_P_q + (i - r + s - 1)*X_i[i,j]
    }
  }

  ### steady state number customers waiting in the auxilary queue
  L_A_q <- 0

  for (i in 1:nrow(X_i)){
    for (j in (s+1+1):(K+1)){
      L_A_q <- L_A_q + (j - s - 1)*X_i[i,j]
    }
  }

  ### expected waiting time in queue for primary customers, the mean delay for type p customers
  W_P_q = L_P_q / lambda
  ### expected waiting time in queue for auxiliary customers, the mean delay for type A customers
  W_A_q = (L_P_q / lambda) + (L_A_q/lambda_aux)
  results <- c(r)
  results <- append(results,s)
  results <- append(results,K)
  results <- append(results,mu_p)
  results <- append(results,mu_aux)
  results <- append(results,lambda)
  results <- append(results,p)
  results <- append(results,rho_A)
  results <- append(results,alpha)
  results <- append(results,L_P_q)
  results <- append(results,L_A_q)
  results <- append(results,W_P_q)
  results <- append(results,W_A_q)
  return(results)
}
