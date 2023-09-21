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
#' @param skips := T/F that is True if the queue is unstable
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
#' skip <- FALSE
#' A <- Calc_Am(K,s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p)
#' R <- Calc_R(A, K,s)
#' B <-  Calc_Bmn(K,s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p)
#' X_i <- Calc_X(K,s,r, A,B,R)
#' Calc_results(p,lambda,lambda_aux,r,s,mu_p,mu_aux,K,X_i, skip)

Calc_results <-function(p,lambda,lambda_aux,r,s,mu_p,mu_aux,K,X_i, skip, reduced_by){
  if (skip == FALSE){
    ### L^p_N = mean number of individuals in the primary queue when the truncation is N #skip
    ### L^a_N = mean number of individuals in the auxiliary queue when the truncation is N #skip
    ### Traffic intensity of type A customer
    rho_A <- lambda_aux /(s*mu_aux)
    ### Traffic intensity of all customers
    rho <- lambda / (r*mu_p)
    ### the probability of a Type A customer waiting while a primary server is free
    alpha <- 0
    for(i in 1:((r-s))){
      if (s+2 == ncol(X_i)){
        alpha <- alpha + X_i[i,ncol(X_i)]
      }else{
        for(j in (s+2):ncol(X_i)){
        alpha <- alpha + X_i[i,j]
        }
      }
    }
    ### steady state number customers waiting in the primary queue
    L_P_q <- 0 ### steady state number customers waiting in the primary queue
    p_d_P<- 0  ### the probability of delay for Type P customer
    p_d_A <- 0 ### the probability of delay for Type a customer
    for (j in 1:(s-1+1)){
      for (i in (r-j+1):nrow(X_i)){
        num_in_queue <- (i - r + j - 1)
        L_P_q <- L_P_q + num_in_queue*X_i[i,j]
        p_d_P <- p_d_P+ (1 - ((1-p)^num_in_queue))*X_i[i,j]
        p_d_A <- p_d_A+ (1 - ((p)^num_in_queue))*X_i[i,j]
      }
    }
    for (j in (s+1):(K+1)){
      for (i in (r-s+1+1):nrow(X_i)){
        num_in_queue <- (i - r + s - 1)
        L_P_q <- L_P_q + num_in_queue*X_i[i,j]
        p_d_P <- p_d_P+ (1 - ((1-p)^num_in_queue))*X_i[i,j]
        p_d_A <- p_d_A+ (1 - ((p)^num_in_queue))*X_i[i,j]
      }
    }

    ### steady state number customers waiting in the auxiliary queue
    L_A_q <- 0
    for (i in 1:nrow(X_i)){
      for (j in (s+1+1):(K+1)){
        L_A_q <- L_A_q + (j - s - 1)*X_i[i,j]
        p_d_A <- p_d_A + X_i[i,j]
      }
    }

    ### expected waiting time in queue for primary customers, the mean delay for type p customers
    W_P_q = L_P_q / lambda
    ### expected waiting time in queue for auxiliary customers, the mean delay for type A customers
    W_A_q = (L_P_q / lambda) + (L_A_q/lambda_aux)


    ### the probability regular service occurs if aux queue is too long
    #the probability a type A arrival is the next event, where primary service doesnt matter
    prob <- lambda_aux / (lambda + mu_aux + mu_p)
    beta <- X_i[1,K+1] * (lambda_aux / (lambda + mu_aux ))
    # if there is at least 1 primary call being served
    for(i in 2:(ncol(X_i))){
      beta <- beta + X_i[i,K+1]
    }
  }
  #times the probability of an aux call being the next event

  results <- list(r,s, K,mu_p, mu_aux, lambda, p, reduced_by)

  if (skip == FALSE){
    results <- append(results,rho)
    results <- append(results,rho_A)
    results <- append(results,alpha)
    results <- append(results,L_P_q)
    results <- append(results,L_A_q)
    results <- append(results,p_d_P)
    results <- append(results,p_d_A)
    results <- append(results,W_P_q)
    results <- append(results,W_A_q)
    results <- append(results,beta)
  }else{
    for (i in 1:10){
      results <- append(results, "unstable")
    }
  }
  #reformat results
  #create DF of results
  DF <- data.frame(results)
  # set column names
  metrics <- list("r", "s","K", "mu_p", "mu_aux","lambda", "p", "reduced_by" )
  metrics <- append(metrics, "rho_P" )
  metrics <- append(metrics, "rho_A" )
  metrics <- append(metrics, "alpha_typeA_prob_delay" )
  metrics <- append(metrics, "L_p_q" )
  metrics <- append(metrics, "L_A_q" )
  metrics <- append(metrics, "p_d_P")
  metrics <- append(metrics, "p_d_A")
  metrics <- append(metrics, "W_P_q" )
  metrics <- append(metrics, "W_A_q" )
  metrics <- append(metrics, "beta_typeA_coveredby_Pserver" )
  colnames(DF) <- metrics
  return(DF)
}
