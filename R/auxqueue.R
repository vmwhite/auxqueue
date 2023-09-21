#'A QUEUEING SYSTEM WITH AUXILIARY SERVERS
#'
#' @param c_p := number of primary servers
#' @param c_aux := number of auxiliary servers
#' @param mu_p := service time of primary servers
#' @param mu_aux := service time of auxiliary servers
#' @param lambda := arrival rate of customers
#' @param p := percentage of customers that are of type primary
#' @param K := truncation parameter, default of 13
#' @param Stability:= default TRUE. TRUE is the queue checks for stability of inputs to 0.02 tolerance. False is with K being strictly equal to K.
#'
#' @return results in the form of a data frame
#' @export
#'
#' @examples
#' c_p <- 5
#' c_aux <- 3
#' mu_p <- 4
#' mu_aux <- 2
#' lambda <- 1
#' p <- .8
#' K <- 13
#' results<- auxqueue(c_p,c_aux,mu_p, mu_aux, lambda, p, K=13)
auxqueue <- function(c_p,c_aux,mu_p, mu_aux, lambda, p, K=13, Stability = TRUE){
  #### redefine parameters
  q <- 1-p
  lambda_aux <- q*lambda
  lambda_p <-p*lambda
  s<- c_aux
  r <- c_p
  skip<- FALSE
  if (grepl(  "ERROR",print(Checks(lambda, mu_p, mu_aux, c_p, c_aux,p)), fixed = TRUE) == TRUE){
    skip <- TRUE
    newlist$reduced_by <- 0
  }
  if (skip == FALSE){
    #### Solve for Truncation parameter K ######
    newlist<- Solve_K(s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p, K, Stability)
    X_i <- newlist$X
    reduced_by <- newlist$reduced_by
    #checking for reductions
    if(reduced_by > 0){
      K <- K - s
      s <- max(1,round(s/reduced_by))
      K <- K + s
      r <- max(1,round(r /reduced_by))
      lambda  <- lambda /reduced_by
      lambda_aux <- lambda *(1-p)
      lambda_p  <- lambda*p
    }else{
      K <- ncol(X_i) - 1
    }
    ### solve for A matrix ####
    #a_ij <- calc_a_ij()
    #A <- Calc_Am(K,s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p)

    ### solve for R
    #R <- Calc_R(A, K,s)

    ###solve for B
    #B <- Calc_Bmn(K,s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p)

    ###solve for X
    #X_i <- Calc_X(K,s,r, A,B,R)
    }
    ### calculate results
    results <- Calc_results(p,lambda,lambda_aux,r,s,mu_p,mu_aux,K,X_i, skip, newlist$reduced_by)

return(results)
  }
