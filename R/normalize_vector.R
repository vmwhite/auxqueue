#' Calculate additional stationary transition probabilities above boundary behavior
#'
#' @param vec <- vector of X_i to add transition probabilities to
#' @param tolerance <- tolerance to stop loop at
#' @param time_limit <- time limit in seconds
#'
#' @return returns expanded transtion probabilty vector
#' @export
#'
#' @examples
#' #' K <- 13
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
#' vec <- Calc_X(K,s,r, A,B,R)
#' normalize_vector(vec, matrix_size, R)
normalize_vector <- function(vec, matrix_size, R, tolerance = 0.001, time_limit= 5) {
  # Get the start time
  start_time <- Sys.time()
  while ( 1 - sum(vec, na.rm=TRUE) > tolerance && difftime(Sys.time(), start_time, units = "secs") < time_limit ) {
    new_X <- vec[-(1:(length(vec) - matrix_size))] %*% R
    for (i in 1:ncol(new_X)){
      vec <- rbind(vec, new_X[i] )
    }
  }
  return(vec)
}
