#' Title
#'
#' @param lambda
#' @param mu_p
#' @param mu_aux
#' @param c_p
#' @param c_aux
#' @param n_p
#' @param n_aux
#' @param p
#'
#' @return
#' @export
#'
#' @examples
Checks <- function(lambda, mu_p, mu_aux,c_p, c_aux,p)
{
  q <- 1-p
  lambda_aux <- q*lambda
  lambda_p <-p*lambda
  s<- c_aux
  r <- c_p

  AUX_r_c_warning <- "( lambda/(mu*c) ) has to be less than one!!"
  AUX_class <- "the class of the object x has to be M/M/C (i_AUX)"
  ALL_c_warning <- "the number of servers for each server type needs to be larger than 1"
  ALL_positive <- "Lambda and service times needs to be non-negative"
  ALL_p_warning <- "The probability of a auxiliary call needs to be between 0 and 1"
  AUX_anomalous <- "Some value of lambda, mu, c or n is anomalous. Check the values."
  ALL_integer <- "the  number of servers needs to be integer"

  s <- lambda_aux / mu_aux
  lam_1 <- c_p*mu_p
  lam_2 <- (c_p - c_aux)* mu_p + c_aux*mu_aux
  lam <- min(lam_1,lam_2)

  if ((c_p || c_s) < 1)
    stop(ALL_c_warning)

  if (lambda < 0)
    stop(ALL_positive)

  if ((mu_p || mu_s) <= 0)
    stop(ALL_positive)

  if (p < 0 || p > 1)
    stop(ALL_p_warning)

  if (lambda_aux >= c_aux*mu_aux ||  lambda >= lam )
  {
    ro <- max(lambda/(r*mu_p), lambda/((c_p - c_aux)* mu_p + c_aux*mu_aux))
    ro_a <- lambda_aux/(c_aux*mu_aux)
    #cat(paste("Throughput is: ", mu * c, "\n", sep=""))
    cat(paste("Utilization exceeds 100% use!! primary servers:", ro * 100, "%\n auxillary servers:", ro_a * 100, "%\n", sep=""))
    stop(AUX_r_c_warning)
  }

}