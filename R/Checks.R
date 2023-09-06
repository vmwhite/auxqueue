#' Checking for queue necessary conditions
#' @param lambda := arrival rate of all customers
#' @param mu_p := primary server rate
#' @param mu_aux := auxiliary server rate
#' @param c_p := number of primary servers
#' @param c_aux := number of auxiliary servers
#' @param p := percentage of customers that are of type primary
#'
#' @return nothing or error string if necessary conditions are not met
#' @export
#'
#' @examples
#' r <- 10
#' s <- 5
#' c_p <- r
#' c_aux <- s
#' lambda <- .2
#' p <-.9
#' mu_p <-.4
#' mu_aux <-.5
#' Checks(lambda, mu_p, mu_aux,c_p, c_aux,p)
Checks <- function(lambda, mu_p, mu_aux,c_p, c_aux,p)
{
  q <- 1-p
  lambda_aux <- q*lambda
  lambda_p <-p*lambda
  s<- c_aux
  r <- c_p

  AUX_r_c_warning <- "ERROR: ( lambda/(mu*c) ) has to be less than one!!"
  AUX_class <- "ERROR: the class of the object x has to be M/M/C (i_AUX)"
  ALL_c_warning <- "ERROR: the number of servers for each server type needs to be larger than 1"
  ALL_positive <- "ERROR: Lambda and service times needs to be non-negative"
  ALL_p_warning <- "ERROR: The probability of a auxiliary call needs to be between 0 and 1"
  AUX_anomalous <- "ERROR: Some value of lambda, mu, c or n is anomalous. Check the values."
  ALL_integer <- "ERROR: the  number of servers needs to be integer"

  s <- lambda_aux / mu_aux
  lam_1 <- c_p*mu_p
  lam_2 <- (c_p - c_aux)* mu_p + c_aux*mu_aux
  lam <- min(lam_1,lam_2)

  if ((c_p || c_s) < 1)
    paste0(ALL_c_warning)

  if (lambda < 0)
    paste0(ALL_positive)

  if ((mu_p || mu_s) <= 0)
    paste0(ALL_positive)

  if (p < 0 || p > 1)
    paste0(ALL_p_warning)

  if (lambda_aux >= c_aux*mu_aux ||  lambda >= lam )
  {
    ro <- max(lambda/(r*mu_p), lambda/((c_p - c_aux)* mu_p + c_aux*mu_aux))
    ro_a <- lambda_aux/(c_aux*mu_aux)
    #cat(paste("Throughput is: ", mu * c, "\n", sep=""))
    paste0(AUX_r_c_warning, " Utilization exceeds 100% use!! primary servers rho:", ro * 100, "%\n auxillary servers rho:", ro_a * 100, "%\n", sep="")
  }

}
