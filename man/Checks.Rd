% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Checks.R
\name{Checks}
\alias{Checks}
\title{Checking for queue necessary conditions}
\usage{
Checks(lambda, mu_p, mu_aux, c_p, c_aux, p)
}
\arguments{
\item{lambda}{:= arrival rate of all customers}

\item{mu_p}{:= primary server rate}

\item{mu_aux}{:= auxiliary server rate}

\item{c_p}{:= number of primary servers}

\item{c_aux}{:= number of auxiliary servers}

\item{p}{:= percentage of customers that are of type primary}
}
\value{
nothing or error string if necessary conditions are not met
}
\description{
Checking for queue necessary conditions
}
\examples{
r <- 10
s <- 5
c_p <- r
c_aux <- s
lambda <- .2
p <-.9
mu_p <-.4
mu_aux <-.5
Checks(lambda, mu_p, mu_aux,c_p, c_aux,p)
}
