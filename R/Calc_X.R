#' Solve for stationary probability vector X
#'
#' @param K := truncation parameter
#' @param s := number of auxiliary servers
#' @param r := number of primary servers
#' @param A :=  A_m matrix that represents non-boundary behavior
#' @param B := matrix B of the states with boundary behavior
#' @param R := R matrix
#'
#' @return stationary probability vector X
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
Calc_X <- function(K,s,r, A,B,R){
  #Calculate the truncated r by r Generator matrix
  G <- trunc_G(K,s,r, A,B,R)
  # Transpose matrix for solving for X
  G <- t(G)
  ## append final equation
  #x_0*e + x_1*e +.. x_r(I-R)^{-1}*e = 1
  # e is the column vector with all its components equal to 1
  e <- matrix(1, nrow = r, ncol=1)
  matrix_size = K+1
  # dummy x variable
  X_r <- matrix(1, nrow = 1, ncol=matrix_size)
  row_list <- c()
  I <- diag(1,matrix_size)
  a_Xr <- X_r %*% solve(I - R)
  for (i in 1:((r + 1) * matrix_size)){
    if (i <(((r * matrix_size)+1))){
      row_list<- append(row_list, 1)
      num <- 1
    }else{
      row_list<- append(row_list, a_Xr[num])
      num <- num+1
    }
  }
  # length(row_list) should == ncol(G)
  G <- rbind(G,row_list)

  b <- matrix(0, nrow=((r+1)*matrix_size),ncol=1)
  b <- rbind(b,1)


  # Use QR decomposition to solve the system
  t <- try(X <- qr.solve(qr(as.numeric(factor(G))), b))
  if("try-error" %in% class(t)){
    X <- lsfit(G, b)
    X <- X$coefficients
    #X <- ginv(G) %*% b ## takes longer
  }




  # add additional rows to X
  X <- normalize_vector(X,matrix_size,R)

  ## reformat to x_ij
  X_i <- matrix(0, nrow = max(length(X), nrow(X))/matrix_size, ncol=matrix_size) #depending on solution method use nrow or length
  count <- 1
  for (i in 1:(max(length(X), nrow(X))/matrix_size)){
    for(j in 1:matrix_size){
      X_i[i,j] <- X[count]
      count <- count+1
    }
  }
  return(X_i)
}
