# Solve for stationary probability vector X
#' Title
#'
#' @param K
#' @param s
#' @param r
#' @param A
#' @param B
#' @param R
#'
#' @return
#' @export
#'
#' @examples
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
  # length(row_list) should == ncol(a)
  G <- rbind(G,row_list)

  b <- matrix(0, nrow=((r+1)*matrix_size),ncol=1)
  b <- rbind(b,1)


  # Use QR decomposition to solve the system
  qr_decomp <- qr(G)
  X <- qr.solve(qr_decomp, b)
  best_fit <- lsfit(G, b)

  # normalize X
  X <- normalize_vector(X)

  ## reformat to x_ij
  X_i <- matrix(0, nrow = nrow(X)/matrix_size, ncol=matrix_size)
  count <- 1
  for (i in 1:(nrow(X)/matrix_size)){
    for(j in 1:matrix_size){
      X_i[i,j] <- X[count]
      count <- count+1
    }
  }
  return(X_i)
}
