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
  matrix_size <- (K) + 1
  aij <- matrix(0, nrow=1,ncol=matrix_size)
  a_irow_icol_jrow_jcol <- array(c(aij), c(r+1, r+1, matrix_size, matrix_size))
  inner_col2 <- 0
  A_count <- 0

  for(i in (0:r-1)){
    n <- i
    for (j in (0:(r-1))){
      m <- n - j + 1
      if (m > 0){
        a_irow_icol_jrow_jcol[i+1,j+1,,] <- B[m+1,n+1,,]
      }
    }
  }
  ### summation for final row
  row <- r+1
  inner_col <- 1
  A_count <- 0
  for (col in 1:(r+1)){
    r_row_sum <- 0
    ## first summation
    if (col == (r-s-1)+1){
      if (K-2*s > 0){
        for (i in 0:(K-2*s)){
          r_row_sum <- r_row_sum + (matrix_power(R,i) %*% B[(s+i+2)+1,(r+i)+1,,])
        }
      }else if (K-2*s == 0){
        r_row_sum <- r_row_sum + (matrix_power(R,0) %*% B[(s+0+2)+1,(r+0)+1,,])
      }
      ## last summation
    }else if(col == (r)+1 && (K-s+2-1) >= 0){
      for (i in 0:(K-s+2-1)){
        r_row_sum <- r_row_sum + (matrix_power(R,i) %*% A[(i+1)+1,,])
      }
      ### (r-s-1) + 1 to (r) + 1 summation ####
    }else if (col < ((r) +1) && (col>((r-s-1)+1) )){
      ## sum of the B's
      B_count <- max(K-s+1 + 1 - A_count, 2+1)
      if (B_count == 3){
        if (K-2*s == 0){
          r_row_sum <- r_row_sum + (matrix_power(R,0) %*% B[2+1,(r+0)+1,,])
        }else if (K-2*s > 0){
          for (i in 0:(K-2*s)){
            r_row_sum <- r_row_sum + (matrix_power(R,i) %*% B[2+1,(r+i)+1,,])
          }
        }
      }else{
        if (K-2*s == 0){
          r_row_sum <- r_row_sum + (matrix_power(R,0) %*% B[(s+0+1)+1,(r+0)+1,,])
        }else if (K-2*s > 0){
          for (i in 0:(K-2*s)){
            r_row_sum <- r_row_sum + (matrix_power(R,i) %*% B[(s+i+1)+1,(r+i)+1,,])
          }
        }
      }
      ## sum the A's
      max <- (K-s-r)
      if (((K-s-r) <= (K-2*s +1)) && A_count == 0 && (K-s+2) >= 0 ){
        r_row_sum <- r_row_sum + (matrix_power(R,(K-2*s+1)) %*% A[(K-s+2)+1,,])
      }else if ((A_count) == 0 && (K-s+2) >= 0){
        r_row_sum <- r_row_sum + (matrix_power(R,(K-2*s+1)) %*% A[(K-s+2)+1,,])
      }else if ((K-s+2) >= 0){
        for (i in (K-2*s +1):(K-2*s +1 + A_count)){
          if (i <= K-s-r){
            r_row_sum <- r_row_sum + (matrix_power(R,i) %*% A[(r+i +2)+1,,])
          }
        }
      }
      A_count <- A_count + 1
    }
    a_irow_icol_jrow_jcol[row,col,,] <- r_row_sum
  }


  ### restructure the matrix to solve
  # Reshape the matrix to a 2-dimensional matrix of size (r+1)*matrix_size by (r+1)*matrix_size
  # Assuming r and matrix_size are already defined
  a <- matrix(0, nrow = (r + 1) * matrix_size, ncol = (r + 1) * matrix_size)
  # Calculate the number of rows and columns in a_irow_icol_jrow_jcol
  a_rows <- nrow(a_irow_icol_jrow_jcol)
  a_cols <- ncol(a_irow_icol_jrow_jcol)

  # Calculate the number of blocks in the reshaped matrix
  blocks_rows <- matrix_size
  blocks_cols <-  matrix_size

  for (i in 1:(r+1)) {
    for (j in 1:(r+1)) {
      # Calculate the indices for block assignment
      start_row <- (i - 1) * blocks_rows + 1
      end_row <- i * blocks_rows
      start_col <- (j - 1) * blocks_cols + 1
      end_col <- j * blocks_cols

      # Assign the block from a_irow_icol_jrow_jcol to a
      a[start_row:end_row, start_col:end_col] <- a_irow_icol_jrow_jcol[i, j, , ]
    }
  }
  a <- t(a)
  # e is the column vector with all its components equal to 1
  e <- matrix(1, nrow = r, ncol=1)


  X_r <- matrix(1, nrow = 1, ncol=matrix_size)


  ## append final equation
  #x_0*e + x_1*e +.. x_r(I-R)^{-1}*e = 1
  row_list <- c()
  I <- diag(1,matrix_size)
  a_Xr <- X_r %*% solve(I - R)
  a_Xr[2]
  for (i in 1:((r + 1) * matrix_size)){
    if (i <(((r * matrix_size)+1))){
      row_list<- append(row_list, 1)
      num <- 1
    }else{
      row_list<- append(row_list, a_Xr[num])
      num <- num+1
    }
  }
  length(row_list)
  ncol(a)
  a <- rbind(a,row_list)

  b <- matrix(0, nrow=((r+1)*matrix_size),ncol=1)
  b <- rbind(b,1)


  # Use QR decomposition to solve the system
  qr_decomp <- qr(a)
  X <- qr.solve(qr_decomp, b)
  best_fit <- lsfit(a, b)

  # Set the time limit in seconds
  time_limit <- 5
  # Function to normalize a vector
  normalize_vector <- function(vec, tolerance = 0.00, time_limit= 5) {
    # Get the start time
    start_time <- Sys.time()
    while ( 1 - sum(vec) > tolerance && difftime(Sys.time(), start_time, units = "secs") < time_limit ) {
      new_X <- vec[-(1:(length(vec) - matrix_size))] %*% R
      for (i in 1:ncol(new_X)){
        vec <- rbind(vec, new_X[i] )
      }
    }
    return(vec)
  }
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
