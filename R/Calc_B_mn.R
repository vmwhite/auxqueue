#' Calculating B matrices
#'
#' @param K := truncation parameter
#' @param s := number of auxiliary servers
#' @param r := number of primary servers
#' @param lambda := arrival rate of all customers
#' @param lambda_aux := arrival rate of customers of type auxiliary
#' @param lambda_p := arrival rate of customers of type primary
#' @param mu_p := primary server rate
#' @param mu_aux := auxiliary server rate
#' @param p := percentage of customers that are of type primary
#'
#' @return matrix B of the states with boundary behavior
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
#' Calc_Bmn(K,s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p)
Calc_Bmn <- function(K,s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p){
  q <- 1 - p
  # + 1 since r indexes at 1
  matrix_size <- (K) + 1
  m_bmax <- max(K-s+2+1,r+1)
  n_bmax<- max( matrix_size,(r-(2*s)+ K+1)+1, r)
  B_mn <- matrix(0, nrow=matrix_size, ncol=matrix_size)
  B <- array(c(B_mn), c(m_bmax , n_bmax , matrix_size, matrix_size))
  for (m in 1:m_bmax){
    for (n in 1:n_bmax){
      i <- n
      #transition from state (n,) to  state(n - m + 1,)
      i_two <- n - m + 1 +1 #add one since R indexes at 1
      for (j in 1:matrix_size){
        #### if at least one P server is free #####
        if((i-1) + min((j-1),s) < r ){
          # for transitions from state  (i,) to (i+1,)
          if(i_two  == i + 1 ){
            if (j-1 < K){
              B[m,n,j,j] = lambda_p
            }else{
              B[m,n,j,j] = lambda
            }
            # for transitions from state  (i,) to (i,)
          }else if (i_two == i){
            B[m,n,j,j-1] = s*mu_aux
            if((j-1 +1) <= K){
              B[m,n,j,j+1] = lambda_aux
            }
            # for transitions from state  (i,) to (i-1,)
          } else if (i_two == i - 1){
            B[m,n,j,j] = (i-1)*mu_p# -1 since R indexes i at 1
          }
          #### if all p servers are busy, at least 1 aux server is free, and there is a primary queue #####
        }else if(((i-1) + min((j-1),s)) > r && (j-1) < s){
          # for transitions from state  (i,) to (i+1,)
          if(i_two  == i + 1 ){
            B[m,n,j,j] = lambda
            # for transitions from state  (i,) to (i,)
          }else if (i_two == i){
            B[m,n,j,j-1] = (j-1)*mu_aux*p
            # for transitions from state  (i,) to (i-1,)
          } else if (i_two == i - 1){
            B[m,n,j,j] = (r-(j-1))*mu_p*p + (j-1)*mu_aux*q  # -1 since R indexes j at 1
            # for transitions from state  (i,) to (i-2,)
          }else if (i_two == i - 2){
            B[m,n,j,j+1] =(r-(j-1))*mu_p*q # -1 since R indexes j at 1
          }
          #### if all p servers are busy, at least 1 aux server is free, and there is NO primary queue #####
        }else if((i-1) + min((j-1),s) == r && (j-1) < s){
          # for transitions from state  (i,) to (i+1,)
          if(i_two  == i + 1 ){
            B[m,n,j,j] = lambda
            # for transitions from state  (i,) to (i,)
          }else if (i_two == i){
            B[m,n,j,j-1] = (j-1)*mu_aux
            # for transitions from state  (i,) to (i-1,)
          } else if (i_two == i - 1){
            B[m,n,j,j] = (i-1)*mu_p  # -1 since R indexes i at 1
          }
          #### if all p and aux servers are busy #####
        }else if((i-1) + min((j-1),s) >= r && (j-1) >= s){
          # for transitions from state  (i,) to (i+1,)
          if(i_two  == i + 1 ){
            if( (i - 1) <= r - s - 1 && (j-1) < K ){
              B[m,n,j,j] = lambda_p
            }else{
              B[m,n,j,j] = lambda
            }
            # else transitions
          }else{
            #### AND NO primary queue #####
            if( (i-1) == r-s && (j-1)>=s ){
              # for transitions from state  (i,) to (i,)
              if (i_two == i ){
                  B[m,n,j,j-1] = s*mu_aux
              # for transitions from state  (i,) to (i-1,)
              } else if (i_two == i - 1){
                B[m,n,j,j] = (r-s)*mu_p
              }
              #### AND primary queue but NO aux queue#####
            }else if( (i-1)>= r-s && (j-1)==s ){
              # for transitions from state  (i,) to (i,)
              if (i_two == i){
                B[m,n,j,j-1] = s*mu_aux
                # for transitions from state  (i,) to (i-1,)
              } else if (i_two == i - 1){
                B[m,n,j,j] = (r-s)*mu_p*p + s*mu_aux*q
                # transition from (i,j) to  (i - k -1, j+k), k ==0 is to (i-1) state
              }else{
                for (k in 1:((i-1)-r+s)){
                  if( i_two == i - k - 1 && (((j-1)+ k == K) || k == (i-1)-r+s )){
                    # if max aux queue is reached or if primary queue is exhausted
                    if ((j-1) +k <= K){
                      B[m,n,j,j+k] =  (r-s)*mu_p*q^k
                    }
                  }else if ( i_two == i - k - 1 && (((j-1)+k < K))){
                    B[m,n,j,j+k] = (r-s)*mu_p*(q^k)*p
                  }
                }
              }
              #### AND primary queue AND aux queue#####
            }else if( (i-1)> r-s && (j-1)>s ){
              # for transitions from state  (i,) to (i,)
              if (i_two == i){
                B[m,n,j,j-1] =s*mu_aux
                # transition from (i,j) to  (i - k -1, j+k), k > 0 since there is a primary queue
              }else{
                for (k in 0:((i-1)-r+s)){
                  if( i_two == i - k - 1 && (((j-1)+ k == K) || k == (i-1)-r+s )){
                    # if max aux queue is reached or if primary queue is exhausted
                    if ((j-1) +k <= K){
                      B[m,n,j,j+k] =  (r-s)*mu_p*q^k
                     }
                    }else if ( i_two == i - k - 1 && (((j-1)+k < K))){
                    B[m,n,j,j+k] = (r-s)*mu_p*(q^k)*p
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  m_i_i = 2
  for (row in 1:matrix_size){
    for (col in 1:matrix_size){
      if (row == col){
        for (n in 1:n_bmax){
          for (m in 1:m_bmax){
            for (col_2 in 1:matrix_size){
              if ((m_i_i!= m || col!=col_2)){
                B[m_i_i,n,row,col] =  B[m_i_i,n,row,col]- B[m,n,row,col_2]
              }
            }
          }
        }
      }
    }
  }
  return(B)
}
