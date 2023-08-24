#' Calculating A matrix
#'
#' @param K
#' @param s
#' @param r
#' @param lambda
#' @param lambda_aux
#' @param lambda_p
#' @param mu_p
#' @param mu_aux
#' @param p
#'
#' @return
#' @export
#'
#' @examples
Calc_Am <- function(K,s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p){
  #define q
  q <- 1-p

  # Auxillary queue is never larger than K-s, therefore j <= K
  matrix_size <- K+1 #add one since R indexes at 1

  # A go from 0 to K-s+2
  A_m_lim <- K - s + 2 + 1 #add one since R indexes at 1

  #Create Empty A blocks
  #A_m is the transition from (i,) to (i + -m + 1,)
  A_m <- matrix(0, nrow=matrix_size, ncol=matrix_size)
  A <- array(c(A_m), c(A_m_lim, matrix_size, matrix_size))

  #A blocks are the same for all i >= r-2s+ K + 1
  i <- r - 2*s + K +1 +1  #add one since R indexes at 1

  #fill A blocks
  for (m in 1:A_m_lim ){
    #transition to  state(i + -m + 1,)
    i_two <- i - m + 1 + 1 #add one since R indexes at 1
    for (j in 1:matrix_size){
      #### if at least one P server is free #####
      if((i-1) + (j-1) < r ){
        # for transitions from state  (i,) to (i+1,)
        if(i_two  == i + 1 ){
          A[m,j,j] = lambda_p
          # for transitions from state  (i,) to (i,)
        }else if (i_two == i){
          A[m,j,j] = lambda
          A[m,j,j-1] = s*mu_aux
          # for transitions from state  (i,) to (i-1,)
        } else if (i_two == i - 1){
          A[m,j,j] = (i-1)*mu_p# -1 since R indexes i at 1
        }
        #### if all p servers are busy, at least 1 aux server is free, and there is a primary queue #####
      }else if((i-1) + (j-1) > r && (j-1) < s){
        # for transitions from state  (i,) to (i+1,)
        if(i_two  == i + 1 ){
          A[m,j,j] = lambda
          # for transitions from state  (i,) to (i,)
        }else if (i_two == i){
          A[m,j,j-1] = (j-1)*mu_aux*p
          # for transitions from state  (i,) to (i-1,)
        } else if (i_two == i - 1){
          A[m,j,j] = (r-(j-1))*mu_p*p + (j-1)*mu_aux*q  # -1 since R indexes j at 1
          # for transitions from state  (i,) to (i-2,)
        }else if (i_two == i - 2){
          A[m,j,j+1] =(r-(j-1))*mu_p*q # -1 since R indexes j at 1
        }
        #### if all p servers are busy, at least 1 aux server is free, and there is NO primary queue #####
      }else if((i-1) + (j-1) == r && (j-1) < s){
        # for transitions from state  (i,) to (i+1,)
        if(i_two  == i + 1 ){
          A[m,j,j] = lambda
          # for transitions from state  (i,) to (i,)
        }else if (i_two == i){
          A[m,j,j-1] = (j-1)*mu_aux
          # for transitions from state  (i,) to (i-1,)
        } else if (i_two == i - 1){
          A[m,j,j] = (i-1)*mu_p  # -1 since R indexes i at 1
        }
        #### if all p and q servers are busy #####
      }else if((i-1) + (j-1) >= r && (j-1) >= s){
        # for transitions from state  (i,) to (i+1,)
        if(i_two  == i + 1 ){
          A[m,j,j] = lambda
          # else transitions
        }else{
          #### AND NO primary queue #####
          if( (i-1) == r-s && (j-1)>=s){
            # for transitions from state  (i,) to (i,)
            if (i_two == i){
              A[m,j,j-1] = s*mu_p
              # for transitions from state  (i,) to (i-1,)
            } else if (i_two == i - 1){
              A[m,j,j] = (r-s)*mu_p
            }
            #### AND primary queue but NO aux queue#####
          }else if( (i-1)> r-s && (j-1)==s ){
            # for transitions from state  (i,) to (i,)
            if (i_two == i){
              A[m,j,j-1] = s*mu_aux*p
              # for transitions from state  (i,) to (i-1,)
            } else if (i_two == i - 1){
              A[m,j,j] = (r-s)*mu_p*p + s*mu_aux*q
              # transition from (i,j) to  (i - k -1, j+k), k ==0 is to (i-1) state
            }else{
              for (k in 1:((i-1)-r+s)){
                if( k < (i-1) - r + s &&  i_two == i - k - 1 && (((j-1)+k < K))){
                  A[m,j,j+k] = (r-s)*mu_p*(q^k)*p
                }else if (k <= (i-1) - r + s && i_two == i - k - 1 && (j-1)+ k <= K){
                  A[m,j,j+k] = (r-s)*mu_p*q^k
                }
              }
            }
            #### AND primary queue AND aux queue#####
          }else if( (i-1)> r-s && (j-1)>s ){
            # for transitions from state  (i,) to (i,)
            if (i_two == i){
              A[m,j,j-1] =s*mu_aux
              # transition from (i,j) to  (i - k -1, j+k), k > 0 since there is a primary queue
            }else{
              for (k in 0:((i-1)-r+s)){
                if( k < (i-1) - r + s &&  i_two == i - k - 1 && (((j-1)+k < K) )){
                  A[m,j,j+k] = (r-s)*mu_p*(q^k)*p
                }else if (k <= (i-1) - r + s && i_two == i - k - 1 && (j-1)+ k <= K){
                  A[m,j,j+k] =  (r-s)*mu_p*q^k
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
        for (m in 1:A_m_lim){
          for (col_2 in 1:matrix_size){
            if ((m_i_i!= m || col!=col_2) && A_m_lim >= m_i_i){
              A[m_i_i,row,col] = A[m_i_i,row,col]-A[m,row,col_2]
            }
          }
        }
      }
    }
  }
  return(A)
}
