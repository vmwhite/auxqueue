#' Title
#'
#' @param A
#' @param n
#'
#' @return
#' @export
#'
#' @examples
matrix_power <- function(A, n) {
  if (n == 0){
    A_power <- diag(nrow(A))
  }else if (n==1){
    A_power <- A
  }else{
    A_power <- A
    for (i in 2:n) {
      A_power<- A_power %*% A
    }
  }
  return(A_power)
}
