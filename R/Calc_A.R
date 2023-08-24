### Calculating single A matrix
#' Title
#'
#' @param A
#' @param matrix_size
#' @param A_m_lim
#'
#' @return
#' @export
#'
#' @examples
Calc_A <- function(A, matrix_size, A_m_lim){
  A_matrix<- matrix(0, nrow=matrix_size, ncol=matrix_size)
  for (m in 1:A_m_lim){
    A_matrix <- A_matrix + A[m,,]
  }
}
