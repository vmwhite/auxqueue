# Function to normalize a vector
#' Title
#'
#' @param vec
#' @param tolerance
#' @param time_limit
#'
#' @return
#' @export
#'
#' @examples
normalize_vector <- function(vec, matrix_size, R, tolerance = 0.001, time_limit= 5) {
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
