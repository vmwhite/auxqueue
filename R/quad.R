#' quadratic function
#'
#' @param a coefficient of X^2
#' @param b coefficient of X
#' @param c constant
#'
#' @return Two X values from resulting quadratic equation
#' @export
#'
#' @examples
#' a <- 2
#' b <- -1
#' c <- -3
#' quad(a,b,c)
quad <- function(a, b, c)
{
  a <- as.complex(a)
  answer <- c((-b + sqrt(b^2 - 4 * a * c)) / (2 * a),
              (-b - sqrt(b^2 - 4 * a * c)) / (2 * a))
  if(all(Im(answer) == 0)) answer <- Re(answer)
  if(answer[1] == answer[2]) return(answer[1])
  answer
}
