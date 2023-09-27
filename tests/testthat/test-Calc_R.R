test_that("solving for R works", {
  ## test params
  K <- 2
  r <- 2
  s <- 1
  q <- .05
  lam <- 1
  mu_P <- 3
  mu_A <- 6
  p <- 1-q
  lam_P <- lam*p
  lam_A <- lam*q
  matrix_size <- (K) + 1
  A <- Calc_Am(K,s,r,lam,lam_P,lam_A,mu_P,mu_A, p)

  ## check determinate
  A_1 <- matrix(c(-(lam+2*mu_P),0,0,
                  mu_A*p,-(lam+mu_P +mu_A),  0,
                  0,mu_A,-(lam+mu_P+s*mu_A))
                , nrow = 3, ncol = 3, byrow = TRUE)
  I <- diag(1,nrow=matrix_size)
  det <- det(I - A[2,,])

  # Minor and cofactor
  minor <- function(A, i, j) det( A[-i,-j] )
  cofactor <- function(A, i, j) (-1)^(i+j) * minor(A,i,j)

  # With a loop
  adjoint1 <- function(A) {
    n <- nrow(A)
    B <- matrix(NA, n, n)
    for( i in 1:n )
      for( j in 1:n )
        B[j,i] <- cofactor(A, i, j)
    return(B)
  }
  adj <- adjoint1(I - A[2,,])
  temp_2 <- adj / det

  #check  det(I - A[2,,])
  temp <- (I - A[2,,])
  temp <- solve(temp)
  expect_equal(temp_2, temp)

  # Check matrix multiplication
  R_N1 <- A[1,,] %*% temp
  R_N2 <- crossprod(A[1,,], temp)
  expect_equal(R_N1, R_N2)

  # Check R calculation
  tol <- 0.01
  R_1 <- calc_R(A,K,s,tol)

  R_2 <- (A[1,,] + (temp %*% temp %*% A[3,,]) + (temp %*% temp %*% temp %*% A[4,,]) ) %*% temp
  expect_equal(round(R_1,6), round(R_2,6))

  })
