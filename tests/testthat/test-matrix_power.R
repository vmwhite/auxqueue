test_that("multiplication works", {
  K <- 13
  r <- 10
  s <- 5
  lambda <- .2
  p <-.9
  lambda_aux <- (p-1) * lambda
  lambda_p <- p * lambda
  mu_p <-.4
  mu_aux <-.5
  A <- Calc_Am(K,s,r,lambda,lambda_aux,lambda_p,mu_p,mu_aux, p)

  test <- matrix_power(A, 5)
  expect_equal(test, A %*% A %*% A %*% A %*% A)
  test <- matrix_power(A, 3)
  expect_equal(test, A %*% A %*% A)
  test <- matrix_power(A, 1)
  expect_equal(test, A )
  test <- matrix_power(A, 0)
  expect_equal(test, diag(nrow(A)) )
})
