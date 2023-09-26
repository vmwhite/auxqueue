test_that("check A_k matrix", {
  ## test params
  K <- 4
  r <- 3
  s <- 1
  q <- .2
  lam <- 1
  mu_P <- 5
  mu_A <- 7
  p <- 1-q
  lam_P <- lam*p
  lam_A <- lam*q

  ##### hard coded matrices from Green 1994 page 1212-1213
  A_0 <-    matrix(c(lam,0,0,0,0,
                     0,lam,0,0, 0,
                     0, 0,lam,0,  0,
                     0, 0,0,lam,  0,
                     0, 0,0,0,lam)
                   , nrow = 5, ncol = 5, byrow = TRUE)
  A_1 <-    matrix(c(-(lam+3*mu_P),0,0,0,0,
                     mu_A*p,-(lam+calc_mu_mn(2,1,mu_P,mu_A)),0,0, 0,
                     0, mu_A*p,-(lam+calc_mu_mn(2,1,mu_P,mu_A)),0,  0,
                     0, 0,mu_A,-(lam+calc_mu_mn(2,1,mu_P,mu_A)),  0,
                     0,0,0,mu_A,-(lam+calc_mu_mn(2,1,mu_P,mu_A)))
                   , nrow = 5, ncol = 5, byrow = TRUE)
  A_2 <-    matrix(c(3*mu_P*p,0,0,0,0,
                     0, 2*mu_P*p + mu_A*q,0,0, 0,
                     0, 0,2*mu_P*p,0,  0,
                     0, 0,0,2*mu_P*p,  0,
                     0, 0,0,0,2*mu_P)
                   , nrow = 5, ncol = 5, byrow = TRUE)
  A_3 <-    matrix(c(0,3*mu_P*q,0,0,0,
                     0, 0,2*mu_P*calc_alpha_i(p,1),0, 0,
                     0, 0,0,2*mu_P*calc_alpha_i(p,1),  0,
                     0, 0,0,0, 2*mu_P*q, ## check if this should be q or calc_alpha_i(p,1)
                     0,0,0,0, 0)
                   , nrow = 5, ncol = 5, byrow = TRUE)

  A_4 <-    matrix(c(0,0,0,0,0,
                     0, 0,0, 2*mu_P*calc_alpha_i(p,2), 0,
                     0, 0,0,0,2*mu_P*(q^2),
                     0,0,0, 0, 0,
                     0,0,0,0, 0)
                   , nrow = 5, ncol = 5, byrow = TRUE)
  A_5 <-    matrix(c(0,0,0,0,0,
                     0, 0,0,0, 2*mu_P*(q^3),
                     0,0,0,  0, 0,
                     0,0,0,0, 0,
                     0,0,0,0,0)
                   , nrow = 5, ncol = 5, byrow = TRUE)

  #### calc Am
  A <- Calc_Am(K,s,r,lam,lam_P,lam_A,mu_P,mu_A, p)
  expect_equal(A[1,,], A_0)
  expect_equal(A[2,,], A_1)
  expect_equal(A[3,,], A_2)
  expect_equal(A[4,,], A_3)
  expect_equal(A[5,,], A_4)
  expect_equal(A[6,,], A_5)
})
