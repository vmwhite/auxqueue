test_that("Test B matrix calculations", {
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
  B_00 = matrix(c(lam_P,0,0, 0, 0,
                  0,lam_P,0,  0, 0,
                  0,0,lam_P,  0, 0,
                  0,0,0,lam_P, 0,
                  0,0,0,0,lam)
                , nrow = 5, ncol = 5, byrow = TRUE)
  B_02 = matrix(c(lam_P,0,0, 0, 0,
                  0,lam,0,  0, 0,
                  0,0,lam,  0, 0,
                  0,0,0,lam, 0,
                  0,0,0,0,lam)
                , nrow = 5, ncol = 5, byrow = TRUE)
  B_10 = matrix(c(-(lam+(0*mu_P)),lam_A,0, 0, 0,
                  mu_A,-(lam +calc_mu_mn(0,1, mu_P,mu_A)),lam_A,  0, 0,
                  0,mu_A,-(lam +calc_mu_mn(0,1, mu_P,mu_A)),  lam_A, 0,
                  0,0, mu_A,-(lam +calc_mu_mn(0,1, mu_P,mu_A)),  lam_A,
                  0,0,0,mu_A,-(lam +calc_mu_mn(0,1, mu_P,mu_A)))
                , nrow = 5, ncol = 5, byrow = TRUE)
  B_11 = matrix(c(-(lam+(1*mu_P)),lam_A,0, 0, 0,
                  mu_A,-(lam +calc_mu_mn(1,1, mu_P,mu_A)),lam_A,  0, 0,
                  0,mu_A,-(lam +calc_mu_mn(1,1, mu_P,mu_A)),  lam_A, 0,
                  0,0, mu_A,-(lam +calc_mu_mn(1,1, mu_P,mu_A)),  lam_A,
                  0,0,0,mu_A,-(lam +calc_mu_mn(1,1, mu_P,mu_A)))
                , nrow = 5, ncol = 5, byrow = TRUE)
  ## Error in Green Text, transition from (2,2+x) to (2,3+x) where 0<=x<=K-3 should be equal to 0 not lam_A
   B_12 = matrix(c(-(lam+(2*mu_P)),lam_A,0, 0, 0,
                  mu_A,-(lam +calc_mu_mn(2,1, mu_P,mu_A)),0,  0, 0,
                  0,mu_A,-(lam +calc_mu_mn(2,1, mu_P,mu_A)),  0, 0,
                  0,0, mu_A,-(lam +calc_mu_mn(2,1, mu_P,mu_A)),  0,
                  0,0,0,mu_A,-(lam +calc_mu_mn(2,1, mu_P,mu_A)))
                , nrow = 5, ncol = 5, byrow = TRUE)
  B_21 = matrix(c(mu_P,0,0, 0, 0,
                  0,mu_P,0,  0, 0,
                  0,0,mu_P,  0, 0,
                  0,0,0,mu_P, 0,
                  0,0,0,0,mu_P)
                , nrow = 5, ncol = 5, byrow = TRUE)
  B_22 = matrix(c(2*mu_P,0,0, 0, 0,
                  0,2*mu_P,0,  0, 0,
                  0,0,2*mu_P,  0, 0,
                  0,0,0,2*mu_P, 0,
                  0,0,0,0, 2*mu_P)
                , nrow = 5, ncol = 5, byrow = TRUE)
  B_33 = matrix(c(0,0,0,0, 0,
                  0,0, 2*mu_P*q, 0, 0,
                  0,0,0, 2*mu_P*q, 0,
                  0,0,0, 0, 2*mu_P*q,
                  0,0,0,0, 0)
                , nrow = 5, ncol = 5, byrow = TRUE)
  B_34 = matrix(c(0,3*mu_P*q,0,0, 0,
                  0,0,2*mu_P*calc_alpha_i(p,1),0, 0,
                  0,0,0, 2*mu_P*calc_alpha_i(p,1), 0,
                  0,0,0,  0, 2*mu_P*q,
                  0,0,0,0, 0)
                , nrow = 5, ncol = 5, byrow = TRUE)
  B_44 = matrix(c(0,0,0,0, 0,
                  0,0,0, 2*mu_P*(q^2), 0,
                  0,0,0, 0, 2*mu_P*(q^2),
                  0,0,0,  0,0, 0,
                  0,0,0,0, 0, 0)
                , nrow = 5, ncol = 5, byrow = TRUE)
  B_45 = matrix(c(0,0,0,0, 0,
                  0,0,0,  2*mu_P*calc_alpha_i(p,2), 0,
                  0,0,0,   0, 2*mu_P*(q^2),
                  0,0,0,  0, 0,
                  0,0,0,0, 0)
                , nrow = 5, ncol = 5, byrow = TRUE)
  B_55 = matrix(c(0,0,0,0, 0,
                  0,0,0, 0, 2*mu_P*(q^3),
                  0,0,0,  0, 0,
                  0,0,0,  0,0,
                  0,0,0,0, 0)
                , nrow = 5, ncol = 5, byrow = TRUE)
  B_56 = matrix(c(0,0,0,0, 0,
                  0,0,0, 0, 2*mu_P*(q^3),
                  0,0,0,  0, 0,
                  0,0,0,  0,0,
                  0,0,0,0, 0)
                , nrow = 5, ncol = 5, byrow = TRUE)

  B <- Calc_Bmn(K,s,r,lam,lam_A,lam_P,mu_P,mu_A, p)

  #Check values -- Note: R indexes at 1
  expect_equal(B[1,1,,], B_00)
  expect_equal(B[1,2,,], B_00)
  expect_equal(B[1,3,,], B_02)
  expect_equal(B[2,1,,], B_10)
  expect_equal(B[2,2,,], B_11)
  expect_equal(B[2,3,,], B_12)
  expect_equal(B[3,2,,], B_21)
  expect_equal(B[3,3,,], B_22)
  expect_equal(B[4,4,,], B_33)
  expect_equal(B[4,5,,], B_34)
  expect_equal(B[5,5,,], B_44)
  expect_equal(B[5,6,,], B_45)
  expect_equal(B[6,6,,], B_55)
  expect_equal(B[6,7,,], B_56)

  ##### hard coded matrices smaller case test
  K <- 2
  r <- 2
  s <- 1

  B_00 = matrix(c(lam_P,0,0,
                  0,lam_P,0,
                  0,0,lam)
                , nrow = 3, ncol = 3, byrow = TRUE)
  B_01 = matrix(c(lam_P,0,0, 0, 0,
                  0,lam_P,0,  0, 0,
                  0,0,lam_P,  0, 0,
                  0,0,0,lam_P, 0,
                  0,0,0,0,lam)
                , nrow = 3, ncol = 3, byrow = TRUE)
  B_10 = matrix(c(-(lam+(0*mu_P)),lam_A,0, 0, 0,
                  mu_A,-(lam +calc_mu_mn(0,1, mu_P,mu_A)),lam_A,  0, 0,
                  0,mu_A,-(lam +calc_mu_mn(0,1, mu_P,mu_A)),  lam_A, 0,
                  0,0, mu_A,-(lam +calc_mu_mn(0,1, mu_P,mu_A)),  lam_A,
                  0,0,0,mu_A,-(lam +calc_mu_mn(0,1, mu_P,mu_A)))
                , nrow = 3, ncol = 3, byrow = TRUE)
  B_11 = matrix(c(-(lam+(1*mu_P)),lam_A,0, 0, 0,
                  mu_A,-(lam +calc_mu_mn(1,1, mu_P,mu_A)),lam_A,  0, 0,
                  0,mu_A,-(lam +calc_mu_mn(1,1, mu_P,mu_A)),  lam_A, 0,
                  0,0, mu_A,-(lam +calc_mu_mn(1,1, mu_P,mu_A)),  lam_A,
                  0,0,0,mu_A,-(lam +calc_mu_mn(1,1, mu_P,mu_A)))
                , nrow = 3, ncol = 3, byrow = TRUE)
  B_21 = matrix(c(mu_P,0,0, 0, 0,
                  0,mu_P,0,  0, 0,
                  0,0,mu_P,  0, 0,
                  0,0,0,mu_P, 0,
                  0,0,0,0,mu_P)
                , nrow = 3, ncol = 3, byrow = TRUE)
  B_22 = matrix(c(2*mu_P,0,0, 0, 0,
                  0,2*mu_P,0,  0, 0,
                  0,0,2*mu_P,  0, 0,
                  0,0,0,2*mu_P, 0,
                  0,0,0,0, 2*mu_P)
                , nrow = 3, ncol = 3, byrow = TRUE)
  B_32 = matrix(c(0,3*mu_P*q,0,0, 0,
                  0,0,2*mu_P*calc_alpha_i(p,1),0, 0,
                  0,0,0, 2*mu_P*calc_alpha_i(p,1), 0,
                  0,0,0,  0, 2*mu_P*q,
                  0,0,0,0, 0)
                , nrow = 3, ncol = 3, byrow = TRUE)

  B <- Calc_Bmn(K,s,r,lam,lam_A,lam_P,mu_P,mu_A, p)

  #Check values -- Note: R indexes at 1
  expect_equal(B[1,1,,], B_00)
  expect_equal(B[1,2,,], B_01)
  expect_equal(B[2,1,,], B_10)
  expect_equal(B[2,2,,], B_11)
  expect_equal(B[3,2,,], B_21)
  expect_equal(B[3,3,,], B_22)
  expect_equal(B[4,3,,], B_32)

})
