test_that("multiplication works", {
  expect_that( object = is.numeric(temp_C), condition = equals(TRUE) )
  expect_equal(2 * 2, 4)
})
