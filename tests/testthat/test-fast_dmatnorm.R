test_that("dmatnorm works", {
  require(matrixNormal)
  A <- datasets::CO2[1:10, 4:5]
  #A <- matrix(rnorm(10000*2, 5, 2), nrow = 10000, ncol = 2)
  M <- cbind(stats::rnorm(10, 435, 296), stats::rnorm(10, 27, 11))
  #V <- matrix(rnorm(1000*1000))
  V <- matrix(c(87, 13, 13, 112), nrow = 2, ncol = 2, byrow = TRUE)
  #' V # Right covariance matrix (2 x 2), say the covariance between parameters.
  U <- diag(10) # Block of left-covariance matrix, say the covariance between subjects.

  # check output matches below
  # have to import bc functions below have checks removed
  expect_equal(fast_dmatnorm(A, M, U, V), dmatnorm(A, M, U, V))
  expect_equal(fast_dmatnorm(A, M, U, V, log = FALSE), dmatnorm(A, M, U, V, log = FALSE))
})
