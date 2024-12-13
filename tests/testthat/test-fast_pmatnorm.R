test_that("CDF works", {
  require(matrixNormal)
  n <- 3
  p <- 2
  #' # Create the mean matrix
  M <- matrix(0, nrow = n, ncol = p)
  matrixU <- matrix(c(1, 2, 0, 2, -1, -2, 1, 3, 0), nrow = n, ncol = n)
  matrixV <- matrix(c(2, 0, 1, 4), nrow = p, ncol = p)
  U_cov <- crossprod(matrixU)  # Row covariance matrix
  V_cov <- crossprod(matrixV)  # Column covariance matrix
  Lower <- matrix(-10, nrow = n, ncol = p)
  Upper <- matrix(10, nrow = n, ncol = p)
  fast_pnormmat(Lower = Lower, Upper = Upper, M = M, U_cov = U_cov, V_cov = V_cov, method = "naive_monte_carlo", N = 1000)
})
