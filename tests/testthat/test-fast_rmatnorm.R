test_that("sample works", {
  MSE = function(true_val, est_val) {
    mean((true_val - est_val)^2)
  }
  res_gen = function(num_samp, true_M, true_U, true_V, multi_fast) {
    n = nrow(true_U)
    p = nrow(true_V)
    sampled_mean = apply(multi_fast, c(1, 2), mean)
    mse_mean = MSE(true_M, sampled_mean)
    flattened_samples = matrix(NA, nrow = num_samp, ncol = n * p)
    for (i in 1:num_samp) {
      flattened_samples[i, ] = as.vector(multi_fast[,,i])
    }
    empirical_cov = cov(flattened_samples)
    theoretical_cov = kronecker(V, U)
    mse_cov = MSE(empirical_cov, theoretical_cov)
    return(list(mse_mean = mse_mean, mse_cov = mse_cov))
  }
  M <- cbind(stats::rnorm(10, 435, 296), stats::rnorm(10, 27, 11))  # Mean matrix
  V <- matrix(c(87, 13, 13, 112), nrow = 2, ncol = 2, byrow = TRUE)  # Right covariance
  U <- diag(10)  # Left covariance
  n <- nrow(M)  # Rows of M
  p <- ncol(M)  # Columns of M
  multi_fast_10000 = fast_rmatnorm(num_samp = 10000, n = n, p = p, M = M, U_cov = U, V_cov = V)
  mse_res_10000_cov = res_gen(10000, M, U, V, multi_fast_10000)
  expect_equal(mse_res_10000_cov$mse_mean < 0.1, TRUE)
  expect_equal(mse_res_10000_cov$mse_cov < 2, TRUE)
  U_prec = solve(U)
  V_prec = solve(V)
  multi_fast_10000_prec = fast_rmatnorm(num_samp = 10000, n = n, p = p, M = M, U_prec = U_prec, V_prec = V_prec, useCov = FALSE)
  mse_res_10000_prec = res_gen(10000, M, U, V, multi_fast_10000_prec)
  expect_equal(mse_res_10000_prec$mse_mean < 0.1, TRUE)
  expect_equal(mse_res_10000_prec$mse_cov < 2, TRUE)
})
