#' Random generation for the matrix normal distribution accepting either covariance or precision
#'
#' @param num_samp The number of observations to be sampled
#' @param n The row dimension of the sampled matrix normal variable
#' @param p The column dimension of the sampled matrix normal variable
#' @param M The Mean matrix
#' @param U_cov The row covariance matrix with dimension n * n
#' @param V_cov The column covariance matrix with dimension p * p
#' @param U_prec The row precision matrix with dimension n * n
#' @param V_prec The column precision matrix with dimension p * p
#' @param useCov logical; if TRUE the U_cov and V_cov will be used otherwise U_prec and V_prec will be used
#' @param small_tol a small number with default number 1e-6 which is used if the the full rank precision issue encountered
#' @return If num_samp = 1, it will return a matrix. Otherwise a 3-D array will be returned
#' @examples
#' n = 100
#' p = 100
#' U <- matrix(0.5,nrow=n,ncol=n) + 0.5*diag(n)
#' V <- matrix(0.8,nrow=p,ncol=p) + 0.2*diag(p)
#' Y<- sample_matrix_normal(U_cov = U, V_cov=V)
#' @export
fast_rmatnorm <- function(num_samp = 1,
                          n = 10,
                          p = 5,
                          M = NULL,
                          U_cov = NULL,
                          V_cov = NULL,
                          U_prec = NULL,
                          V_prec = NULL,
                          useCov = TRUE, small_tol = 1e-6) {
  if (!is.null(M)) {
    n <- nrow(M)
    p <- ncol(M)
  }
  if (is.null(M)) {
    M <- matrix(0, nrow = n, ncol = p)
  }

  if (useCov) {
    # Default covariance matrices
    if (is.null(U_cov)) {
      U_cov <- diag(n)
    }
    if (is.null(V_cov)) {
      V_cov <- diag(p)
    }

    # Cholesky decomposition of covariance matrices
    Ru <- chol(U_cov)
    Rv <- chol(V_cov)
    if (num_samp == 1) {
      Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
      Y <- M + crossprod(Ru, Z) %*% Rv
      return(Y)
    } else {
      return_res <- array(NA, dim = c(n, p, num_samp))
      Z <- array(rnorm(n * p * num_samp), dim = c(n, p, num_samp))
      for (i in 1:num_samp) {
        return_res[,,i] <- M + crossprod(Ru, Z[,,i]) %*% Rv
      }
      return(return_res)
    }

  } else {
    # Default precision matrices
    if (is.null(U_prec)) {
      U_prec <- diag(n)
    }
    if (is.null(V_prec)) {
      V_prec <- diag(p)
    }

    # Cholesky decomposition and inversion of precision matrices
    Ru <- chol(U_prec)
    Rv <- chol(V_prec)
    Ru_inv <- backsolve(Ru, diag(nrow(Ru)))
    Rv_inv <- backsolve(Rv, diag(nrow(Rv)))
    if (num_samp == 1) {
      Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
      Y <- M + tcrossprod(Ru_inv %*% Z, Rv_inv)
      return(Y)
    } else {
      return_res <- array(NA, dim = c(n, p, num_samp))
      Z <- array(rnorm(n * p * num_samp), dim = c(n, p, num_samp))
      for (i in 1:num_samp) {
        return_res[,,i] <- M + tcrossprod(Ru_inv %*% Z[,,i], Rv_inv)
      }
      return(return_res)
    }
  }
}
