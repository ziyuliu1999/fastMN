## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
## usethis namespace: start
#' @useDynLib fastMN, .registration = TRUE
## usethis namespace: end
NULL

dim_check <- function(M, U, V) {
  if(!is.null(M)){
    n <- nrow(M)
    p <- ncol(M)
  }
  if(!is.null(U)){
    n <- ncol(U)
  }
  if(!is.null(V)){
    p <- ncol(V)
  }
  if(is.null(U)){
    U <- diag(n)
  }
  if(is.null(V)){
    V <- diag(p)
  }
  if(is.null(M)){
    M = matrix(0,nrow=n,ncol=p)
  }
  return(list(M,U,V, n, p))
}


is.sym <- function(X, tol=1e-8) {
  if (dim(X)[1] != dim(X)[2]) {FALSE}

  else if (all(abs(X-t(X)) < tol)) {TRUE}
  else {FALSE}
}



is.pos.def <- function(X, tol=1e-8) {
  eigenvalues <- eigen(X, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) > tol) {TRUE}
  else {FALSE}
}

# MATRIX CHECK
check_matnorm <- function(Z = NULL, # For sampling function, there's no Z to check
                          M,
                          U,
                          V,
                          tol=1e-8) {
  if (!is.null(Z) & anyNA(Z)) {
    stop("Z contains missing values.", call. = FALSE)
  }
  if (anyNA(M)) {
    stop("M contains missing values.", call. = FALSE)
  }
  if (anyNA(U)) {
    stop("U contains missing values.")
  }
  if (anyNA(V)) {
    stop("V contains missing values.")
  }
  if (nrow(M) != nrow(U)) {
    stop("The mean matrix M has different sample size than the scale sample size
         matrix U. M has ", dim(M)[[1]], "rows, and U has ", dim(U)[[1]], ".")
  }
  if (ncol(M) != nrow(V)) {
    stop("The mean matrix M has different number of parameters than scale
         parameter matrix V: M  -- ", dim(M)[2], "; V -- ", dim(V)[1], ".")
  }
  if (!is.sym(U, tol)) {
    stop("U is not symmetric.")
  }
  if (!is.sym(V, tol)) {
    stop("V is not symmetric.")
  }

  if (!is.pos.def(U, tol)) {
    stop("U is not positive definite. Calculation may not be accurate.
         Possibly lower tolerance.")
  }
  if (!is.pos.def(V, tol)) {
    stop("V is not positive definite. Calculation may not be accurate.
         Possibly lower tolerance.")
  }
  return(invisible())
}

# function checking Lower and Upper for vectorization
vectorized_check <- function(Lower, Upper, n, p) {
  Lower <- if (is.matrix(Lower)) {
    if (dim(Lower)[1] != n | dim(Lower)[2] != p) {
      stop("Lower must have dimensions ", n, "x", p, ".")
    }
    as.vector(Lower)
  } else if (Lower == -Inf) {
    rep(-Inf, n * p)
  } else {
    stop("Lower must be a numeric matrix or -Inf.")
  }

  Upper <- if (is.matrix(Upper)) {
    if (dim(Upper)[1] != n | dim(Upper)[2] != p) {
      stop("Upper must have dimensions ", n, "x", p, ".")
    }
    as.vector(Upper)
  } else if (Upper == Inf) {
    rep(Inf, n * p)
  } else {
    stop("Upper must be a numeric matrix or Inf.")
  }

  return(list(Lower = Lower, Upper = Upper))
}


#' Random generation for the matrix normal distribution accepting both covariance and precision
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
#' Y<- fast_rmatnorm(U_cov = U, V_cov=V)
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
    if (!is.pos.def(U_cov)) {
      U_cov = U_cov + diag(small_tol)
    }
    if (!is.pos.def(V_cov)) {
      V_cov = V_cov + diag(small_tol)
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
    if (!is.pos.def(U_prec)) {
      U_prec = U_prec + diag(small_tol)
    }
    if (!is.pos.def(V_prec)) {
      V_prec = V_prec + diag(small_tol)
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






#' Calculate PDF for the matrix normal distribution accepting both covariance and precision
#'
#' @param Z The number of observations to be sampled
#' @param M The Mean matrix
#' @param U The row covariance/precision matrix with dimension n * n
#' @param V The column covariance/precision matrix with dimension p * p
#' @param log logical; if TRUE, probabilities p are given as log(p)
#' @param Precision logical; if TRUE U and V will be considered as precision
#' @param tol a small number with default number 1e-6 which is used if the the full rank precision issue encountered
#' @return Return the pdf of Z
#' @examples
#' n = 100
#' p = 100
#' U <- matrix(0.5,nrow=n,ncol=n) + 0.5*diag(n)
#' V <- matrix(0.8,nrow=p,ncol=p) + 0.2*diag(p)
#' Y<- sample_matrix_normal(U_cov = U, V_cov=V)
#' M = matrix(0, nrow = n, ncol = p)
#' p_val = fast_dmatnorm(Y, M, U, V)
#' @export

fast_dmatnorm <- function(Z, M, U, V, log=TRUE, Precision=FALSE, tol=1e-8) {
  dc = dim_check(M,U,V)
  M <- dc[[1]]
  U <- dc[[2]]
  V <- dc[[3]]
  n <- dc[[4]]
  p <- dc[[5]]
  # U,v square
  #  M nxp
  # cholesky check - pos definite
  check_matnorm(Z,M,U,V,tol)

  Ru <- chol(U)
  Rv <- chol(V)

  # solving for A
  D <- Z - M

  # calculate log density first
  # RCPP for sum(A^2)
  if (!Precision) {
    A2 = forwardsolve(Ru, D, upper.tri=TRUE,transpose=TRUE)
    A = t(forwardsolve(t(Rv), t(A2)))
    # use Rcpp function
    numerator = -0.5 * sum_of_squares(A)#sum(A^2)
    denom = (n*p/2)*log(2*pi) + sum(n*log(diag(Rv))) + sum(p*log(diag(Ru)))
  } else {
    B = Ru %*% D %*% t(Rv)
    # use Rcpp function
    numerator = -0.5 * sum_of_squares(B)#sum(B^2)
    denom = (n*p/2)*log(2*pi) - sum(n*log(diag(Rv))) - sum(p*log(diag(Ru)))
  }

  log.dens = numerator - denom

  if (log) {
    return(log.dens)
  } else {
    return(exp(log.dens))
  }
}

#' Approximate the CDF of a Matrix Normal Distribution
#'
#' This function calculates the cumulative distribution function (CDF) of a matrix normal distribution
#' It supports three methods: naive Monte Carlo approximation, Sobol sequence-based quasi-Monte Carlo, and direct computation using mvtnorm::pmvnorm.
#'
#' @param Lower A matrix of lower bounds for the integration region. Defaults to \code{-Inf}.
#' @param Upper A matrix of upper bounds for the integration region. Defaults to \code{Inf}.
#' @param M The mean of the matrix normal distribution.
#' @param U_cov The row covariance matrix. Used when \code{useCov = TRUE}.
#' @param V_cov The column covariance matrix. Used when \code{useCov = TRUE}.
#' @param U_prec A matrix specifying the row precision matrix (inverse of covariance). Used when \code{useCov = FALSE}.
#' @param V_prec A matrix specifying the column precision matrix (inverse of covariance). Used when \code{useCov = FALSE}.
#' @param useCov A logic indicating whether to use covariance matrices (\code{TRUE}) or precision matrices (\code{FALSE}). Defaults to \code{TRUE}.
#' @param method A character string specifying the method to use for the CDF approximation. Defaults to \code{"naive_monte_carlo"}. Options are:
#' \itemize{
#'   \item \code{"naive_monte_carlo"}: Uses naive Monte Carlo integration.
#'   \item \code{"sobol"}: Uses Sobol sequence-based quasi-Monte Carlo integration (requires the `randtoolbox` package).
#'   \item \code{"pmvnorm"}: Uses the `mvtnorm::pmvnorm` function for direct computation.
#' }
#' @param N An integer specifying the number of samples to draw for Monte Carlo or Sobol methods. Defaults to \code{NULL}, in which case an adaptive sample size is used.
#' @param tol A numeric value specifying the tolerance for precision checks in covariance/precision matrices. Defaults to \code{1e-8}.
#' @param max_iter An integer specifying the maximum number of iterations for iterative methods. Defaults to \code{1000}.
#' @param algorithm A function specifying the integration algorithm to use with mvtnorm::pmvnorm. Defaults to \code{mvtnorm::GenzBretz()}.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{method}: The method used for the computation (e.g., \code{"naive monte carlo"}, \code{"sobol monte carlo"}, \code{"mvnorm computation"}).
#'   \item \code{cdf}: The estimated CDF value.
#'   \item \code{log_cdf}: The log of the CDF value.
#' }
#'
#' @examples
#' # Define the dimensions and M, U, V
#' n <- 3
#' p <- 2
#' # Create the mean matrix
#' M <- matrix(0, nrow = n, ncol = p)
#' matrixU <- matrix(c(1, 2, 0, 2, -1, -2, 1, 3, 0), nrow = n, ncol = n)
#' matrixV <- matrix(c(2, 0, 1, 4), nrow = p, ncol = p)
#' U_cov <- crossprod(matrixU)  # Row covariance matrix
#' V_cov <- crossprod(matrixV)  # Column covariance matrix
#'
#' # Set the lower and upper bounds for the integration
#' Lower <- matrix(-10, nrow = n, ncol = p)
#' Upper <- matrix(10, nrow = n, ncol = p)
#'
#' # Compute the CDF using the naive Monte Carlo method
#' fast_pnormmat(Lower = Lower, Upper = Upper, M = M,
#'               U_cov = U_cov, V_cov = V_cov, method = "naive_monte_carlo", N = 1000)
#'
#' # Compute the CDF using the Sobol method
#' fast_pnormmat(Lower = Lower, Upper = Upper, M = M,
#'               U_cov = U_cov, V_cov = V_cov, method = "sobol", N = 1000)
#'
#' # Compute the CDF using the mvtnorm package
#' fast_pnormmat(Lower = Lower, Upper = Upper, M = M,
#'               U_cov = U_cov, V_cov = V_cov, method = "pmvnorm")
#'
#' @export
fast_pnormmat <-function(Lower = -Inf, # Lower bound matrix
                         Upper = Inf, # Upper bound matrix
                         M = NULL, # Mean matrix
                         U_cov = NULL,
                         V_cov = NULL,
                         U_prec = NULL,
                         V_prec = NULL,
                         useCov = TRUE,
                         method = "naive_monte_carlo", # naive_monte_carlo, sobol, pmvnorm
                         N = NULL, # optional number of samples for monte carlo
                         tol=1e-8,
                         max_iter = 1000,
                         algorithm = mvtnorm::GenzBretz() # default algorithm for mvtnorm::pmvnorm
) {

  # Matrix dimension checks
  dc <- if (useCov) {
    dim_check(M = M, U = U_cov, V = V_cov)
  } else {
    dim_check(M = M, U = U_prec, V = V_prec)
  }
  M <- dc[[1]]; n <- dc[[4]]; p <- dc[[5]]

  # Check Lower and Upper bounds using vectorized_check
  bounds <- vectorized_check(Lower, Upper, n, p)
  Lower <- bounds$Lower
  Upper <- bounds$Upper

  # Validate covariance or precision matrices
  if (useCov) {
    check_matnorm(Z = NULL, M = M, U = U_cov, V = V_cov, tol = tol)
    U <- U_cov; V <- V_cov
  } else {
    check_matnorm(Z = NULL, M = M, U = U_prec, V = V_prec, tol = tol)
    U <- solve(U_prec); V <- solve(V_prec)
  }

  if (method == "pmvnorm") {
    # Vectorized method using mvtnorm
    cdf <- mvtnorm::pmvnorm(
      lower = Lower,
      upper = Upper,
      mean = as.vector(M),
      sigma = kronecker(V, U),
      algorithm = algorithm
    )
    method <- "mvnorm computation"
  } else if (method == "sobol") {

    if (!requireNamespace("randtoolbox", quietly = TRUE)) {
      stop("The 'randtoolbox' package is required for the Sobol method. Please install it using install.packages('randtoolbox').")
    }
    if (is.null(N)) {
      N <- max(2000, 10*n*p)  # Adaptive sample size
    }
    sobol_points <- randtoolbox::sobol(n = N, dim = n * p)

    # Transform Sobol to normal distribution (vectorized) because that's how sobol works
    Z <- qnorm(sobol_points)

    # Reshape Z to n x p x N
    Z <- array(Z, dim = c(n, p, N))

    Ru <- chol(U)
    Rv <- chol(V)
    samples <- array(0, dim = c(n, p, N))
    for (i in 1:N) { #back in matrix form
      samples[,,i] <- M + crossprod(Ru, Z[,,i]) %*% Rv
    }

    # Check bounds
    within_bounds <- sweep(samples, c(1, 2), Lower, `>=`) & sweep(samples, c(1, 2), Upper, `<=`)
    valid_samples <- apply(within_bounds, 3, all)
    count <- sum(valid_samples)
    cdf <- count / N
    method <- "sobol monte carlo"

  } else if (method == "naive_monte_carlo"){
    # Naive monte Carlo method
    if (is.null(N)) {
      N <- max(2000, 10*n*p)  # Adaptive sample size
    }
    # Use our sampler
    samples <- fast_rmatnorm(num_samp = N,
                             n = n,
                             p = p,
                             M = M,
                             U_cov = U,
                             V_cov = V,
                             U_prec = U,
                             V_prec = V,
                             useCov = useCov)
    within_bounds <- sweep(samples, c(1, 2), Lower, `>=`) & sweep(samples, c(1, 2), Upper, `<=`)
    valid_samples <- apply(within_bounds, 3, all)
    count <- sum(valid_samples)
    cdf <- count / N
    method <- "naive monte carlo"
  } else {
    stop("Invalid method. Please choose one of 'naive_monte_carlo', 'sobol', or 'pmvnorm'.")
  }

  df <- data.frame(method = method, cdf = cdf, log_cdf = log(cdf))
  return(df)
}

