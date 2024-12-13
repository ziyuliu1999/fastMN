
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastMN

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/ziyuliu1999/fastMN/graph/badge.svg)](https://app.codecov.io/gh/ziyuliu1999/fastMN)
<!-- badges: end --> <!-- badges: start -->
[![R-CMD-check](https://github.com/ziyuliu1999/fastMN/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ziyuliu1999/fastMN/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of fastMN is to efficiently compute the PDF and CDF of the
matrix normal distribution and to rapidly generate random samples from
it.

## Installation

You can install the development version of fastMN from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ziyuliu1999/fastMN")
```

You may need to add in your personal github token.

## Functions

Three main functions introduced in the fastMN package

`fast_rmatnorm`: Sample random variables following matrix normal
distributions

`fast_dmatnorm`: Calculate PDF of the matrix normal distribution

`fast_pnormmat`: Calculate CDF of the matrix normal distribution

## Example

Followings aer basic example which shows you how to use the functions in
general:

*Example for the fast_rmatnorm*

``` r
library(fastMN)
set.seed(1)
n = 100
p = 100
U <- matrix(0.5,nrow=n,ncol=n) + 0.5*diag(n)
V <- matrix(0.8,nrow=p,ncol=p) + 0.2*diag(p)
Y <- fast_rmatnorm(n = n, p = p,U_cov = U, V_cov=V)
```

*Example for the fast_dmatnorm*

``` r
n = 100
p = 100
U <- matrix(0.5,nrow=n,ncol=n) + 0.5*diag(n)
V <- matrix(0.8,nrow=p,ncol=p) + 0.2*diag(p)
Y<- fast_rmatnorm(n = n, p = p, U_cov = U, V_cov=V)
M = matrix(0, nrow = n, ncol = p)
p_val = fast_dmatnorm(Y, M, U, V)
```

*Example for the fast_pnormmat*

``` r
n <- 3
p <- 2
# Create the mean matrix
M <- matrix(0, nrow = n, ncol = p)
matrixU <- matrix(c(1, 2, 0, 2, -1, -2, 1, 3, 0), nrow = n, ncol = n)
matrixV <- matrix(c(2, 0, 1, 4), nrow = p, ncol = p)
U_cov <- crossprod(matrixU)  # Row covariance matrix
V_cov <- crossprod(matrixV)  # Column covariance matrix
# Set the lower and upper bounds for the integration
Lower <- matrix(-10, nrow = n, ncol = p)
Upper <- matrix(10, nrow = n, ncol = p)
#'
# Approximate the CDF using the naive Monte Carlo method
fast_pnormmat(Lower = Lower, Upper = Upper, M = M, U_cov = U_cov, V_cov = V_cov, method = "naive_monte_carlo", N = 1000)
#>              method   cdf   log_cdf
#> 1 naive monte carlo 0.287 -1.248273
```

``` r
# Compute the CDF using the mvnorm::pmvnorm
fast_pnormmat(Lower = Lower, Upper = Upper, M = M, U_cov = U_cov, V_cov = V_cov, method = "pmvnorm")
#>               method       cdf log_cdf
#> 1 mvnorm computation 0.2653779 -1.3266
```
