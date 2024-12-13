
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastMN

<!-- badges: start -->

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

## Functions

Three main functions introduced in the fastMN package

`fast_rmatnorm`: Sample random variables following matrix normal
distributions

`fast_dmatnorm`: Calculate PDF of the matrix normal distribution

`fast_rmatnorm`: Calculate CDF of the matrix normal distribution

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(fastMN)
set.seed(1)
```
