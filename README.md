
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Overview

<!-- badges: start -->

<!-- badges: end -->

User-friendly and flexible algorithm for time series **segmentation**
with a Regression model with a Hidden Logistic Process (RHLP).

## Installation

You can install the development version of RHLP from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fchamroukhi/RHLP")
```

To build *vignettes* for examples of usage, type the command below
instead:

``` r
# install.packages("devtools")
devtools::install_github("fchamroukhi/RHLP", 
                         build_opts = c("--no-resave-data", "--no-manual"), 
                         build_vignettes = TRUE)
```

Use the following command to display vignettes:

``` r
browseVignettes("RHLP")
```

## Usage

``` r
library(RHLP)

data("simulatedtimeserie")

K <- 5 # number of regimes (mixture components)
p <- 3 # dimension of beta (order of the polynomial regressors)
q <- 1 # dimension of w (order of the logistic regression: to be set to 1 for segmentation)
variance_type <- variance_types$heteroskedastic

n_tries <- 1
max_iter = 1500
threshold <- 1e-6
verbose <- TRUE
verbose_IRLS <- FALSE

solution <- emRHLP(simulatedtimeserie$X, t(simulatedtimeserie$Y), K, p, q, 
                   variance_type, n_tries, max_iter, threshold, verbose, verbose_IRLS)

solution$plot()
```

<img src="man/figures/README-unnamed-chunk-5-1.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-2.png" style="display: block; margin: auto;" />
