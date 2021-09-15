# marketingtools

<!-- badges: start -->
<!-- badges: end -->

The goal of `marketingtools` is to collect and make accessible code for conducting empirical (marketing) research.

## Features

- Computation of Gaussian Copula correction terms, using their ECDF (Park & Gupta 2012) (`make_copula`)
- Estimation of Feasible Generalized Least Squares (FGLS) with unbalanced panels (`itersur`)
- Estimation of market share attraction models (Fok 2001)
  o Transformation of wide panels into their base-brand representation
  o Support for homogenous and heterogenous coefficients
  o Restricted competition; fully extended model not implemented
  o Estimation via `itersur()` (FGLS)
  
## Installation

You can install the released version of musiclabels from GitHub with:

``` r
install.packages("devtools")
devtools::install_github("hannesdatta/marketingtools")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(musicMetadata)

### to be done.

```

## Contributions

For a general guide on how to contribute to this repository/package, see https://tilburgsciencehub.com/learn/git-collaborate.

Curious to see what you can do? Head over to the repository's Issue page to find out more about it.
