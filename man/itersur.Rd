\name{itersur}
\alias{itersur}
\title{
  Iterative Seamingly-Unrelated-Regression Estimation
}
\description{
  Estimates a SUR model on a possibly unbalanced panel, as layed out in Zellner (1962) and Bronnenberg, Mahajan, and Vanhonacker (2000).

}

\references{
  Bronnenberg, B. J., Mahajan, V., & Vanhonacker, W. R. (2000). The emergence of market structure in new repeat-purchase categories: The interplay of market share and retailer distribution. Journal of Marketing Research, 37(1), 16-31.
  Zellner, A. (1962). An efficient method of estimating seemingly unrelated regressions and tests for aggregation bias. Journal of the American statistical Association, 57(298), 348-368.
}

\examples{
  X=runif(1000)
  X_star = make_copula(X)
}
