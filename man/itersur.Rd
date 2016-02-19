\name{itersur}
\alias{itersur}
\title{
  Seamingly Unrelated Regression (SUR) Estimation using FGLS
}
\description{
  Estimates a SUR model on a possibly unbalanced panel, as layed out in Zellner (1962) and Bronnenberg, Mahajan, and Vanhonacker (2000). Supports Praise-Winsten correction for auto-correlation in the residuals.

}
\usage{
itersur(X, Y, index, model="FGLS", maxiter=1000)
}
\arguments{
  \item{X}{matrix of explanatory variables, possibly constructed using \code{attraction_data}.}
  \item{Y}{vector of dependent variable, possibly constructed using \code{attraction_data}.}
  \item{index}{\code{data.frame} with two columns, describing the index of the data set. The first column needs to list an indicator for time (e.g., weeks); the second column lists an indicator for cross-sectional units (e.g., factor, character, or numeric). The columns need not be named. Index can be constructed using \code{attraction_data}.}
  \item{model}{Choice of \code{"FGLS"} (for feasible generalized least squares), or \code{"FGLS-Praise-Winsten"} for feasible generalized least squares with Praise-Winsten transformation of the dependent and explanatory variables (to correct for autocorrelation in the residuals) (see van Heerde, Leeflang, and Wittink 2004; Greene 2003.}
  \item{maxiter}{Maximum number of FGLS iterations (set to 1000 by default)}
  }
\value{
\code{itersur} returns on object of class "itersur".

The functions \code{show} are used to obtain and print a summary of the results. The generic accessor function \code{coefficients} extracts the coefficients as a \code{data.frame}.

An object of class \code{"itersur"} is an S4-class containing at least the following components, which can be accessed using the \code{@}-operator:

  \item{coefficients}{a \code{data.frame} with coefficients}
  \item{X}{the predictor matrix to estimate the final iteration of the model}
  \item{Y}{the dependent variable to estimate the final iteration of the model (unchanged for FGLS, transformed for Praise-Winsten.)}
  \item{index}{The indexing \code{data.frame}, with column 1 being the indicators for time periods, and column 2 being the indicators for cross-sectional units.}
  \item{iterations}{The number of iterations of the algorithm until convergence.}
  \item{bic}{The BIC for the model.}
  \item{predicted}{the fitted values}
  \item{resid}{the residuals, that is response minus fitted values.}
  \item{varcovar}{The variance-covariance matrix of the estimates.}
  \item{method}{The method used to estimat the model; one of code{"FGLS"} or \code{"FGLS-Praise-Winsten"}.}
  \item{rho}{A numeric vector, indicating the auto-correlation of residuals for all cross-sectional units.}

}
\references{
  Bronnenberg, B. J., Mahajan, V., & Vanhonacker, W. R. (2000). The Emergence of Market Structure in New Repeat-purchase Categories: The Interplay of Market Share and Retailer Distribution. Journal of  Marketing Research, 37(1), 16-31.
  
  Greene, W. H. (2003). Econometric Analysis. Pearson Education.

  van Heerde, H. J., Leeflang, P. S., & Wittink, D. R. (2004). Decomposing the Sales Promotion Bump with Store Data. Marketing Science, 23(3), 317-334.
  
  Zellner, A. (1962). An Efficient Method of Estimating seemingly Unrelated Regressions and Tests for Aggregation Bias. Journal of the American Statistical Association, 57(298), 348-368.
}

\examples{
See examples listed under \link{attraction_data}.

}
