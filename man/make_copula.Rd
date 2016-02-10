\name{make_copula}
\alias{make_copula}
\title{
  Copula transformation function  
}
\description{
  This functions transforms a variable (condition: non-normally distributed) to its copula-correction term which can be entered in model estimation. 
  See Park and Gupta (2012, Marketing Science)
}

\references{
  Park, Sungho, and Sachin Gupta. "Handling endogenous regressors by joint estimation using copulas." Marketing Science 31.4 (2012): 567-586.
  Hannes Datta, Bram Foubert, and Harald J. Van Heerde (2015) The Challenge of Retaining Customers Acquired with Free Trials. Journal of Marketing Research: April 2015, Vol. 52, No. 2, pp. 217-234.
}

\examples{
  X=runif(1000)
  X_star = make_copula(X)
}
