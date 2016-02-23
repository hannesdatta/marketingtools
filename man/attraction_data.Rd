\name{attraction_data}
\title{Prepare data to estimate a market share attraction model}
\alias{attraction_data}
\concept{ LMM }
\description{
  Prepare a data set (in long format) to its base-brand representation, so that it can be estimated using SUR or maximum likelihood.
}
\usage{
attraction_data(formula, data = NULL, subset = NULL, heterogenous = NULL, index = NULL, model = "MCI", benchmark = NULL, \dots)
}
\arguments{
  \item{formula}{a two-sided linear formula object describing the attraction model to be estimated, with the market share on the left of a \code{~} operator (must sum to 1 per period) and the explanatory variables, separated
    by \code{+} operators, on the right. All explanatory variables are estimated as homogenous, unless they are added as one-sided formulas under \code{heterogenous}.}

  \item{data}{an optional data frame containing the variables named in
    \code{formula}, \code{heterogenous}, and \code{index}.  By default the variables are taken from the environment from which \code{attraction_data} is called.}

  \item{heterogenous}{a one-sided formula object, with variables on the right of a \code{~} sign describing all heterogenous variables, separated by \code{+} operators.}

  \item{index}{a one-sided formula object, with two variables on the right of a \code{~} sign, describing (1) the brand index, and (2) the time index of the data, seperated by a \code{+} sign. The package authors strongly recommend to sort the data by brand and time before calling \code{attraction_date}.}
  
  \item{model}{specifies the transformation type that will be applied to the data. Can be one of \code{"MCI"} (standard) or \code{"MNL".}. For details, see Fok 2001.}

  \item{benchmark}{character vector, specifying the benchmark brand to be chosen for normalization; default of \code{NULL} implies that the algorithm automatically determines the benchmark brand (chosing the brand with most available observations for all variables and all time periods).
  }
  \item{subset}{an optional expression indicating the subset of the rows
    of \code{data} that should be used in the fit. This can be a logical
    vector, or a numeric vector indicating which observation numbers are
    to be included, or a character vector of the row names to be
    included. All observations are included by default.}

  \item{\dots}{other potential arguments. Currently not implemented.}
}

\value{
  An object of class \code{\link[=XXXmerMod-class]{attr.data}}, for which some methods
  are available (e.g., show).}
}

\details{
  \itemize{
    \item{...}
      }
}
\seealso{
  \code{\link[marketingtools]{itersur}} for details on how to estimate the model.
}
\examples{
# Load package
  require(marketingtools) 

# Simulate raw data
  rawdata <- do.call('cbind', attraction_simulate_data())
  head(rawdata)

# Transform to base-brand representation
  dtbb <- attraction_data(formula = y ~ X.var_3 + X.var_4 + X.var_5 + X.var_1 + X.var_2, heterogenous = ~ X.var_1 + X.var_2, index = ~ index.var + index.t, data=rawdata, model="MCI")

# verify object
  validObject(dt)

# show a summary about the transformed data
  show(dt)

# Estimate model using itersur (FGLS)
  m <- itersur(X=dtbb@X,Y=as.matrix(dtbb@y), index=data.frame(date=dtbb@period,brand=dtbb@individ), method = "FGLS")
  
  show(m) # summarize m
  print(m@sigma) # show contemporaneous variance-covariance matrix
  
# Estimate model using itersur with correction for autocorrelation in the residuals (FGLS-Praise-Winsten)
  m <- itersur(X=dtbb@X,Y=as.matrix(dtbb@y), index=data.frame(date=dtbb@period,brand=dtbb@individ), method = "FGLS-Praise-Winsten")
  
  show(m) # summarize m
  print(m@sigma) # show contemporaneous variance-covariance matrix
    

}
\keyword{models}
\references{
   Cooper, L. G., & Nakanishi, M. (1988). Market-share analysis: Evaluating competitive marketing
   effectiveness. Boston: Kluwer Academic Publishers.
   
   Fok, D., Franses, P. H., & Paap, R. (2001). Econometric analysis of the market share 
   attraction model.

}
}