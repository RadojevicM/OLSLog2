\name{ols}
\alias{ols}

\title{
  linear OLS-regression
}
\description{
  The function estimates a linear model by using Ordinary-Least-Squares regression. 
  OLS-coefficients and common summary statistics of the estimation are provided.
}
\usage{
  ols(X, Y, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
    A data matrix with a column of 1 for each observation and one column for each independent variable must be specified.
  }
  \item{Y}{
    A data matrix with a column for the dependent variable must be specified.
  }
  \item{data}{
    The name of the used data-frame must be specified, so that the names of the variables can be extracted.
  }
 
}
\details{
  Given the matrix of the dependent variabe Y and the matrix of the independent variables X, the function ols returns an estimate 
  of the respective OLS-coefficients for each independent variable and the estimate of the intercept for the dependent variable.
  Additionally, ols also returns common summary statistics for the model.
}
\value{
  The output is a k+1x6 matrix with the name of the estimated coefficients, estimate, standard error, t-value and p-value; 
  and a summary of the residuals, residual standard error, degrees of freedom, multiple-R-squared and F-statistics
}

\author{
  Johannes Besch, \email{besch@ipz.uzh.ch}, Marco Radojevic \email{radojevic@ipz.uzh.ch}
}

\examples{
require(MASS)
model1 = glm(Sex ~ Height + Smoke + Pulse, data=MASS::survey, family=binomial(link=logit))
summary(model1)
}

\keyword{ models }