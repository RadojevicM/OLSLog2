\name{logit}
\alias{logit}

\title{
  Logistic Regression 
}
\description{
  The function estimates a logistic regression model through MLE. Method: Quasi-Newton 
}
\usage{
  logit(fml, data)
}

\arguments{
  \item{fml}{
    Formula-object of the disired independent and dependent variables. E.g. as.formular("sex ~ rural + age")
  }
  \item{data}{
    data frame which contains the variables specified in fml 
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



