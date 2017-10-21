\name{logit}
\alias{logit}

\title{
  Logistic Regression 
}
\description{
  The function estimates a logistic regression model with the Maximum-Likelihood Method. 
}
\usage{
  logit(fml, dat)
}

\arguments{
  \item{fml}{
    Formula-object of the desired independent and dependent variables: formula.object <- as.formula("y ~ x1 + x2 + ... xn").
  }
  \item{dat}{
    Data frame which contains the variables specified in the formula-object. 
  }
 
}
\details{
 Estimates a logistic regression for a binary independent variable where the data generating process can be described with a Bernoulli-distribution. Parameters are maximised with the Maximum-Likelihood-Method using a Quasi-Newton-Algorithm.
}
\value{
  Output are the coefficient estimates of a logistic regression and its respective standard errors.
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



