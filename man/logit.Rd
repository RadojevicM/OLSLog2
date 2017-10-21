\name{logit}
\alias{logit}

\title{
  Logistic Regression 
}
\description{
  The function estimates a logistic regression model through MLE. Method: Quasi-Newton 
}
\usage{
  logit(fml, dat)
}

\arguments{
  \item{fml}{
    Formula-object of the disired independent and dependent variables. E.g. as.formular("sex ~ rural + age")
  }
  \item{dat}{
    data frame which contains the variables specified in fml 
  }
 
}
\details{
 Estimates a logistic regression
}
\value{
  Output is logistic regression
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



