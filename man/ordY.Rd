\name{ordY}
\alias{ordY}
\title{Ordinalization of a Continuous Variable}
\description{This functions creates an ordinalized form of a continuous variable.}
\usage{
ordY(mp, cat, y)
}
\arguments{
  \item{mp}{A vector of marginal probabilities defining the ordinalized variable.}
  \item{cat}{A numeric vector containing the categories for each respective marginal probability in \code{mp}.}
  \item{y}{A continuous variable to be ordinalized into categories in \code{cat} as defined by \code{mp}.} 
}
\value{A data frame containing the given continuous variable and the ordinalized variable with names y and x, respectively.}
\seealso{\code{\link{mps2cps}}}
\examples{
y<-rnorm(100000)
dat<-ordY(mp=c(0.25, 0.5, 0.25), cat=c(1,2,3), y=y)
}
