\name{ophi2corrZ}
\alias{ophi2corrZ}
\title{Computation of the Correlation of Bivariate Standard Normal Variables from the Ordinal Phi Coefficient}
\description{This is an intermediate function that transforms marginal probabilities into cumulative probabilities and uses the \code{ordcont} function in the \code{GenOrd} package to compute the correlation of bivariate standard normal variables from the ordinal phi coefficient.}
\usage{
ophi2corrZ(ophi, p1, p2)
}
\arguments{
  \item{ophi}{The ordinal phi coefficient.}
  \item{p1}{A numeric vector containing marginal probabilities defining categories for the first ordinal variable.}
  \item{p2}{A numeric vector containing marginal probabilities defining categories for the second ordinal variable.}
}
\value{The correlation of standard normal variables.}
\references{
Demirtas, H., Ahmadian, R., Atis, S., Can, F.E., and Ercan, I. (2016). A nonnormal look at polychoric correlations: modeling the change in correlations before and after discretization. Forthcoming in \emph{Computational Statistics}.

Ferrari, P.A. and Barbiero, A. (2012). Simulating ordinal data. \emph{Multivariate Behavioral Research}, \bold{47(4)}, 566-589.
}
\seealso{\code{\link{mps2cps}}, \code{\link{ophi2poly}}}
\examples{
set.seed(567)
library(moments)

y1<-rweibull(n=100000, scale=1, shape=3.6)
y1.skew<-round(skewness(y1), 5)
y1.exkurt<-round(kurtosis(y1)-3, 5)

gaussmix <- function(n,m1,m2,s1,s2,pi) {
  I <- runif(n)<pi
  rnorm(n,mean=ifelse(I,m1,m2),sd=ifelse(I,s1,s2))
}
y2<-gaussmix(n=100000, m1=0, s1=1, m2=2, s2=1, pi=0.3)
y2.skew<-round(skewness(y2), 5)
y2.exkurt<-round(kurtosis(y2)-3, 5)

ophi2corrZ(ophi=-0.7, p1=c(0.4, 0.3, 0.2, 0.1), p2=c(0.2, 0.2, 0.6))
}
