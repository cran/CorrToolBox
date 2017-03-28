\name{corrZ2ophi}
\alias{corrZ2ophi}
\title{Computation of the Ordinal Phi Coefficient from the Correlation of Bivariate Standard Normal Variables}
\description{This is an intermediate function that utilizes \code{\link{mps2cps}} to transform the specified marginal probabilities into cumulative probabilities and uses the \code{contord} function in the \code{GenOrd} package to compute the ordinal phi coefficient derived from discretizing bivariate standard normal variables.}
\usage{
corrZ2ophi(corrZ, p1, p2)
}
\arguments{
  \item{corrZ}{The correlation of two standard normal variables.}
  \item{p1}{A numeric vector containing marginal probabilities defining categories for the first ordinal variable.}
  \item{p2}{A numeric vector containing marginal probabilities defining categories for the second ordinal variable.}
}
\value{The ordinal phi coefficient.}
\references{
Demirtas, H., Ahmadian, R., Atis, S., Can, F.E., and Ercan, I. (2016). A nonnormal look at polychoric correlations: modeling the change in correlations before and after discretization. \emph{Computational Statistics}, \bold{31(4)}, 1385-1401.

Ferrari, P.A. and Barbiero, A. (2012). Simulating ordinal data. \emph{Multivariate Behavioral Research}, \bold{47(4)}, 566-589.
}
\seealso{\code{\link{mps2cps}}, \code{\link{poly2ophi}}}
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

corrZ2ophi(corrZ=0.502, p1=c(0.4, 0.3, 0.2, 0.1), p2=c(0.2, 0.2, 0.6))
}
