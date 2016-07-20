\name{corrY2corrZ}
\alias{corrY2corrZ}
\title{Computation of the Correlation of Bivariate Standard Normal Variables from the Correlation of Bivariate Nonnormal Variables}
\description{This is an intermediate function that computes the correlation of bivariate standard normal variables from the correlation of continuous nonnormal variables. Fleishman coefficients for each nonnormal variable with the specified skewness and excess kurtosis are found. The Fleishman coefficients and correlation of nonnormal variables are used to find the correlation of the two respective standard normal variables as seen in Demirtas, Hedeker, and Mermelstein (2012).}
\usage{
corrY2corrZ(corrY, skew.vec, kurto.vec)
}
\arguments{
  \item{corrY}{The correlation of two continuous nonnormal variables.}
  \item{skew.vec}{The skewness vector for continuous variables.}
  \item{kurto.vec}{The kurtosis vector for continuous variables.}
}
\value{The correlation of the two respective standard normal variables.}
\references{
Demirtas, H., Hedeker, D., and Mermelstein, R. J. (2012). Simulation of massive public health data by power polynomials. \emph{Statistics in Medicine}, \bold{31(27)}, 3337-3346.

Fleishman A.I. (1978). A method for simulating non-normal distributions. \emph{Psychometrika}, \bold{43(4)}, 521-532.
}
\seealso{\code{\link{tet2phi}}, \code{\link{poly2ophi}}}
\examples{
set.seed(987)
library(moments)

y1<-rweibull(n=100000, scale=1, shape=1)
y1.skew<-round(skewness(y1), 5)
y1.exkurt<-round(kurtosis(y1)-3, 5)

gaussmix <- function(n,m1,m2,s1,s2,pi) {
  I <- runif(n)<pi
  rnorm(n,mean=ifelse(I,m1,m2),sd=ifelse(I,s1,s2))
}
y2<-gaussmix(n=100000, m1=0, s1=1, m2=3, s2=1, pi=0.5)
y2.skew<-round(skewness(y2), 5)
y2.exkurt<-round(kurtosis(y2)-3, 5)

corrY2corrZ(corrY=-0.4, skew.vec=c(y1.skew, y2.skew), kurto.vec=c(y1.exkurt, y2.exkurt))
}
