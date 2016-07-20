\name{corrZ2phi}
\alias{corrZ2phi}
\title{Computation of the Phi Coefficient from the Correlation of Bivariate Standard Normal Variables}
\description{This function computes the phi coefficient derived from dichotomizing bivariate standard normal variables.}
\usage{
corrZ2phi(corrZ, p1, p2)
}
\arguments{
  \item{corrZ}{The correlation of two standard normal variables.}
  \item{p1}{The expected value of the first variable after dichotomization.}
  \item{p2}{The expected value of the second variable after dichotomization.}
}
\value{The phi coefficient.}
\references{
Demirtas, H. (2016). A note on the relationship between the phi coefficient and the tetrachoric correlation under nonnormal underlying distributions. \emph{The American Statistician}, \bold{70(2)}, 143-148.
}
\seealso{\code{\link{tet2phi}}}
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

corrZ2phi(corrZ=-0.456, p1=0.85, p2=0.15)
}
