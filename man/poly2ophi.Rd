\name{poly2ophi}
\alias{poly2ophi}
\title{Computation of the Ordinal Phi Coefficient from the Polychoric Correlation}
\description{This function computes the ordinal phi coefficient between two variables after both of the variables are ordinalized given the correlation before ordinalization (polychoric correlation) as seen in Demirtas et al. (2016). Before computation of the ordinal phi coefficient, the specified polychoric correlation is compared to the lower and upper correlation bounds of the two continuous variables as defined by the respective skewness and excess kurtosis using the generate, sort and correlate (GSC) algorithm in Demirtas and Hedeker (2011).}
\usage{
poly2ophi(polycorr, dist1, dist2)
}
\arguments{
  \item{polycorr}{The polychoric correlation.}
  \item{dist1}{A list of length 3 containing the skewness, excess kurtosis, and a numeric vector of marginal probabilities for the first continuous variable with names skewness, exkurtosis, and p, respectively.}
  \item{dist2}{A list of length 3 containing the skewness, excess kurtosis, and a numeric vector of marginal probabilities for the second continuous variable with names skewness, exkurtosis, and p, respectively.}
}
\value{The ordinal phi coefficient.}
\references{
Demirtas, H., Ahmadian, R., Atis, S., Can, F.E., and Ercan, I. (2016). A nonnormal look at polychoric correlations: modeling the change in correlations before and after discretization. \emph{Computational Statistics}, \bold{31(4)}, 1385-1401.

Demirtas, H. and Hedeker, D. (2011). A practical way for computing approximate lower and upper correlation bounds. \emph{The American Statistician}, \bold{65(2)}, 104-109.

Ferrari, P.A. and Barbiero, A. (2012). Simulating ordinal data. \emph{Multivariate Behavioral Research}, \bold{47(4)}, 566-589.
}
\seealso{\code{\link{corrY2corrZ}}, \code{\link{corrZ2ophi}}, \code{\link{mps2cps}}}
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

poly2ophi(polycorr=0.5, 
          dist1=list(skewness=y1.skew, exkurtosis=y1.exkurt, p=c(0.4, 0.3, 0.2, 0.1)),
          dist2=list(skewness=y2.skew, exkurtosis=y2.exkurt , p=c(0.2, 0.2, 0.6)))

poly2ophi(polycorr=0.5, 
          dist1=list(skewness=y1.skew, exkurtosis=y1.exkurt, p=c(0.1, 0.1, 0.1, 0.7)),
          dist2=list(skewness=y2.skew, exkurtosis=y2.exkurt , p=c(0.8, 0.1, 0.1)))
}
