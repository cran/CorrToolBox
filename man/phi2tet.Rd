\name{phi2tet}
\alias{phi2tet}
\title{Computation of the Tetrachoric Correlation from the Phi Coefficient}
\description{This function computes the tetrachoric correlation between two continuous variables given the correlation after dichotomization of both variables (phi coefficient) as seen in Demirtas (2016). Before computation of the tetrachoric correlation, the specified phi coefficient is compared to the lower and upper correlation bounds for the two binary variables using the generate, sort and correlate (GSC) algorithm in Demirtas and Hedeker (2011).}
\usage{
phi2tet(phicoef, dist1, dist2)
}
\arguments{
  \item{phicoef}{The phi coefficient.}
  \item{dist1}{A list of length 3 containing the skewness, excess kurtosis, and expected value after dichotomization for the first continuous variable with names skewness, exkurtosis, and p, respectively.}
  \item{dist2}{A list of length 3 containing the skewness, excess kurtosis, and expected value after dichotomization for the second continuous variable with names skewness, exkurtosis, and p, respectively.}
}
\value{The tetrachoric correlation.}
\references{
Demirtas, H. (2016). A note on the relationship between the phi coefficient and the tetrachoric correlation under nonnormal underlying distributions. \emph{The American Statistician}, \bold{70(2)}, 143-148.

Demirtas, H. and Hedeker, D. (2011). A practical way for computing approximate lower and upper correlation bounds. \emph{The American Statistician}, \bold{65(2)}, 104-109.
}
\seealso{\code{\link{corrZ2corrY}}}
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

phi2tet(phicoef=0.1, 
        dist1=list(skewness=y1.skew, exkurtosis=y1.exkurt, p=0.85), 
        dist2=list(skewness=y2.skew, exkurtosis=y2.exkurt, p=0.15))

phi2tet(phicoef=0.5, 
        dist1=list(skewness=y1.skew, exkurtosis=y1.exkurt, p=0.10), 
        dist2=list(skewness=y2.skew, exkurtosis=y2.exkurt, p=0.30))
}
