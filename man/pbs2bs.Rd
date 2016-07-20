\name{pbs2bs}
\alias{pbs2bs}
\title{Computation of the Biserial Correlation from the Point-Biserial Correlation}
\description{This function computes the biserial correlation between two continuous variables given the correlation after dichotomization of one of the variables (point-biserial correlation) as seen in Demirtas and Hedeker (2016). Before computation of the biserial correlation, the specified point-biserial correlation is compared to the lower and upper correlation bounds of the continuous variable and binary variable using the generate, sort and correlate (GSC) algorithm in Demirtas and Hedeker (2011).}
\usage{
pbs2bs(pbs, bin.var, cont.var, p=NULL, cutpoint=NULL)
}
\arguments{
  \item{pbs}{The point-biserial correlation.}
  \item{bin.var}{A numeric vector of the continuous variable before dichotomization.}
  \item{cont.var}{A numeric vector of the the continuous variable that is not transformed.}
  \item{p}{The expected value of the numeric vector \code{bin.var} after dichotomization. Either \code{p} or \code{cutpoint} should be specified.}
  \item{cutpoint}{The value at which the vector \code{bin.var} should be dichotomized. Either \code{p} or \code{cutpoint} should be specified.}
}
\value{The biserial correlation.}
\references{
Demirtas, H. and Hedeker, D. (2011). A practical way for computing approximate lower and upper correlation bounds. \emph{The American Statistician}, \bold{65(2)}, 104-109.

Demirtas, H. and Hedeker, D. (2016). Computing the point-biserial correlation under any underlying continuous distribution. \emph{Communications in Statistics-Simulation and Computation}, \bold{45(8)}, 2744-2751.}
\examples{
set.seed(123)
y1<-rweibull(n=100000, scale=1, shape=1.2)

gaussmix <- function(n,m1,m2,s1,s2,pi) {
  I <- runif(n)<pi
  rnorm(n,mean=ifelse(I,m1,m2),sd=ifelse(I,s1,s2))
}
y2<-gaussmix(n=100000, m1=0, s1=1, m2=3, s2=1, pi=0.6)

pbs2bs(pbs=0.25, bin.var=y1, cont.var=y2, p=0.55)
pbs2bs(pbs=0.25, bin.var=y1, cont.var=y2, cutpoint=0.65484)
}
