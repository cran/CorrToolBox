\name{ps2pps}
\alias{ps2pps}
\title{Computation of the Point-Polyserial Correlation from the Polyserial Correlation}
\description{This function computes the point-polyserial correlation between two variables after one of the variables is ordinalized given the correlation before ordinalization (polyserial correlation) as seen in Demirtas and Hedeker (2016). Before computation of the point-polyserial correlation, the specified polyserial correlation is compared to the lower and upper correlation bounds of the two continuous variables using the generate, sort and correlate (GSC) algorithm in Demirtas and Hedeker (2011).}
\usage{
ps2pps(ps, ord.var, cont.var, cats, p=NULL, cutpoint=NULL)
}
\arguments{
  \item{ps}{The polyserial correlation.}
  \item{ord.var}{A numeric vector of the continuous variable before ordinalization.}
  \item{cont.var}{A numeric vector of the the continuous variable that is not transformed.}
  \item{cats}{A numeric vector of the categories in the ordinalization of \code{ord.var}.}
  \item{p}{A numeric vector of the marginal probabilities corresponding to each category in \code{cats}. The marginal probabilities must sum to 1. Either \code{p} or \code{cutpoint} should be specified.}
  \item{cutpoint}{A numeric vector of the cutpoints used to define the categories in \code{cats}. Either \code{p} or \code{cutpoint} should be specified.}
}
\value{The point-polyserial correlation.}
\references{
Demirtas, H. and Hedeker, D. (2011). A practical way for computing approximate lower and upper correlation bounds. \emph{The American Statistician}, \bold{65(2)}, 104-109.

Demirtas, H. and Hedeker, D. (2016). Computing the point-biserial correlation under any underlying continuous distribution. \emph{Communications in Statistics-Simulation and Computation}, \bold{45(8)}, 2744-2751.
}
\seealso{\code{\link{ordY}}, \code{\link{mps2cps}}}
\examples{
set.seed(234)
y1<-rweibull(n=100000, scale=1, shape=25)

gaussmix <- function(n,m1,m2,s1,s2,pi) {
  I <- runif(n)<pi
  rnorm(n,mean=ifelse(I,m1,m2),sd=ifelse(I,s1,s2))
}
y2<-gaussmix(n=100000, m1=0, s1=1, m2=2, s2=1, pi=0.5)

ps2pps(ps=0.6, ord.var=y1, cont.var=y2, cats=c(1,2,3,4), p=c(0.4, 0.3, 0.2, 0.1))
ps2pps(ps=0.6, ord.var=y1, cont.var=y2, cats=c(1,2,3,4), cutpoint=c(0.97341, 1.00750, 1.03421))
}
