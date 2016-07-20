\name{mps2cps}
\alias{mps2cps}
\title{Computation of Cumulative Probabilities Given a Set of Marginal Probabilities}
\description{This function computes cumulative probabilities for each ordinal variable as defined by marginal probabilities provided in a list.}
\usage{
mps2cps(mps)
}
\arguments{
  \item{mps}{A list of marginal probability vectors corresponding to each ordinal variable. Each vector within the list \code{mps} must sum to 1.}
}
\value{A list of vectors containing cumulative probabilities for each set of marginal probabilities specified in \code{mps}. The i-th element of the list is a vector of the cumulative probabilities defining the marginal distribution of the i-th element of \code{mps}. If the i-th variable has k categories (k > 2), the i-th vector in the output will contain (k-1) probability values. The k-th element is implicitly 1.}
\examples{
mps2cps(list(c(0.4, 0.3, 0.2, 0.1), c(0.2, 0.2, 0.6)))
}
