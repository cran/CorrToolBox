\name{CorrToolBox-package}
\alias{CorrToolBox-package}
\alias{CorrToolBox}
\docType{package}
\title{
Modeling Correlational Magnitude Transformations in Discretization Contexts
}
\description{
This package implements the computational algorithms for modeling the correlation transitions under specified distributional assumptions within the realm of discretization in the context of the latency and threshold concepts. Functions that compute the correlational magnitude changes in both directions (identification of the pre-discretization correlation value in order to attain a specified post-discretization magnitude, and the other way around) are provided.

This package consists of eight main functions. Computing the tetrachoric correlation from the phi coefficient and vice versa are done in \code{\link{phi2tet}} and \code{\link{tet2phi}}, respectively. Computing the polychoric correlation from the ordinal phi coefficient and vice versa are done in \code{\link{ophi2poly}} and \code{\link{poly2ophi}}, respectively. Computing the biserial correlation from the point-biserial correlation and vice versa are done in \code{\link{pbs2bs}} and \code{\link{bs2pbs}}, respectively. Computing the polyserial correlation from the point-polyserial correlation and vice versa are done in \code{\link{pps2ps}} and \code{\link{ps2pps}}, respectively. 

Auxiliary functions are also provided. \code{\link{corrY2corrZ}}, \code{\link{corrZ2corrY}}, \code{\link{corrZ2ophi}}, \code{\link{corrZ2phi}}, and \code{\link{ophi2corrZ}} are intermediate functions utilized within the main functions but can be used as stand-alone functions. \code{\link{ordY}} discretizes a continuous variable, and \code{\link{mps2cps}} provides cumulative probabilities for each set of marginal probabilities in a list. Additional intermediate functions from imported packages include \code{phi2tetra} from the \code{psych} package, \code{ordcont} and \code{contord} from the \code{GenOrd} package, and \code{pmvnorm} from the \code{mvtnorm} package. 

Within each correlation transition function, the correlation boundaries for the given marginal distributions are compared to the specified input correlation to ensure there are no violations according to Demirtas and Hedeker (2011). The function \code{valid.limits.BinOrdNN} in the package \code{BinOrdNonNor} is utilized for this step. Additionally, \code{Fleishman.coef.NN} in the package \code{BinOrdNonNor} is used wherever Fleishman coefficients need to be calculated for a continuous variable.
}
\details{
\tabular{ll}{
Package: \tab CorrToolBox\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-07-19\cr
License: \tab GPL-2 | GPL-3\cr
}}
\author{
Rawan Allozi, Hakan Demirtas

Maintainer: Rawan Allozi <ralloz2@uic.edu>
}
\references{
Demirtas, H. and Hedeker, D. (2011). A practical way for computing approximate lower and upper correlation bounds. \emph{The American Statistician}, \bold{65(2)}, 104-109.

Demirtas, H., Hedeker, D., and Mermelstein, R. J. (2012). Simulation of massive public health data by power polynomials. \emph{Statistics in Medicine}, \bold{31(27)}, 3337-3346.

Demirtas, H. (2016). A note on the relationship between the phi coefficient and the tetrachoric correlation under nonnormal underlying distributions. \emph{The American Statistician}, \bold{70(2)}, 143-148.

Demirtas, H. and Hedeker, D. (2016). Computing the point-biserial correlation under any underlying continuous distribution. \emph{Communications in Statistics-Simulation and Computation}, \bold{45(8)}, 2744-2751.

Demirtas, H., Ahmadian, R., Atis, S., Can, F.E., and Ercan, I. (2016). A nonnormal look at polychoric correlations: modeling the change in correlations before and after discretization. Forthcoming in \emph{Computational Statistics}.

Ferrari, P.A. and Barbiero, A. (2012). Simulating ordinal data. \emph{Multivariate Behavioral Research}, \bold{47(4)}, 566-589.

Fleishman A.I. (1978). A method for simulating non-normal distributions. \emph{Psychometrika}, \bold{43(4)}, 521-532.

Vale, C.D. and Maurelli, V.A. (1983). Simulating multivariate nonnormal distributions. 
\emph{Psychometrika}, \bold{48(3)}, 465-471.
}
