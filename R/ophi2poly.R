#find polychoric correlation given ordinal phi coefficient
ophi2poly <- function (ophicoef, dist1, dist2) {
  if(min(dist1$p)<=0 | max(dist1$p)>=1) {
    stop("Elements of p for distribution 1 must be between 0 and 1.")
  }
  if(min(dist2$p)<=0 | max(dist2$p)>=1) {
    stop("Elements of p for distribution 2 must be between 0 and 1.")
  }
  if(sum(dist1$p)>(1+.Machine$double.eps^0.5) | sum(dist1$p)<(1-.Machine$double.eps^0.5)) { #tolerance added for use across platforms   
    stop('Marginal probabilities for distribution 1 must sum to 1.') 
  }
  if(sum(dist2$p)>(1+.Machine$double.eps^0.5) | sum(dist2$p)<(1-.Machine$double.eps^0.5)) {
    stop('Marginal probabilities for distribution 2 must sum to 1.') 
  }
  if(is.null(dist1$skewness) | is.null(dist1$exkurtosis) | is.null(dist1$p)) {
    stop("Skewness, excess kurtosis, and p for distribution 1 must be specified.")
  }
  if(is.null(dist2$skewness) | is.null(dist2$exkurtosis) | is.null(dist2$p)) {
    stop("Skewness, excess kurtosis, and p for distribution 2 must be specified.")
  }
  
  #ensure skewness and excess kurtosis values are within possible range
  check.values<-try(validation.skewness.kurtosis(n.NN=2,
                                                 skewness.vec = c(dist1$skewness, dist2$skewness), 
                                                 kurtosis.vec = c(dist1$exkurtosis, dist2$exkurtosis)))
  if(check.values!=TRUE) {
    stop
  }

  #ensure that ordinal phi coefficient is within feasible range
  #create cumulative probabilities
  cps<-mps2cps(mps=list(dist1$p, dist2$p))
  corr.limits<-valid.limits.BinOrdNN(plist=cps, no.bin=0, no.ord=2, no.NN=0)
  
  if(ophicoef<corr.limits$lower[2,1] | ophicoef>corr.limits$upper[2,1]) {
    stop(paste('Specified ordinal phi coefficient is not within the feasible correlation range of [',
               corr.limits$lower[2,1],', ', corr.limits$upper[2,1], '] for the given distributional characteristics.', sep=''))
  }  
  
  #find polychoric correlation of Z from ordinal phi coefficient
  polyZ<-ophi2corrZ(ophi=ophicoef, p1=dist1$p, p2=dist2$p)
  
  #find polychoric correlation of Y from polychoric correlation of Z
  polyY<-corrZ2corrY(corrZ=polyZ, 
                     skew.vec=c(dist1$skewness, dist2$skewness), 
                     kurto.vec=c(dist1$exkurtosis, dist2$exkurtosis))
  return(polyY)
}

