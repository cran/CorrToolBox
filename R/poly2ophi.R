#find ordinal phi coefficient given polychoric correlation
poly2ophi <- function(polycorr, dist1, dist2) {
  if(min(dist1$p)<=0 | max(dist1$p)>=1) {
    stop("Elements of p for distribution 1 must be between 0 and 1.")
  }
  if(min(dist2$p)<=0 | max(dist2$p)>=1) {
    stop("Elements of p for distribution 2 must be between 0 and 1.")
  }
  if(sum(dist1$p)!=1) {
    stop('Marginal probabilities for distribution 1 must sum to 1.') 
  }
  if(sum(dist2$p)!=1) {
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
  
  #ensure that polychoric correlation is within feasible range
  corr.limits<-valid.limits.BinOrdNN(plist=NULL, skew.vec=c(dist1$skewness, dist2$skewness), kurto.vec=c(dist1$exkurtosis, dist2$exkurtosis), no.bin=0, no.ord=0, no.NN=2)
  
  if(polycorr<corr.limits$lower[2,1] | polycorr>corr.limits$upper[2,1]) {
    stop(paste('Specified polychoric correlation is not within the feasible correlation range of [',
               corr.limits$lower[2,1],', ', corr.limits$upper[2,1], '] for the given distributional characteristics.', sep=''))
  }

  
  #find polychoric correlation of Z
  polyZ<-corrY2corrZ(corrY=polycorr, 
                     skew.vec=c(dist1$skewness, dist2$skewness), 
                     kurto.vec=c(dist1$exkurtosis, dist2$exkurtosis))

  #solve for ordinal phi coefficient
  ophicoef<-corrZ2ophi(corrZ=polyZ, 
                       p1=dist1$p, 
                       p2=dist2$p)
  return(ophicoef)
}
