#find tetrachoric correlation given phi coefficent
phi2tet <- function(phicoef, dist1, dist2) {
  if(dist1$p<=0 | dist1$p>=1) {
    stop("p for distribution 1 must be between 0 and 1.")
  }
  if(dist2$p<=0 | dist2$p>=1) {
    stop("p for distribution 2 must be between 0 and 1.")
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

  #ensure that phi coefficient is within feasible range
  corr.limits<-valid.limits.BinOrdNN(plist=list(dist1$p, dist2$p), no.bin=2, no.ord=0, no.NN=0)
  
  if(phicoef<corr.limits$lower[2,1] | phicoef>corr.limits$upper[2,1]) {
    stop(paste('Specified phi coefficient is not within the feasible correlation range of [',
               corr.limits$lower[2,1],', ', corr.limits$upper[2,1], '] for the given distributional characteristics.', sep=''))
  }  
  
  #find tetrachoric for Z1, Z2 given phi coefficient
  tetZ<-phi2tetra(ph=phicoef, m=c(dist1$p, dist2$p))
  
  #find tetrachoric correlation of Y from tetrachoric correlation of Z
  tetY<-corrZ2corrY(corrZ=tetZ, 
                    skew.vec=c(dist1$skewness, dist2$skewness), 
                    kurto.vec=c(dist1$exkurtosis, dist2$exkurtosis))
  
  return(tetY)
}

