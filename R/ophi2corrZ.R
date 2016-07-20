ophi2corrZ<-function(ophi, p1, p2) {
  if(min(p1)<=0 | max(p1)>=1) {
    stop("Elements of p for distribution 1 must be between 0 and 1.")
  }
  if(min(p2)<=0 | max(p2)>=1) {
    stop("Elements of p for distribution 2 must be between 0 and 1.")
  }
  if(sum(p1)!=1) {
    stop('Marginal probabilities for distribution 1 must sum to 1.') 
  }
  if(sum(p2)!=1) {
    stop('Marginal probabilities for distribution 2 must sum to 1.') 
  }
  
  #create cumulative probabilities
  cps<-mps2cps(mps=list(p1, p2))
  
  #ensure that ordinal phi coefficient is within feasible range
  corr.limits<-valid.limits.BinOrdNN(plist=cps, no.bin=0, no.ord=2, no.NN=0)
  
  if(ophi<corr.limits$lower[2,1] | ophi>corr.limits$upper[2,1]) {
    stop(paste('Specified ordinal phi coefficient is not within the feasible correlation range of [',
               corr.limits$lower[2,1],', ', corr.limits$upper[2,1], '] for the given distributional characteristics.', sep=''))
  }  
  
  #find polychoric for Z1, Z2 given ordinal phi coefficient
  corrmat.OO<- diag(2)
  corrmat.OO[lower.tri(corrmat.OO)] <- ophi
  corrmat.OO[upper.tri(corrmat.OO)] <- ophi
  
  corrZ<-ordcont(marginal=cps, Sigma=corrmat.OO)[[1]][2,1]
  return(corrZ)
}

