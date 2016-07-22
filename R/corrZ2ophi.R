#function for finding ordinal phi coefficient from tetrachoric/polychoric correlation of Z
corrZ2ophi<-function(corrZ, p1, p2) {
  if(min(p1)<=0 | max(p1)>=1) {
    stop("Elements of p for distribution 1 must be between 0 and 1.")
  }
  if(min(p2)<=0 | max(p2)>=1) {
    stop("Elements of p for distribution 2 must be between 0 and 1.")
  }
  if(sum(p1)>(1+.Machine$double.eps^0.5) | sum(p1)<(1-.Machine$double.eps^0.5)) { #tolerance added for use across platforms
    stop('Marginal probabilities for distribution 1 must sum to 1.') 
  }
  if(sum(p2)>(1+.Machine$double.eps^0.5) | sum(p2)<(1-.Machine$double.eps^0.5)) {
    stop('Marginal probabilities for distribution 2 must sum to 1.') 
  }
  
  #create cumulative probabilities
  cps<-mps2cps(mps=list(p1, p2))
  
  #ensure that corrZ is within feasible range
  corr.limits<-valid.limits.BinOrdNN(plist=NULL, skew.vec=c(0,0), kurto.vec=c(0,0), no.bin=0, no.ord=0, no.NN=2)
  
  if(corrZ<corr.limits$lower[2,1] | corrZ>corr.limits$upper[2,1]) {
    stop(paste('Specified correlation is not within the feasible correlation range of [',
               corr.limits$lower[2,1],', ', corr.limits$upper[2,1], '] for bivariate standard normal variables.', sep=''))
  }
  
  #find ordinal phi coefficient
  corrmat.NN<- diag(2)
  corrmat.NN[lower.tri(corrmat.NN)] <- corrZ
  corrmat.NN[upper.tri(corrmat.NN)] <- corrZ
  
  ophicoef<-contord(marginal=cps, Sigma=corrmat.NN)[2,1]
  return(ophicoef)
}






