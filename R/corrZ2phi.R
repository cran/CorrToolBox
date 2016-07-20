#function for finding phi coefficient from tetrachoric/polychoric correlation of Z 
corrZ2phi<-function(corrZ, p1, p2) {
  if(p1<=0 | p1>=1) {
    stop("p for distribution 1 must be between 0 and 1.")
  }
  if(p2<=0 | p2>=1) {
    stop("p for distribution 2 must be between 0 and 1.")
  }
  
  #ensure that corrZ is within feasible range
  corr.limits<-valid.limits.BinOrdNN(plist=NULL, skew.vec=c(0,0), kurto.vec=c(0,0), no.bin=0, no.ord=0, no.NN=2)
  
  if(corrZ<corr.limits$lower[2,1] | corrZ>corr.limits$upper[2,1]) {
    stop(paste('Specified correlation is not within the feasible correlation range of [',
               corr.limits$lower[2,1],', ', corr.limits$upper[2,1], '] for bivariate standard normal variables.', sep=''))
  }
  
  #solve for phi coefficient
  zp1<-qnorm(p=p1, mean=0, sd=1)
  zp2<-qnorm(p=p2, mean=0, sd=1)
  
  corrmat <- diag(2)
  corrmat[lower.tri(corrmat)] <- corrZ
  corrmat[upper.tri(corrmat)] <- corrZ
  
  meanvec<-rep(0,2)
  
  cdfprob <- as.numeric(pmvnorm(lower=-Inf, upper=c(zp1,zp2), mean=meanvec, corr=corrmat))
  phicoef<-(cdfprob-p1*p2)/(sqrt((p1)*(1-p1)*(p2)*(1-p2)))
  
  return(phicoef)
}

