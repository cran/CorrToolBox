#function for finding tetrachoric/polychoric correlation of Z from tetrachoric/polychoric correlation of Y
corrY2corrZ<-function(corrY, skew.vec, kurto.vec) {
  #ensure skewness and excess kurtosis values are within possible range
  check.values<-try(validation.skewness.kurtosis(n.NN=2,
                                                 skewness.vec = skew.vec, 
                                                 kurtosis.vec = kurto.vec))
  if(check.values!=TRUE) {
    stop
  }
  
  #ensure that corrY is within feasible range
  corr.limits<-valid.limits.BinOrdNN(plist=NULL, skew.vec=skew.vec, kurto.vec=kurto.vec, no.bin=0, no.ord=0, no.NN=2)
  
  if(corrY<corr.limits$lower[2,1] | corrY>corr.limits$upper[2,1]) {
    stop(paste('Specified correlation is not within the feasible correlation range of [',
               corr.limits$lower[2,1],', ', corr.limits$upper[2,1], '] for the given distributional characteristics.', sep=''))
  }
  
  #solve for fleishman polynomials
  coef.mat<-Fleishman.coef.NN(skew.vec=skew.vec, 
                              kurto.vec=kurto.vec)
  
  #solve for tetrachoric/polychoric correlation of Z
  F1<-coef.mat[1,]
  F2<-coef.mat[2,]
  
  A<-corrY
  B<-(F1['b']*F2['b'])+3*(F1['b']*F2['d'])+3*(F1['d']*F2['b'])+9*(F1['d']*F2['d'])
  C<-2*(F1['c']*F2['c'])
  D<-6*(F1['d']*F2['d'])
  
  corrZ<-polyroot(c(-A, B, C, D))[1]
  corrZ<-Re(corrZ)[abs(Im(corrZ)) < 1e-6]
  
  return(corrZ)
}
