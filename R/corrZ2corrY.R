#function for finding tetrachoric/polychoric correlation of Y from tetrachoric/polychoric correlation of Z
corrZ2corrY<-function(corrZ, skew.vec, kurto.vec) {
  #ensure skewness and excess kurtosis values are within possible range
  check.values<-try(validation.skewness.kurtosis(n.NN=2,
                                                 skewness.vec = skew.vec, 
                                                 kurtosis.vec = kurto.vec))
  if(check.values!=TRUE) {
    stop
  }
  
  #ensure that corrZ is within feasible range
  corr.limits<-valid.limits.BinOrdNN(plist=NULL, skew.vec=c(0,0), kurto.vec=c(0,0), no.bin=0, no.ord=0, no.NN=2)
  
  if(corrZ<corr.limits$lower[2,1] | corrZ>corr.limits$upper[2,1]) {
    stop(paste('Specified correlation is not within the feasible correlation range of [',
               corr.limits$lower[2,1],', ', corr.limits$upper[2,1], '] for bivariate standard normal variables.', sep=''))
  }
  
  #solve for fleishman polynomials
  coef.mat<-Fleishman.coef.NN(skew.vec=skew.vec, 
                              kurto.vec=kurto.vec)
  
  F1<-coef.mat[1,]
  F2<-coef.mat[2,]
  
  #solve for tetrachoric correlation of Y1, Y2
  B<-(F1['b']*F2['b'])+3*(F1['b']*F2['d'])+3*(F1['d']*F2['b'])+9*(F1['d']*F2['d'])
  C<-2*(F1['c']*F2['c'])
  D<-6*(F1['d']*F2['d'])
  
  corrY<-(corrZ)*B+((corrZ)^2)*C+((corrZ)^3)*D
  names(corrY)<-NULL
  
  return(corrY)
}


