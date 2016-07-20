#function for ordinalizing continuous variable
ordY<-function(mp, cat, y) {
  if(min(mp)<=0 | max(mp)>=1) {
    stop("Marginal probabilities must be between 0 and 1.")
  }   
  if(sum(mp)!=1) {
    stop('Marginal probabilities must sum to 1.') 
  }
  if(length(mp)!=length(cat)) {
    stop('There must be a corresponding probability for each given category.')
  }
  #get cumulative probabilities
  cp<-mps2cps(mps=list(mp))
  
  #ordinalize y1
  ocats<-data.frame(n=cat, pmin=c(0,cp[[1]]), pmax=c(cp[[1]],1))
  ocats<-ocats[which(ocats$pmin!=ocats$pmax),]
  
  ocats$min<-quantile(y, ocats$pmin)
  ocats$max<-quantile(y, ocats$pmax)
  
  y.df<-data.frame(y=y, x=rep(NA, length(y)))
  
  for(i in ocats$n) {
    min<-ocats[which(ocats$n==i),'min']
    max<-ocats[which(ocats$n==i),'max']
    if(i==min(ocats$n)) {
      y.df[which(y.df[,'y']==min), 'x']<-i     
    }
    y.df[which(y.df$y>min & y.df$y<=max), 'x']<-i
  }
  
  return(y.df)
}

