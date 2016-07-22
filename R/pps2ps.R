#find polyserial corelation given point-polyserial correlation
pps2ps <- function (pps, ord.var, cont.var, cats, p=NULL, cutpoint=NULL) {
  if(!is.null(cutpoint) & !is.null(p)) {
    stop("Must specify either p or cutpoint, not both.")
  }
  if(is.null(cutpoint) & is.null(p)) {
    stop("Must specify p or cutpoint.")
  }
  if(!is.null(p)) {
    if(min(p)<=0 | max(p)>=1) {
      stop("Elements of p for distribution 1 must be between 0 and 1.")
    }   
    if(sum(p)>(1+.Machine$double.eps^0.5) | sum(p)<(1-.Machine$double.eps^0.5)) { #tolerance added for use across platforms
      stop('Marginal probabilities must sum to 1.') 
    }
    if(length(p)!=length(cats)) {
      stop('There must be a corresponding probability for each category.')
    }
    cps<-mps2cps(mps=list(p))[[1]]
  }
  if(!is.null(cutpoint)) {
    if((length(cutpoint)+1)!=length(cats)) {
      stop('The number of categories must be equal to the number of cutpoints + 1.')
    }
    cutpoint<-sort(cutpoint)
    cats<-sort(cats)
    y1.cdf<-ecdf(ord.var)
    if(min(cutpoint)<min(ord.var)) {
      min.remove<-which(cutpoint<min(ord.var))
      n.minr<-length(min.remove)
      cutpoint<-cutpoint[-min.remove]
      cats<-cats[(n.minr+1):length(cats)]
    }
    if(max(cutpoint)>max(ord.var)) {
      max.remove<-which(cutpoint>max(ord.var))
      n.maxr<-length(max.remove)
      cutpoint<-cutpoint[-max.remove]
      cats<-cats[1:(length(cats)-n.maxr)]
    }
    cps<-y1.cdf(cutpoint)
    p<-c(cps, 1)-c(0,cps)
  }
  
  
  #ensure that pps is within feasible range
  y2.skew<-skewness(cont.var)
  y2.exkurt<-kurtosis(cont.var)-3
  
  corr.limits<-valid.limits.BinOrdNN(plist=list(cps), skew.vec=y2.skew, kurto.vec=y2.exkurt, no.bin=0, no.ord=1, no.NN=1)
  
  if(pps<corr.limits$lower[2,1] | pps>corr.limits$upper[2,1]) {
    stop(paste('Specified point-polyserial correlation is not within the feasible correlation range of [',
               corr.limits$lower[2,1],', ', corr.limits$upper[2,1], '] for the given distributional characteristics.', sep=''))
  }  
  
  #discretize Y
  discY<-ordY(mp=p, cat=cats, y=ord.var)
  
  #sample correlation between discretized y1 and continuous y1
  x1.y1.c<-cor(discY$x, discY$y)
  
  #get polyserial correlation
  ps<-pps/x1.y1.c
  
  return(ps)
}



