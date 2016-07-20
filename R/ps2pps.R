#find point-polyserial correlation given polyserial correlation
ps2pps <- function(ps, ord.var, cont.var, cats, p=NULL, cutpoint=NULL) {
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
    if(sum(p)!=1) {
      stop('Marginal probabilities must sum to 1.') 
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


  #ensure that polyserial correlation is within a feasible range
  y1.skew<-skewness(ord.var)
  y1.exkurt<-kurtosis(ord.var)-3  
  
  y2.skew<-skewness(cont.var)
  y2.exkurt<-kurtosis(cont.var)-3
  
  corr.limits<-valid.limits.BinOrdNN(plist=NULL, skew.vec=c(y1.skew, y2.skew), kurto.vec=c(y1.exkurt, y2.exkurt), no.bin=0, no.ord=0, no.NN=2)
  
  if(ps<corr.limits$lower[2,1] | ps>corr.limits$upper[2,1]) {
    stop(paste('Specified polyserial correlation is not within the feasible correlation range of [',
               corr.limits$lower[2,1],', ', corr.limits$upper[2,1], '] for the given distributional characteristics.', sep=''))
  }
  
  #discretize Y
  discY<-ordY(mp=p, cat=cats, y=ord.var)
  
  #sample correlation between discretized y1 and y1
  x1.y1.c<-cor(discY$x, discY$y)
  
  #get point-polyserial correlation
  pps<-ps*x1.y1.c
  
  return(pps)
}






