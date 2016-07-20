#find point-biserial correlation given biserial corelation
bs2pbs <- function (bs, bin.var, cont.var, p=NULL, cutpoint=NULL) {
  if(!is.null(cutpoint) & !is.null(p)) {
    stop("Must specify either p or cutpoint, not both.")
  }
  if(is.null(cutpoint) & is.null(p)) {
    stop("Must specify p or cutpoint.")
  }
  if(!is.null(cutpoint)) {
    y1.cdf<-ecdf(bin.var)
    q<-y1.cdf(cutpoint)
    p<-1-q
  }
  if(!is.null(p)) {
    if(p<=0 | p>=1) {
      stop("p must be between 0 and 1.")
    }
    q<-1-p
    cutpoint<-quantile(bin.var, prob=q)
  }

  
  #ensure that biserial correlation is within a feasible range
  y1.skew<-skewness(bin.var)
  y1.exkurt<-kurtosis(bin.var)-3  
  
  y2.skew<-skewness(cont.var)
  y2.exkurt<-kurtosis(cont.var)-3
  
  corr.limits<-valid.limits.BinOrdNN(plist=NULL, skew.vec=c(y1.skew, y2.skew), kurto.vec=c(y1.exkurt, y2.exkurt), no.bin=0, no.ord=0, no.NN=2)
  
  if(bs<corr.limits$lower[2,1] | bs>corr.limits$upper[2,1]) {
    stop(paste('Specified biserial correlation is not within the feasible correlation range of [',
               corr.limits$lower[2,1],', ', corr.limits$upper[2,1], '] for the given distributional characteristics.', sep=''))
  }
  
  #dichotomize Y1
  x1<-ifelse(test=bin.var>cutpoint, yes=1, no=0)
  
  #sample correlation between dichotomized y1 and y1
  x1.y1.c<-cor(x1, bin.var)
  
  #get point-biserial correlation
  pbs<-bs*x1.y1.c
  
  return(pbs)
}

