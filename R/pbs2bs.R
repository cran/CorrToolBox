#find biserial corelation given point-biserial correlation
pbs2bs <- function (pbs, bin.var, cont.var, p=NULL, cutpoint=NULL) {
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
    q<-1-p
    cutpoint<-quantile(bin.var, prob=q)
  }
  if(p<=0 | p>=1) {
    stop("p must be between 0 and 1.")
  }
  
  #ensure that pbs is within feasible range
  y2.skew<-skewness(cont.var)
  y2.exkurt<-kurtosis(cont.var)-3
    
  corr.limits<-valid.limits.BinOrdNN(plist=list(p), skew.vec=y2.skew, kurto.vec=y2.exkurt, no.bin=1, no.ord=0, no.NN=1)
  
  if(pbs<corr.limits$lower[2,1] | pbs>corr.limits$upper[2,1]) {
    stop(paste('Specified point-biserial correlation is not within the feasible correlation range of [',
               corr.limits$lower[2,1],', ', corr.limits$upper[2,1], '] for the given distributional characteristics.', sep=''))
  }  

  x1<-ifelse(test=bin.var>cutpoint, yes=1, no=0)

  #sample correlation between dichotomized y1 and y1
  x1.y1.c<-cor(x1, bin.var)

  #get biserial correlation
  bs<-pbs/x1.y1.c
  
  return(bs)
}



