#function that returns cumulative probabilities for each vector of marginal probabilities
mps2cps<-function(mps) {
  if(is.list(mps)==FALSE) {
    stop('mps must be a list.')
  }
  if(all(unlist(lapply(mps, sum))==1)!=TRUE) {
    stop('Marginal probabilities in each respective element of the given list must sum to 1.')
  }
  if(min(unlist(lapply(mps, min)))<0 | max(unlist(lapply(mps, min)))>1) {
    stop("Marginal probabilities must be between 0 and 1.")
  }
  cps<-list()
  for(j in 1:length(mps)) {
    cps[[j]]<-mps[[j]][1]
    for(i in 2:(length(mps[[j]])-1)) {
      cps[[j]][i]<-cps[[j]][i-1]+mps[[j]][i]
    }
  }
  return(cps)
}




