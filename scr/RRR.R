
rrr_JRRS = function(Y,X,rank)
{
  n = nrow(Y)
  q = ncol(Y)
  p = ncol(X)
  
  const <- 4
  JRRS <- rep(0, rank)
  CCC <- array(0,dim=c(rank,p,q))
  
  for (ii in 1:rank){
    cat("rank = ",ii,"\n")
    ffiitt <- rrr(Y=Y, X=X, maxrank = rank,ic.type = "GIC")
    CCC[ii,,] <- ffiitt$coef
    r.max <- min(q, p/2, n-1)
    p0.max = min(p/2, n/2)
    S <- norm(Y-X%*%CCC[ii,,], type="F")^2/(n*q-(q + p0.max- r.max)*r.max)
    sigma <- sqrt(S)
    ###JRRS
    MSE <- norm(Y-X%*%CCC[ii,,], type="F")^2
    J <- sum(apply(CCC[ii,,],1,sum)!=0)
    pen <- const*sigma*ii*(2*q+log(2*exp(1))*J+J*log(exp(1)*p/J))
    JRRS[ii] <- MSE + pen
  }
  
  idx <- which(JRRS == min(JRRS), arr.ind = TRUE)
  cat("idx = ",idx,"\n")
  
  
  C0 <- rrr(Y=Y, X=X, penaltySVD = "rank", maxrank  = idx,ic.type = "GIC")$coef
  C.rank = idx
  
  return(list( coef = C0, r= C.rank))
}
