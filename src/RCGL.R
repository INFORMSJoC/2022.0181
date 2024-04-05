
require(secure)
library("rrpack")


srrr.rank <- function(Y, X, sigma, method = c("glasso", "adglasso"), r.max = min(p,q)){
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  
  const <- 4
  JRRS <- matrix(0, r.max, 100)
  colnames(JRRS) <- 1:100
  rownames(JRRS) <- 1:r.max
  for(r in 1:r.max){
    fit <- srrr(Y=Y, X=X, nrank = r, ic.type = "BIC", method = method)
    C1 <- array(0, dim = c(length(fit$lambda), p ,q))
    for(i in 1:length(fit$lambda)){
      C1[i,,] <- fit$A.path[i,,]%*% t(fit$V.path[i,,])
    }
    MSE <- apply(C1, 1, function(C0) norm(Y-X%*%C0, type="F")^2)
    J <- apply(C1, 1, function(C0) sum(apply(C0,1,sum)!=0))
    pen <- const*sigma*r*(2*q+log(2*exp(1))*J+J*log(exp(1)*p/J))
    JRRS[r,] <- MSE + pen
  }
  idx <- which(JRRS == min(JRRS), arr.ind = TRUE)
  
  fit1 <- srrr(Y=Y, X=X, nrank = idx[1], modstr = list(lamA = fit$lambda[idx[2]]), method = method)
  C <- fit1$coef
  out <- list(Y = Y, X = X, C = C, JRRS = JRRS, nrank = idx[1], lambda = fit$lambda[idx[2]])
  return(out)
}





  
  ##################################  RCGL  ######################################## 
 

  ##########################################################################
# here the rank means max.rank
RCGL = function(Y,X,rank)
{
  n = nrow(Y)
  q = ncol(Y)
  p = ncol(X)
  
  const <- 4
  JRRS <- rep(0, rank)
  CCC <- array(0,dim=c(rank,p,q))
  
  for (ii in 1:rank){
    cat("rank = ",ii,"\n")
    ffiitt <- srrr(Y=Y, X=X, nrank = ii, method = "glasso",ic.type = "GIC")
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
  
  
  C0 <- srrr(Y=Y, X=X, nrank = idx,method = "glasso",ic.type = "GIC")$coef
  C.rank = idx
  
  return(list( coef = C0, r= C.rank))
}
  
  
 


