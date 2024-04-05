
##############################################################################

library("lars")
library("MASS")
library("glmnet")


################## sess  ######################
RS1 = function(Y,X,r=NULL,outlier = TRUE,thres = NULL)
{
  # new model
  n = nrow(Y)
  #r2 = qr(X)$rank
  if (outlier){
    X1 = X
    X = cbind(X1, diag(n)) 
  }
  
  p = ncol(X)
  q = ncol(Y) #Y=XB_new+E
  
  # step1 svd on Y
  r1 = qr(Y)$rank
  
  
  if(is.null(r)){
    rmax = min(r1,n,p,q)
    # strong if else
  }else  
    rmax = min(r1,r,n,p,q)
  
  Q = svd(Y,nu = rmax, nv= rmax) 
  Z = Q$u  
  S = diag(Q$d [1:rmax]) 
  V = Q$v 
  
  if(rmax ==1 )
  { ## r_star =1
    S = diag(1)*Q$d [1]
    V.hat = S %*% matrix(V[,1],nrow=1)
  }else 
    V.hat = S%*% t(V)
  
  
  # step 2 choose rank
  cn=rep(0,rmax)
  lambda=0
  # sign=0
  stop=(10^(-2))*((norm(Y,"F"))^(2))*((n*q)^(-1))/10000  # Frobenius norm,stop for interruption
  # accu=rep(0,(r+4))
  for (i in 1 : rmax){
    if((lambda>stop) || (lambda==0)){
      ZSV = as.matrix(Z[,1:i]) %*% matrix(V.hat[1:i,],nrow = i)
      layer = Z[,i] %*% t(V.hat[i,])
      lambda = ( norm( layer, "F" ) ) * ( sqrt ( n*q ) ^ (-1) )
      # accu[i]=norm(Y - X %*% B.hat,"F")
      cn[i] = sqrt(n) * log( (norm(Y - ZSV,"F") ^ 2) / (n*q) ) +  log(n)* i 
    }
  }
  #print(cn)
  r_star = which.min(cn)
  
  
  # step3 generate U.hat based on lasso with sparsity tunned by BIC
  U.hat=matrix(0,p,r_star)                       
  for(i in 1:r_star){
    fit = glmnet(X,Z[,i],intercept = FALSE)
    tLL = fit$nulldev - deviance(fit)
    k = fit$df
    n = fit$nobs
  
    GIC = log(log(n))*log(p)-tLL
    step.GIC = which.min(GIC)
    lam = fit$lambda[step.GIC]
    U.hat[,i] = coef(fit,lam)[-1]
    
  }
  
  
  ## coeffience with best rank r_star
  if(r_star ==1 )
  { ## r_star =1
    coeff = as.matrix(U.hat[,1]) %*% matrix(V.hat[1,],nrow=1)
  }else 
    coeff = U.hat[,1:r_star] %*% V.hat[1:r_star,]
  
  
  
  if(outlier){
    B.est= coeff[1:(p-n),]
    C.est = coeff[-(1:(p-n)),]
    r.est = qr(B.est)$rank
    if(is.null(thres)){
      thres = sqrt(n) + sqrt(q)
      for (kk in 1:n) {
        if(norm(C.est[kk,],"2")< thres) C.est[kk,]=0
      }
    }
  }else{
    B.est = coeff
    C.est = matrix(rep(0,n*q),n,q)
    r.est = r_star}


  
  
  return(list(rank = r.est, coef=coeff, B  = B.est,C = C.est))
}
