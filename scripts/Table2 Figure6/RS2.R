RS2 = function(Y,X,rfix,outlier = TRUE,thres = NULL,alpha)
{
  # new model
  n = nrow(Y)
  #r2 = qr(X)$rank
  if (outlier){ 
    X1 = X
    X = cbind(X1, diag(n)) # ????X=(XI)
  }
  
  p = ncol(X)
  q = ncol(Y) #Y=XB_new+E
  
  # step1 svd on Y 
  rmax = rfix
  
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
  
   r_star = rfix
  
  
  # step3 generate U.hat based on lasso with sparsity tunned by GIC
  U.hat=matrix(0,p,r_star)                       
  for(i in 1:r_star){
    fit = glmnet(X,Z[,i],alpha =alpha, intercept = FALSE)
    tLL = fit$nulldev - deviance(fit)
    k = fit$df
    n = fit$nobs
    
    GIC = log(log(n))*log(p)-tLL
    step.GIC = which.min(GIC)
    lam = fit$lambda[step.GIC]
    U.hat[,i] = coef(fit,lam)[-1]
    #AICc = tLL + 2*k +2*k*(k+1)/(n-k-1)
    #BIC = log(n)*k -tLL
    #aaa = lars(X,Z[,i],type = "lasso",normalize = FALSE, intercept = FALSE,
    #           use.Gram = FALSE)
    #BIC = log(n)*aaa$df + n*log(as.vector(aaa$RSS)/n)
    #step.BIC = which.min(BIC)
    #cat("parameter = ", step.BIC, "\n")
    
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
      thres = 0.2*norm(C.est,"F")
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