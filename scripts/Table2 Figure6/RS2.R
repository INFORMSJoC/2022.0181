library("lars")
library("MASS")
library("glmnet")


################## sess  ######################
# ??????r fixed
RS2 = function(Y,X,rmax=NULL,r = NULL,outlier = TRUE,thres = NULL)
{
  # new model
  n = nrow(Y)
  #r2 = qr(X)$rank
  if (outlier){#??????outlier??????XI??
    X1 = X
    X = cbind(X1, diag(n)) # ????X=(XI)
  }
  
  p = ncol(X)
  q = ncol(Y) #Y=XB_new+E
  
  # step1 svd on Y
  r1 = qr(Y)$rank
  
  # step 2 choose rank
  if(is.null(r)){
    if(is.null(rmax)){
      rmax = min(r1,n,p,q)
      # strong if else
    }else  
      rmax = min(r1,rmax,n,p,q)
    
    Q = svd(Y,nu = rmax, nv= rmax) 
    Z = Q$u  # Y??????????量????Z=XU?械?Z
    S = diag(Q$d [1:rmax]) 
    V = Q$v # Y??????????量
    
    if(rmax ==1 )
    { ## r_star =1
      S = diag(1)*Q$d [1]
      V.hat = S %*% matrix(V[,1],nrow=1)
    }else 
      V.hat = S%*% t(V)
    
    
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
  }else{
    r_star = r
    Q = svd(Y,nu = r_star, nv= r_star) 
    Z = Q$u  # Y??????????量????Z=XU?械?Z
    S = diag(Q$d [1:r_star]) 
    V = Q$v # Y??????????量
    
    if(r_star ==1 )
    { ## r_star =1
      S = diag(1)*Q$d [1]
      V.hat = S %*% matrix(V[,1],nrow=1)
    }else 
      V.hat = S%*% t(V)
    
  }
  
  
  
  
  # step3 generate U.hat based on lasso with sparsity tunned by BIC
  U.hat=matrix(0,p,r_star)                       
  for(i in 1:r_star){
    fit = glmnet(X,Z[,i], intercept = FALSE)
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
      thres = 0.1 * norm(C.est,"F")
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

# RS2 = function(Y,X,rfix,outlier = TRUE,thres = NULL,alpha)
# {
#   # new model
#   n = nrow(Y)
#   #r2 = qr(X)$rank
#   if (outlier){ 
#     X1 = X
#     X = cbind(X1, diag(n)) # ????X=(XI)
#   }
#   
#   p = ncol(X)
#   q = ncol(Y) #Y=XB_new+E
#   
#   # step1 svd on Y 
#   rmax = rfix
#   
#   Q = svd(Y,nu = rmax, nv= rmax) 
#   Z = Q$u   
#   S = diag(Q$d [1:rmax]) 
#   V = Q$v  
#   
#   if(rmax ==1 )
#   { ## r_star =1
#     S = diag(1)*Q$d [1]
#     V.hat = S %*% matrix(V[,1],nrow=1)
#   }else 
#     V.hat = S%*% t(V)
#   
#    r_star = rfix
#   
#   
#   # step3 generate U.hat based on lasso with sparsity tunned by GIC
#   U.hat=matrix(0,p,r_star)                       
#   for(i in 1:r_star){
#     fit = glmnet(X,Z[,i],alpha =alpha, intercept = FALSE)
#     tLL = fit$nulldev - deviance(fit)
#     k = fit$df
#     n = fit$nobs
#     
#     GIC = log(log(n))*log(p)-tLL
#     step.GIC = which.min(GIC)
#     lam = fit$lambda[step.GIC]
#     U.hat[,i] = coef(fit,lam)[-1]
#     #AICc = tLL + 2*k +2*k*(k+1)/(n-k-1)
#     #BIC = log(n)*k -tLL
#     #aaa = lars(X,Z[,i],type = "lasso",normalize = FALSE, intercept = FALSE,
#     #           use.Gram = FALSE)
#     #BIC = log(n)*aaa$df + n*log(as.vector(aaa$RSS)/n)
#     #step.BIC = which.min(BIC)
#     #cat("parameter = ", step.BIC, "\n")
#     
#   }
#   
#   
#   ## coeffience with best rank r_star
#   if(r_star ==1 )
#   { ## r_star =1
#     coeff = as.matrix(U.hat[,1]) %*% matrix(V.hat[1,],nrow=1)
#   }else 
#     coeff = U.hat[,1:r_star] %*% V.hat[1:r_star,]
#   
#   
#   
#   if(outlier){
#     B.est= coeff[1:(p-n),]
#     C.est = coeff[-(1:(p-n)),]
#     r.est = qr(B.est)$rank
#     if(is.null(thres)){
#       thres = 0.2*norm(C.est,"F")
#       for (kk in 1:n) {
#         if(norm(C.est[kk,],"2")< thres) C.est[kk,]=0
#       }
#     }
#   }else{
#     B.est = coeff
#     C.est = matrix(rep(0,n*q),n,q)
#     r.est = r_star}
#   
#   
#   
#   
#   return(list(rank = r.est, coef=coeff, B  = B.est,C = C.est))
# }