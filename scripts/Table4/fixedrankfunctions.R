## fixed rank methods

RCGL1 = function(Y,X,rank)
{
  n = nrow(Y)
  q = ncol(Y)
  p = ncol(X)
  
  
  C0 <- srrr(Y=Y, X=X, nrank = rank,method = "glasso",ic.type = "GIC")$coef
  C.rank = rank
  
  return(list( coef = C0, r= C.rank))
}


require(secure)
require(glmnet)

#################### only sparse

SSS = function(Y,X)
{
  n = nrow(Y)
  q = ncol(Y)
  p = ncol(X)
  
  B = matrix(nrow = p, ncol = q)
  for(i in 1:q){
    cvob = cv.glmnet(X,Y[,i])
    B[,i]= coef(cvob)[-1]
  }
   
  
  return(list(B = B))
}


require(secure)

############################################################


############################################################
#rank here means the max.rank
secure1 = function(Y,X,rank,nlambda = 100){
  
  # Set largest model to about 25% sparsity
  # See secure.control for setting other parameters
  control <- secure.control(spU=0.25, spV=0.75, elnetAlpha = 1)
  # Complete data case.
  # Fit secure without orthogonality
  fit.orthF <- secure.path(Y,X,nrank=rank,nlambda = nlambda,
                           control=control, orthV = TRUE)
  # fit.orthF <- secure.path(Y,X,nrank=rank.ini,nlambda = nlambda,
  #                          control=control)
  
  CC0 =  fit.orthF$C.est
  rank.C <- qr(CC0)$rank
  
  return(list(coef = CC0, rank = rank.C))
} 


#########################################################
r4.2 <- function(Y,
                 X,
                 maxrank = min(dim(Y), dim(X)),
                 rankinf,
                 method = c("rowl0", "rowl1", "entrywise"),
                 Gamma = NULL,
                 ic.type = c("AIC", "BIC", "PIC"),
                 modstr = list(),
                 control = list())
{
  Call <- match.call()
  
  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)
  
  if (n != nrow(X))
    stop("'nrow(X)' has to be equal to 'nrow(Y)'.")
  
  last.min <- function(a) {
    l <- length(a)
    minid <- which.min(a)
    if (identical(minid, l)) {
      a <- a[l:1]
      minid <-
        l - which(c(a[-length(a)] - a[-1], 0) <= 0)[1] + 1
    }
    minid
  }
  
  control <- do.call("r4.control", control)
  modstr <- do.call("r4.modstr", modstr)
  
  ic <- match.arg(ic.type)
  ic <- switch(ic,
               "AIC" = 1,
               "BIC" = 2,
               "PIC" = 3)
  
  if (!is.null(Gamma)) {
    eigenGm <- eigen(Gamma)
    sqrtGm <- eigenGm$vectors %*%
      diag(sqrt(eigenGm$values)) %*%
      t(eigenGm$vectors)
    sqrtinvGm <- eigenGm$vectors %*%
      diag(sqrt(eigenGm$values ^ (-1))) %*%
      t(eigenGm$vectors)
    Y <- Y %*% sqrtGm
  }
  
  ## xrank <- sum(round(svd(X)$d,4)!=0)
  qrX <- qr(X, tol = control$qr.tol)
  ## C_ls <- qr.coef(qrX, Y)
  ## C_ls <- ifelse(is.na(C_ls), 0, C_ls) ## FIXME
  rX <- qrX$rank
  
  if (is.null(maxrank))
    maxrank <- min(rX, q)
  
  ## number of information criteria
  nIC <- 3
  
  ## Record the Best model for each rank and for each IC
  Barray <- array(dim = c(maxrank, nIC, p, q))
  Carray <- array(dim = c(maxrank, nIC, n, q))
  
  nlam <- modstr$nlam
  
  ## Record all the row l2 norm of the C matrix
  Cmatrix <- array(dim = c(maxrank, nlam, n))
  
  IC <- array(dim = c(maxrank, nlam, nIC))
  smoothIC <- array(dim = c(maxrank, nlam, nIC))
  ## Best one for each rank
  ## and for different information criteria
  Bestid <- matrix(nrow = maxrank, ncol = nIC)
  BestIC <- matrix(nrow = maxrank, ncol = nIC)
  
  SSE <- matrix(nrow = maxrank, ncol = nlam)
  lammatrix <- matrix(nrow = maxrank, ncol = nlam)
  
  ## Compute some matrices
  H0 <- ginv(crossprod(X)) %*% t(X)
  H <- X %*% H0
  
  ## adaptive penalty using leverage. Similar to Owen and She (2011)
  if (modstr$adaptive) {
    weights <- modstr$weights
    if (is.null(weights)) {
      ## Note that even for entrywise penalization,
      ## weights can be a n by 1 vector suppied to hardTH function
      weights <- sqrt(abs(1 - diag(H))) + control$tol
    }
  } else {
    weights <- rep(1, n)
  }
  
  minlam <- modstr$minlam
  maxlam <- modstr$maxlam
  
  
  ## randomly generate some indecies
  delid <- modstr$delid
  if (is.null(delid))
    delid <- sample(seq_len(n), round(n * minlam))
  
  ## Fit model for each rank
  ## check lammax
  lamratio <- 1
  nr <- 1
  time <- vector()
  tj <- 1
  while (nr <= maxrank) {
    ## initial fitting
    ## determine lambda sequence
    ini <- rrr.fit(Y[-delid,], X[-delid,], nrank = nr)
    B0 <- ini$coef
    R <- Y - X %*% B0
    ## l2row  <- hardrowTH(R%*%Gamma,0)$l2
    
    if (method == "rowl0") {
      l2row <-  hardrowTH(R, 0)$l2
      lammax <-  max(l2row / weights) * maxlam * lamratio
      lammin <-  quantile(l2row / weights, 1 - minlam)
    } else if (method == "rowl1") {
      l2row <-  softrowTH(R, 0)$l2
      lammax <-  max(l2row / weights) * maxlam * lamratio
      lammin <-  quantile(l2row / weights, 1 - minlam)
    } else if (method == "entrywise") {
      ## Again, this division should still work for entrywise case
      lammax <-  max(R / weights) * maxlam * lamratio
      lammin <-  quantile(as.vector(R / weights), 1 - minlam)
    }
    lamseq <-  exp(seq(log(lammin), log(lammax), length = nlam))
    
    lammatrix[nr,] <- lamseq
    
    ##Store results for all the lambdas
    Clam <- array(dim = c(nlam, n, q))
    Blam <- array(dim = c(nlam, p, q))
    Clam[1, ,] <- 0
    Blam[1, ,] <- B0
    
    ptm <- proc.time()
    for (lamid in seq_len(nlam)) {
      C0 <- Clam[ifelse(lamid - 1 > 0, lamid - 1, 1), ,]
      B0 <- Blam[ifelse(lamid - 1 > 0, lamid - 1, 1), ,]
      iter <- 1
      diff <- 10 * control$epsilon
      
      while (diff >= control$epsilon & iter < control$maxit) {
        C1 <- C0
        B1 <- B0
        
        R <- Y - X %*% B1
        ## This is based on She(2009)
        ## Cth <- hardrowTH(C1 + (R-C1)%*%Gamma,lamseq[lamid]*weights)
        
        if (method == "rowl0") {
          Cth <- hardrowTH(R, lamseq[lamid] * weights)
          C0 <- Cth$C
        } else if (method == "rowl1") {
          Cth <- softrowTH(R, lamseq[lamid] * weights)
          C0 <- Cth$C
        } else if (method == "entrywise") {
          C0 <- hardThres(R, lamseq[lamid] * weights)
        }
        
        Y0 <- Y - C0
        ## XC <- H %*% Y0 %*% sqrtGm
        XC <- H %*% Y0
        
        ## SVD should be faster then eigen
        svdXC <- tryCatch(
          svd(XC, nu = nr, nv = nr),
          error = function(e)
            2
        )
        if (tryCatch(
          svdXC == 2,
          error = function(e)
            3
        ) == 3) {
          V <- svdXC$v[, seq_len(nr)]
        } else {
          eigenXC <- tryCatch(
            eigen(crossprod(XC)),
            error = function(e)
              2
          )
          if (tryCatch(
            eigenXC == 2,
            error = function(e)
              3
          ) == 3) {
            V <- eigenXC$vectors[, seq_len(nr)]
          }
        }
        
        ## B0 <- H0%*%Y0%*%sqrtGm%*%V%*%t(V)%*%sqrtinvGm
        B0 <- H0 %*% Y0 %*% tcrossprod(V)
        
        C1norm <- sqrt(sum((C1) ^ 2))
        B1norm <- sqrt(sum((B1) ^ 2))
        diff <- ifelse(C1norm == 0, 0,
                       sqrt(sum((C0 - C1) ^ 2)) / C1norm) +
          ifelse(B1norm == 0, 0,
                 sqrt(sum((B0 - B1) ^ 2)) / B1norm)
        iter <- iter + 1
      }
      
      Clam[lamid, ,] <- C0
      Blam[lamid, ,] <- B0
      
      Cmatrix[nr, lamid,] <- apply(C0 ^ 2, 1, sum)
      
      SSE[nr, lamid] <- sum((Y - C0 - X %*% B0) ^ 2)
      
      ## FIXME?
      if (method == "rowl0" | method == "rowl1") {
        nout <- sum(apply(C0 != 0, 1, sum) != 0)
        df <- (rX + q - nr) * nr + nout * q
      } else if (method == "entrywise") {
        nout <- sum(C0 != 0)
        df <- (rX + q - nr) * nr + sum(C0 != 0)
      }
      if (method == "rowl0" | method == "rowl1") {
        nout <- sum(apply(C0 != 0, 1, sum) != 0)
        IC[nr, lamid, 1] <-
          n * q * log(SSE[nr, lamid] / n / q) + 2 * df
        IC[nr, lamid, 2] <-
          n * q * log(SSE[nr, lamid] / n / q) + log(q * n) * df
        IC[nr, lamid, 3] <- n * q * log(SSE[nr, lamid] / n / q) +
          7 * df + ifelse(nout == 0, 0, 2.1 * nout * log(exp(1) *
                                                           n / nout))
      } else if (method == "entrywise") {
        nout <- sum(C0 != 0)
        IC[nr, lamid, 1] <-
          n * q * log(SSE[nr, lamid] / n / q) + 2 * df
        IC[nr, lamid, 2] <-
          n * q * log(SSE[nr, lamid] / n / q) + log(q * n) * df
        ## Need update
        IC[nr, lamid, 3] <- n * q * log(SSE[nr, lamid] / n / q) +
          7 * df + 2 * log(choose(n * q, nout))
      }
    }
    
    time[tj] <- (proc.time() - ptm)[3]
    tj <- tj + 1
    
    ## check whether lammax is too small
    if (sum(apply(Clam[nlam, , ], 1, sum) != 0) > n * 0.05) {
      lamratio <- lamratio + 1
    } else{
      ## lamratio <- 1
      
      ## smooth IC curves then select a model in favor of less outliers
      ## Select the best fit for the given rank
      for (i in seq_len(nIC)) {
        smoothIC[nr, , i] <- loess(IC[nr, , i] ~
                                     c(seq_len(nlam)))$fitted
      }
      Bestid[nr, ] <- apply(smoothIC[nr, , ], 2, last.min)
      ## Bestid[nr,] <- apply(smoothIC[nr,,],2,which.min)
      ## Bestid[nr,] <- apply(IC[nr,,],2,which.min)
      
      for (i in 1:nIC) {
        BestIC[nr, i] <- smoothIC[nr, Bestid[nr, i], i]
        ## BestIC[nr,2] <- smoothIC[nr,Bestid[nr,2],2]
        ## BestIC[nr,3] <- smoothIC[nr,Bestid[nr,3],3]
        if (is.null(Gamma)) {
          Carray[nr, i, , ] <- Clam[Bestid[nr, i], , ]
          Barray[nr, i, , ] <- Blam[Bestid[nr, i], , ]
        } else {
          Carray[nr, i, , ] <- Clam[Bestid[nr, i], , ] %*% sqrtinvGm
          Barray[nr, i, , ] <-
            Blam[Bestid[nr, i], , ] %*% sqrtinvGm
        }
        ## Barray[nr,2,,] <- Blam[Bestid[nr,2],,]
        ## Barray[nr,3,,] <- Blam[Bestid[nr,3],,]
      }
      nr <- nr + 1
    }
  }
  ## Select the best rank
  Bestrank <- apply(BestIC, 2, which.min)
  
  Bestrank = rankinf
  
  coef.b = Barray[Bestrank[ic], ic, , ]
  coef.c = Carray[Bestrank[ic], ic, , ]
  
  out <- list(
    call = Call,
    coef.path = Barray,
    s.path = Carray,
    s.norm.path = Cmatrix,
    ic.path = IC,
    ic.smooth.path = smoothIC,
    lambda.path = lammatrix,
    id.solution = Bestid,
    ic.best = BestIC,
    rank.best = Bestrank,
    ## time = time,
    ## sse = SSE,
    ## lamratio=lamratio
    coef = coef.b,
    s = coef.c,
    rank = Bestrank[ic]
  )
  class(out) <- "r4"
  out
}



##' Fitting reduced-rank ridge regression with given rank and shrinkage penalty
##'
##' Fitting reduced-rank ridge regression with given rank and shrinkage penalty
##'
##' @param Y a matrix of response (n by q)
##' @param X a matrix of covariate (n by p)
##' @param nrank an integer specifying the desired rank
##' @param lambda tunging parameter for the ridge penalty
##' @param coefSVD logical indicating the need for SVD for the
##'   coeffient matrix int the output
##' @return S3 \code{rrr} object, a list consisting of
##'   \item{coef}{coefficient of rrs}
##'   \item{coef.ls}{coefficient of least square}
##'   \item{fitted}{fitted value of rrs}
##'   \item{fitted.ls}{fitted value of least square}
##'   \item{A}{right singular matrix}
##'   \item{Ad}{sigular value vector}
##'   \item{nrank}{rank of the fitted rrr}
##' @examples
##' library(rrpack)
##' Y <- matrix(rnorm(400), 100, 4)
##' X <- matrix(rnorm(800), 100, 8)
##' rfit <- rrs.fit(Y, X)
##' @references
##'
##' Mukherjee, A. and Zhu, J. (2011) Reduced rank ridge regression and its
##' kernal extensions.
##'
##' Mukherjee, A., Chen, K., Wang, N. and Zhu, J. (2015) On the degrees of
##' freedom of reduced-rank estimators in multivariate
##' regression. \emph{Biometrika}, 102, 457--477.
##'
##' @importFrom MASS ginv
##' @export
rrs.fit <- function(Y,
                    X,
                    nrank = min(ncol(Y), ncol(X)),
                    lambda = 1,
                    coefSVD = FALSE)
{
  Call <- match.call()
  
  q <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  
  S_yx <- t(Y) %*% X
  ## This is a key difference
  S_xx <- t(X) %*% X + lambda * diag(p)
  
  ## S_xx_inv <- tryCatch(solve(S_xx+0.1*diag(p)),error=function(e)ginv(S_xx))
  ## S_xx_inv <- ginv(S_xx)
  ## Use the Woodbury matrix identity
  if (lambda != 0) {
    S_xx_inv <- 1 / lambda * diag(p) -
      lambda ^ (-2) * t(X) %*% ginv(diag(n) + lambda ^ (-1) * X %*%
                                      t(X)) %*% X
  } else{
    S_xx_inv <- ginv(S_xx)
    if (sum(is.na(S_xx_inv)) > 0) {
      S_xx_inv <- solve(S_xx + 0.1 * diag(p))
    }
  }
  
  C_ls <- S_xx_inv %*% t(S_yx)
  
  ypy.svd <- TRUE
  ##if(ypy.svd){
  ##This is another key difference
  XC <- rbind(X, sqrt(lambda) * diag(p)) %*% C_ls
  svdXC <- tryCatch(
    svd(XC, nu = nrank, nv = nrank),
    error = function(e)
      2)
  if (tryCatch(
    svdXC == 2,
    error = function(e)
      3) == 3) {
    A <- svdXC$v[, 1:nrank]
    Ad <- (svdXC$d[1:nrank]) ^ 2
  } else{
    ypy.svd <- FALSE
  }
  #}
  if (!ypy.svd) {
    SS <- S_yx %*% C_ls
    SS <- (SS + t(SS)) / 2
    eigenSS <- eigen(SS, symmetric = TRUE)
    A <- as.matrix(eigenSS$vectors[, 1:nrank])
    Ad <- eigenSS$values[1:nrank]
  }
  
  AA <- A %*% t(A)
  C_rr <- C_ls %*% AA
  
  ##    if(c.svd){
  ##      svd_C <- svd(C_rr,nv=nrank,nu=nrank)
  ##      U <- as.matrix(svd_C$u[,1:nrank])
  ##      V <- as.matrix(svd_C$v[,1:nrank])
  ##      D <- diag(svd_C$d[1:nrank],nrow=nrank)
  ##
  ##      ####return ls estimator C_ls, reduced-rank estimator C_rr
  ##      ####return SVD of C_rr
  ##      list(A=A,Ad=Ad,C_ls=C_ls,C_rr=C_rr,U=U,V=V,D=D,C=C_rr,rank=nrank)
  ##    }else{
  ##      list(A=A,Ad=Ad,C_ls=C_ls,C_rr=C_rr,C=C_rr,rank=nrank)
  ##    }
  
  ret <- list(
    call = Call,
    coef = C_rr,
    coef.ls = C_ls,
    fitted = X %*% C_rr,
    fitted.ls = XC,
    A = A,
    Ad = Ad,
    nrank = nrank
  )
  
  if (coefSVD) {
    coefSVD <- svd(C_rr, nrank, nrank)
    coefSVD$d <- coefSVD$d[1:nrank]
    ret <- c(ret, list(coefSVD = coefSVD))
  }
  
  class(ret) <- "rrs.fit"
  ret
}


## Internal function for specifying model parameters
##
## a list of internal model parameters controlling the model fitting
##
## @param nlam parameter in the augmented Lagrangian function
## @param adaptive if TRUE, use leverage values for adaptive penalization
## @param weights user supplied weights for adaptive penalization
## @param minlam maximum proportion of outliers
## @param maxlam maximum proportion of good observations
## @param delid discarded observation indices for initial estimation
##
## @return a list of model parameters.
r4.modstr <- function(nlam = 100,
                      adaptive = TRUE,
                      weights = NULL,
                      minlam = 0.3,
                      maxlam = 1,
                      delid = NULL)
{
  list(
    nlam = nlam,
    adaptive = adaptive,
    weights = weights,
    minlam = minlam,
    maxlam = maxlam,
    delid = delid
  )
}


## Internal function for specifying computation parameters
##
## a list of internal computational parameters controlling optimization
##
## @param epsilon convergence tolerance
## @param maxit maximum number of iterations
## @param qr.tol tolerance for qr decomposition
## @param tol tolerance
##
## @return a list of computational parameters.
r4.control <- function(epsilon = 1e-3,
                       maxit = 100L,
                       qr.tol = 1e-4,
                       tol = 1e-06)
{
  list(
    epsilon = epsilon,
    maxit = maxit,
    qr.tol = qr.tol,
    tol = tol
  )
}

## hard thresholding function
hardThres <- function(x, lambda) {
  ##  ifelse(abs(x) > lambda, x, 0) ## nice but too slow
  xt <- abs(x) - lambda
  x[xt <= 0] <- 0
  x
}


## ##Hard thresholding for rows
## ##lambda can be a vector
hardrowTH <- function(C, lambda, eta = 0) {
  l2row <- apply(C, 1, function(x)
    sqrt(sum(x ^ 2)))
  C <- C / (1 + eta)
  C[l2row <= lambda, ] <- 0
  list(C = C, l2 = l2row)
}


################################################
##############################################################################

library("lars")
library("MASS")
library("glmnet")


################## sess  ######################
# ??????r fixed
RS2 = function(Y,X,rmax=NULL,r = NULL,outlier = TRUE,thres = NULL,alpha)
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
    fit = glmnet(X,Z[,i],alpha = alpha, intercept = FALSE)
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

