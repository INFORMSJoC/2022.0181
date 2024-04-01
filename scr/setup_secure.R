############################################
##' Simulation model 2
##'
##' Generate data from secure + r4(outlier setting)
##'
##' @usage
##' setsecure(n = 400, p = 100, q = 50, nrank = 5, rho_X = 0,
##'          rho_E = 0, nout = 10, vout = NULL, voutsd = 2, nlev = 10,
##'          vlev = 10, vlevsd = NULL, spaB = 0.002, spaC = 0.01, 
##'          SigmaX = CorrCS, SigmaE = CorrCS)
##'
##' @param n sample size
##' @param p number of predictors
##' @param q numbers of responses
##' @param nrank B rank: should be 3 as set in secure
##' @param s = c(s1,s2,s3) sparsity of u:(8,9,9) or (16,18,18)
##' @param rho_X correlation parameter for predictors
##' @param rho_E correlation parameter for errors
##' @param snr signal noise ratio
##' @param nout number of outliers; should be smaller than n
##' @param vout value of outliers
##' @param voutsd control mean-shifted magnitude of outliers, 2 or 4 here
##' @param nlev number of high-leverage outliers
##' @param vlev value of leverage
##' @param vlevsd control magnitude of leverage
##' @return similated model and data
##'
##' @references
##' Mishra et al. (2017) Sequential Co-Sparse Factor Regression. \emph{JCGS}, 26:4, 814-825.
##' She, Y. and Chen, K. (2017) Robust reduced-rank regression. \emph{Biometrika}, 104 (3), 633--647.
##' @importFrom stats sd
##' @export
setsecure =
  function(n = 400,
           p = 500,
           q = 200,
           s = c(8,9,9),
           nrank = 3,
           rho_X = 0.5,
           rho_E = 0.5,
           snr = 0.75,
           nout = 5,
           vout = NULL,
           voutsd = 2,
           nlev = 0,
           vlev = NULL,
           vlevsd = NULL)
{
    require(secure)
    # Simulate data from a sparse factor regression model
    s1 = s[1]; s2=s[2]; s3=s[3]
    U <- matrix(0,ncol=nrank ,nrow=p); V <- matrix(0,ncol=nrank ,nrow=q)
    
    U[,1]<-c(sample(c(1,-1),s1,replace=TRUE),rep(0,(p-s1)))
    U[,2]<-c(rep(0,5),sample(c(1,-1),s2,replace=TRUE),rep(0,(p-5-s2)))
    U[,3]<-c(rep(0,11),sample(c(1,-1),s3,replace=TRUE),rep(0,(p-11-s3)))
    
    V[,1]<-c(sample(c(1,-1),5,replace=TRUE)*runif(5,0.3,1),rep(0,q-5))
    V[,2]<-c(rep(0,5),sample(c(1,-1),5,replace=TRUE)*runif(5,0.3,1),rep(0,q-10))
    V[,3]<-c(rep(0,10),sample(c(1,-1),5,replace=TRUE)*runif(5,0.3,1),rep(0,q-15))
    
    U[,1:3]<- apply(U[,1:3],2,function(x)x/sqrt(sum(x^2)))
    V[,1:3]<- apply(V[,1:3],2,function(x)x/sqrt(sum(x^2)))
    
    D <- diag(c(60,30,10))
    B <- U%*%D%*%t(V)
    
    Xsigma <- rho_X^abs(outer(1:p, 1:p,FUN="-"))
    sim.sample <- secure.sim(U,D,V,n,snr = snr,Xsigma,rho =  rho_E)
    X <- sim.sample$X
    Y = sim.sample$Y
    
    Ymean = X %*% B
    E = Y - Ymean

    
    
    ## set C (sparse)
    if (nout != 0) {
      if (is.null(vout)) {
        Ysd <- apply(Ymean, 2, sd) * sample(c(-1, 1), q, replace = TRUE)
        C <- voutsd * matrix(nrow = nout,
                             ncol = q,
                             byrow = T,
                             Ysd)
        
        Y[1:nout, ] <-  Y[1:nout, ] + C
        
        #Y = Y + C
        
      } else{
        Vout <- vout * sample(c(-1, 1), q, replace = TRUE)
        C <- matrix(nrow = nout,
                    ncol = q,
                    byrow = T,
                    Vout)
        Y[1:nout, ] <- Y[1:nout, ] + C
      }}
    
    if (nlev != 0) {
      if (is.null(vlev)) {
        Xsd <- apply(X, 2, sd) * sample(c(-1, 1), p, replace = TRUE)
        Xlev <- vlevsd * matrix(nrow = nlev,
                                ncol = p,
                                byrow = T,
                                Xsd)
        X[1:nlev, ] <-  Xlev
      } else{
        Vlev <- vlev * sample(c(-1, 1), p, replace = TRUE)
        Xlev <- matrix(nrow = nlev,
                       ncol = p,
                       byrow = T,
                       Vlev)
        X[1:nlev, ] <- Xlev
      }
    }
    
    C1 = C
    C0 = matrix(rep(0 , (n-nout) * q) , (n-nout) , q)
    C = rbind(C1, C0)
    
    return(list(Y=Y,B=B,C=C,E=E,X=X))
    
  
}









