############################################
#                                          #
#         secure--secure setting           #
#                                          #
############################################
require(secure)

#rank here means the max.rank
secure1 = function(Y,X,rank,nlambda = 100){
  
  # Set largest model to about 25% sparsity
  # See secure.control for setting other parameters
  control <- secure.control(spU=0.25, spV=0.25, elnetAlpha = 1)
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

  
 
