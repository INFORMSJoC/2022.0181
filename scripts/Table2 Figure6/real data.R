################# Arabidopsis thaliana data
################# to store results in Results-Sim23_7_5_Arabidopsis

rm(list = ls())
setwd("d:/Downloads/2022.0181/scripts/Table2 Figure6")

path_result <- "Results-Sim23_7_5_Arabidopsis/"
if (!file_test("-d", path_result)) {
  dir.create(path_result)
}
# load package
library(purrr) # for cross3
library(ncvreg) # for Mest computing
library(snowfall) # parallel programming
library(quadprog)
library(parallel)

X <- read.csv("d:/Downloads/2022.0181/data/file1used.csv")
Y <- read.csv("d:/Downloads/2022.0181/data/file2used.csv")

index = which(Y$Pathwayname %in% c("Plastoquinonebiosynthesis",
                                   "Phytosterolbiosynthesis",
                                   "Carotenoidbiosynthesis",
                                   "PorphyrinChlorophyllmetabolism"))
Y0 = Y[index ,]

Xused = t(X[,7:124]) 
Yused = t(Y0[,5:122])
n=nrow(Yused);q=ncol(Yused);p=ncol(Xused)


Y1 = Yused;Xused = Xused

once_function = function(t)
{
  library(data.table)
  library(RhpcBLASctl)
  blas_set_num_threads(1)
  
  rol = ifelse(t%%100==0,100,t%%100)
  
  # output files
  file_name_Stot <- paste( "Stot.txt", sep = "_")
  file_path_Stot <- paste(path_result, file_name_Stot, sep = "")
  
  file_name_Yreal <- paste( "Yreal.txt", sep = "_")
  file_path_Yreal <- paste(path_result, file_name_Yreal, sep = "")
  
  # training set: 110 samples
  set.seed(rol)
  trainSam = sample(seq_len(nrow(Y1)),100)
  X_train =as.matrix(Xused[trainSam,] ) 
  Y_train = as.matrix(Y1[trainSam,] ) 
  
  # test set: 18 samples
  X_test = as.matrix(Xused[-trainSam,])
  Y_test = as.matrix(Y1[-trainSam,])
  
  Y_true = Y_test
  
  
  source("funs.R")
  library("rrpack")
  library("glmnet")
  library(data.table)
  if(0<t & t<(rep_num+1))
  {
    
    file_name_PE_sss <- paste("PE_sss.txt", sep = "_")
    file_path_PE_sss <- paste(path_result, file_name_PE_sss, sep = "")
    
    file_name_time_sss <- paste( "time_sss.txt", sep = "_")
    file_path_time_sss <- paste(path_result, file_name_time_sss, sep = "")
    
    # lasso 
    set.seed(rol) 
    start1 = Sys.time() 
    re1 = SSS(Y_train, X_train) 
    stop1 = Sys.time()
    timecost1 = stop1 - start1 # time record
    B.est1 = re1$B  
    Ytest1 = X_test%*% B.est1 
    PE_sss = norm(Y_true - Ytest1,"F")/norm(Y_true,"F")
    Sys.sleep(runif(5, min=0, max=1)) 
    data.table::fwrite(as.list(PE_sss), file = file_path_PE_sss, append=TRUE, row.names=FALSE, col.names=FALSE)
    data.table::fwrite(as.list(timecost1), file = file_path_time_sss, append=TRUE, row.names=FALSE, col.names=FALSE)
    
  }else if(t>rep_num & t<(2*rep_num+1)){
    # RRR  
    
    file_name_Sres_rrr <- paste( "Sres_rrr3.txt", sep = "_")
    file_path_Sres_rrr <- paste(path_result, file_name_Sres_rrr, sep = "")
    
    file_name_PE_rrr <- paste("PE_rrr3.txt", sep = "_")
    file_path_PE_rrr <- paste(path_result, file_name_PE_rrr, sep = "")
    
    file_name_time_rrr <- paste("time_rrr3.txt", sep = "_")
    file_path_time_rrr <- paste(path_result, file_name_time_rrr, sep = "")
    
    set.seed(rol)
    start2 = Sys.time() 
    re2 <- cv.rrr(Y_train, X_train, nfold = 10, maxrank = 5) 
    stop2 = Sys.time()
    timecost2 = stop2 - start2 # time record
    B.tilde = coef(re2)  
    Ytest2 = X_test %*% B.tilde
    PE_rrr = norm(Y_true - Ytest2,"F")/norm(Y_true,"F")
    Sys.sleep(runif(5, min=0, max=1))
    data.table::fwrite(as.list(PE_rrr), file = file_path_PE_rrr, append=TRUE, row.names=FALSE, col.names=FALSE)
    data.table::fwrite(as.list(timecost2), file = file_path_time_rrr, append=TRUE, row.names=FALSE, col.names=FALSE)
    
  }else if(t>(2*rep_num) & t<(3*rep_num+1)){
    # R4
    file_name_Sres_r4 <- paste( "Sres_r4.txt", sep = "_")
    file_path_Sres_r4 <- paste(path_result, file_name_Sres_r4, sep = "")
    
    file_name_PE_r4 <- paste("PE_r4.txt", sep = "_")
    file_path_PE_r4 <- paste(path_result, file_name_PE_r4, sep = "")
    
    file_name_time_r4 <- paste("time_r4.txt", sep = "_")
    file_path_time_r4 <- paste(path_result, file_name_time_r4, sep = "")
    
    set.seed(rol)
    source("d:/Downloads/2022.0181/src/r4.1.R")
    start3 = Sys.time()
    re3 = r4.1(Y_train,
               X_train,
               maxrank = 5,
               method = "rowl0",
               Gamma = NULL,
               ic.type ="PIC",
               modstr = list(),
               control = list())
    stop3 = Sys.time()
    timecost3 = stop3 - start3 # time record
    rank.est3 = re3$rank
    B.est3 =  coef(re3)
    C.est3 = re3$s 
    
    Ytest3 = X_test %*% B.est3
    PE_r4 = norm(Y_true - Ytest3,"F")/norm(Y_true,"F") 
    Sys.sleep(runif(5, min=0, max=1))
    data.table::fwrite(as.list(PE_r4), file = file_path_PE_r4, append=TRUE, row.names=FALSE, col.names=FALSE)
    data.table::fwrite(as.list(timecost3), file = file_path_time_r4, append=TRUE, row.names=FALSE, col.names=FALSE)
    
    
  }else if(t>(3*rep_num) & t<(4*rep_num+1)){
    # ROP
    file_name_Sres_rop <- paste("Sres_rop.txt", sep = "_")
    file_path_Sres_rop <- paste(path_result, file_name_Sres_rop, sep = "")
    
    file_name_PE_rop <- paste("PE_rop.txt", sep = "_")
    file_path_PE_rop <- paste(path_result, file_name_PE_rop, sep = "")
    
    file_name_time_rop <- paste("time_rop.txt", sep = "_")
    file_path_time_rop <- paste(path_result, file_name_time_rop, sep = "")
    
    set.seed(rol)
    start4 = Sys.time()
    re4 = RS2(Y_train,X_train,rmax = NULL, r =5,outlier = TRUE,alpha = 1) 
    stop4 = Sys.time()
    timecost4 = stop4 - start4 # time record
    rank.est4 = re4 $rank
    B.est4 = re4$B
    C.est4 = re4$C 
    
    Ytest4 = X_test %*% B.est4
    PE_rop = norm(Y_true - Ytest4,"F")/norm(Y_true,"F")
    Sys.sleep(runif(5, min=0, max=1))
    data.table::fwrite(as.list(PE_rop), file = file_path_PE_rop, append=TRUE, row.names=FALSE, col.names=FALSE)
    data.table::fwrite(as.list(timecost4), file = file_path_time_rop, append=TRUE, row.names=FALSE, col.names=FALSE)
    
  }else if(t>(4*rep_num) & t<(5*rep_num+1)){
    # ROPr 
    file_name_Sres_ropr <- paste("Sres_ropr.txt", sep = "_")
    file_path_Sres_ropr <- paste(path_result, file_name_Sres_ropr, sep = "")
    
    file_name_PE_ropr <- paste("PE_ropr.txt", sep = "_")
    file_path_PE_ropr <- paste(path_result, file_name_PE_ropr, sep = "")
    
    file_name_time_ropr <- paste("time_ropr.txt", sep = "_")
    file_path_time_ropr <- paste(path_result, file_name_time_ropr, sep = "")
    
    set.seed(rol)
    start5 = Sys.time()
    re4 = RS2(Y_train,X_train,rmax = NULL, r =5,outlier = TRUE,alpha = 1) 
    rank.est4 = re4 $rank
    B.est4 = re4$B
    C.est4 = re4$C  
    
    Y1 = Y_train - C.est4
    re5 = RS2(Y1,X_train,rmax = NULL, r =5,outlier = FALSE,alpha = 1) 
    B.est5 = re5$B
    stop5 = Sys.time()
    timecost5 = stop5 - start5 # time record
    
    Ytest5 = X_test %*% B.est5  
    PE_ropr = norm(Y_true - Ytest5,"F")/norm(Y_true,"F")
    Sys.sleep(runif(5, min=0, max=1))
    data.table::fwrite(as.list(PE_ropr), file = file_path_PE_ropr, append=TRUE, row.names=FALSE, col.names=FALSE)
    data.table::fwrite(as.list(timecost5), file = file_path_time_ropr, append=TRUE, row.names=FALSE, col.names=FALSE)
    
  }else if(t>(5*rep_num) & t<(6*rep_num+1)){
    # RCGL 
    file_name_Sres_rcgl <- paste("Sres_rcgl.txt", sep = "_")
    file_path_Sres_rcgl <- paste(path_result, file_name_Sres_rcgl, sep = "")
    
    file_name_PE_rcgl <- paste("PE_rcgl.txt", sep = "_")
    file_path_PE_rcgl <- paste(path_result, file_name_PE_rcgl, sep = "")
    
    file_name_time_rcgl <- paste("time_rcgl.txt", sep = "_")
    file_path_time_rcgl <- paste(path_result, file_name_time_rcgl, sep = "")
    
    source("d:/Downloads/2022.0181/src/RCGL.R")
    set.seed(rol)
    start6 = Sys.time()
    re6 = RCGL(Y_train,X_train,r=5) 
    stop6 = Sys.time()
    rank.est6 = re6 $ r
    B.est6 = re6$coef 
    timecost6 = stop6 - start6 # time record
    
    Ytest6 = X_test %*% B.est6
    PE_rcgl = norm(Y_true - Ytest6,"F")/norm(Y_true,"F")
    Sys.sleep(runif(5, min=0, max=1))
    data.table::fwrite(as.list(PE_rcgl), file = file_path_PE_rcgl, append=TRUE, row.names=FALSE, col.names=FALSE)
    data.table::fwrite(as.list(timecost6), file = file_path_time_rcgl, append=TRUE, row.names=FALSE, col.names=FALSE)
    
  }else if(t>(6*rep_num) & t<(7*rep_num+1)){
    # Secure 
    file_name_Sres_secure <- paste("Sres_secure.txt", sep = "_")
    file_path_Sres_secure <- paste(path_result, file_name_Sres_secure, sep = "")
    
    file_name_PE_secure <- paste("PE_secure.txt", sep = "_")
    file_path_PE_secure <- paste(path_result, file_name_PE_secure, sep = "")
    
    file_name_time_secure <- paste("time_secure.txt", sep = "_")
    file_path_time_secure <- paste(path_result, file_name_time_secure, sep = "")
    
    source("d:/Downloads/2022.0181/src/secure.R")
    set.seed(rol)
    start7 = Sys.time()
    re7 = secure1 (Y_train,X_train, rank =5) 
    B.est7 <- re7$coef
    rank.est7 = re7$rank  
    stop7 = Sys.time()
    timecost7 = stop7 - start7 # time record
    
    Ytest7 = X_test %*% B.est7
    PE_secure = norm(Y_true - Ytest7,"F")/norm(Y_true,"F")
    Sys.sleep(runif(5, min=0, max=1))
    data.table::fwrite(as.list(PE_secure), file = file_path_PE_secure, append=TRUE, row.names=FALSE, col.names=FALSE)
    data.table::fwrite(as.list(timecost7), file = file_path_time_secure, append=TRUE, row.names=FALSE, col.names=FALSE)
    
  }else{
    break
  }
  
}


# snowfall
rep_num = 100
sfInit(parallel = TRUE, cpus = detectCores() - 1) 
sfExport("path_result","Xused","Y1","rep_num")
sfLapply(seq_len(7*rep_num), once_function)
sfStop()

 


