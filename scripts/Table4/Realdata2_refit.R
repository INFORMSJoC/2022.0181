# get resuls stored in  "Results-Sim23_7_5_rollingwindow/"

rm(list = ls())
setwd("d:/Downloads/2022.0181/scripts/Table4") 
########### read table
CatAused = read.csv("d:/Downloads/2022.0181/data/itemchosenbyA.csv")
CatBused = read.csv("d:/Downloads/2022.0181/data/itemchosenbyB.csv")
CatCused = read.csv("d:/Downloads/2022.0181/data/itemchosenbyC.csv")
CatDused = read.csv("d:/Downloads/2022.0181/data/itemchosenbyD.csv")
CatEused = read.csv("d:/Downloads/2022.0181/data/itemchosenbyE.csv")

MATtol = cbind(CatAused[,-1],CatBused[,-1],CatCused[,-1],CatDused[,-1],CatEused[,-1]) 
head(MATtol)

dataused = MATtol[-216,]

Yused = dataused[-c(1:5),]
XL1 = dataused[c(5:214),]
XL2 = dataused[c(4:213),]
XL3 = dataused[c(3:212),]
XL4 = dataused[c(2:211),]
XL5 = dataused[c(1:210),]
XX1 = cbind(XL1,XL2,XL3,XL4,XL5)
qsum = ncol(Yused)
Lsum = 5
Xused = matrix(nrow = nrow(Yused),ncol = qsum*Lsum)
for (ll in 1:qsum) {
  for (kk in 1:Lsum) {
    # print(c(ll,kk))
    mm = (ll-1)*Lsum+kk
    tt = qsum*(kk-1)+ll
    # print(c(mm,tt))
    Xused[,mm] = XX1[,tt]
  }
}

Yused = as.matrix(Yused)
cor = cor(Yused)
#print(cor)
is.matrix(Xused)

# ROP
# rolling window 

h = 120
#tt in 1: 210-h
path_result <- "Results-Sim23_7_5_rollingwindow/"
if (!file_test("-d", path_result)) {
  dir.create(path_result)
}



set.seed(2023)
source("fixedrankfunctions.R") 
re1 = RS2(Yused,Xused,rmax = NULL, r =8,outlier = TRUE,alpha = 0.6) 
rank.est1 = re1 $rank
B.est1 = re1 $B
C.est1 = re1 $C
# outlier analysis
# order(apply(abs(C.est1), 1, sum),decreasing=TRUE)[1:10]
index = order(apply(abs(C.est1), 1, sum),decreasing=TRUE)[1:10]
length(which(apply(C.est1,1,sd)!=0))
which(apply(C.est1,1,sd)!=0)%% 24   
P1 = hist(which(apply(C.est1,1,sd)!=0)%% 24,main = "Outlier distribution",xlab = "The hour",ylab = "The number of outliers") 
# refit B and analyze U 
Y.modified1 = Yused - C.est1

Xused = Xused[-index, ]
Yused = Y.modified1[-index, ]

once_function = function(t)
{
  library(data.table)
  library(RhpcBLASctl)
  blas_set_num_threads(1)
  
  hour = ifelse(t%%rol==0,rol,t%%rol)
  
  print(paste0("----------", hour, "-th hour ----------"))
  # output files
  file_name_Stot <- paste( "Stot.txt", sep = "_")
  file_path_Stot <- paste(path_result, file_name_Stot, sep = "")
  
  file_name_Yreal <- paste( "Yreal.txt", sep = "_")
  file_path_Yreal <- paste(path_result, file_name_Yreal, sep = "")
  
  # training set contain past h timestamps
  X_train =as.matrix(Xused[c(hour:(hour+h - 1)),] ) 
  Y_train = as.matrix(Yused[c(hour:(hour+h - 1)),] ) 
  
  # test set 
  X_test = as.matrix(t(Xused[(hour+h),]))
  Y_test = as.matrix(t(Yused[(hour+h),]))
  
  # sample mean
  Y_bar = as.matrix( t(apply(Y_train,2,mean)))
  Y_true = Y_test
  Stot =  (Y_true - Y_bar)^2 
  
  
  source("fixedrankfunctions.R")
  library("rrpack")
  library("glmnet")
  library(data.table)
  if(0<t & t<(rol+1))
  {
    #cat(Stot, "\n", file = file_path_Stot, append = TRUE)
    data.table::fwrite(as.list(Stot), file = file_path_Stot, append=TRUE, row.names=FALSE, col.names=FALSE)
    data.table::fwrite(as.list(Y_true), file = file_path_Yreal, append=TRUE, row.names=FALSE, col.names=FALSE)
    
    file_name_Sres_sss <- paste("Sres_sss.txt", sep = "_")
    file_path_Sres_sss <- paste(path_result, file_name_Sres_sss, sep = "")
    
    file_name_PE_sss <- paste("PE_sss.txt", sep = "_")
    file_path_PE_sss <- paste(path_result, file_name_PE_sss, sep = "")
    
    file_name_time_sss <- paste( "time_sss.txt", sep = "_")
    file_path_time_sss <- paste(path_result, file_name_time_sss, sep = "")
    
    # lasso 
    set.seed(t) 
    start1 = Sys.time() 
    re1 = SSS(Y_train, X_train) 
    stop1 = Sys.time()
    timecost1 = stop1 - start1 # time record
    B.est1 = re1$B  
    Ytest1 = X_test%*% B.est1
    Sres_sss =   (Y_true - Ytest1)^2 
    PE_sss = norm(Y_true - Ytest1,"2")/norm(Y_true,"2")
    Sys.sleep(runif(5, min=0, max=1)) 
    data.table::fwrite(as.list(Sres_sss), file = file_path_Sres_sss, append=TRUE, row.names=FALSE, col.names=FALSE)
    data.table::fwrite(as.list(PE_sss), file = file_path_PE_sss, append=TRUE, row.names=FALSE, col.names=FALSE)
    data.table::fwrite(as.list(timecost1), file = file_path_time_sss, append=TRUE, row.names=FALSE, col.names=FALSE)
    
  }else if(t>rol & t<(2*rol+1)){
    # RRR  
    
    file_name_Sres_rrr <- paste( "Sres_rrr3.txt", sep = "_")
    file_path_Sres_rrr <- paste(path_result, file_name_Sres_rrr, sep = "")
    
    file_name_PE_rrr <- paste("PE_rrr3.txt", sep = "_")
    file_path_PE_rrr <- paste(path_result, file_name_PE_rrr, sep = "")
    
    file_name_time_rrr <- paste("time_rrr3.txt", sep = "_")
    file_path_time_rrr <- paste(path_result, file_name_time_rrr, sep = "")
    
    set.seed(t)
    start2 = Sys.time() 
    re2 <- cv.rrr(Y_train, X_train, nfold = 10, maxrank = 8) 
    stop2 = Sys.time()
    timecost2 = stop2 - start2 # time record
    B.tilde = coef(re2)  
    Ytest2 = X_test %*% B.tilde
    Sres_rrr =   (Y_true - Ytest2)^2
    PE_rrr = norm(Y_true - Ytest2,"2")/norm(Y_true,"2")
    Sys.sleep(runif(5, min=0, max=1))
    #cat(Sres_rrr, "\n", file = file_path_Sres_rrr, append = TRUE)
    data.table::fwrite(as.list(Sres_rrr), file = file_path_Sres_rrr, append=TRUE, row.names=FALSE, col.names=FALSE)
    data.table::fwrite(as.list(PE_rrr), file = file_path_PE_rrr, append=TRUE, row.names=FALSE, col.names=FALSE)
    data.table::fwrite(as.list(timecost2), file = file_path_time_rrr, append=TRUE, row.names=FALSE, col.names=FALSE)
    
  }else if(t>(2*rol) & t<(3*rol+1)){
    # R4
    file_name_Sres_r4 <- paste( "Sres_r4.txt", sep = "_")
    file_path_Sres_r4 <- paste(path_result, file_name_Sres_r4, sep = "")
    
    file_name_PE_r4 <- paste("PE_r4.txt", sep = "_")
    file_path_PE_r4 <- paste(path_result, file_name_PE_r4, sep = "")
    
    file_name_time_r4 <- paste("time_r4.txt", sep = "_")
    file_path_time_r4 <- paste(path_result, file_name_time_r4, sep = "")
    
    set.seed(t)
    source("d:/Downloads/2022.0181/src/r4.1.R")
    start3 = Sys.time()
    re3 = r4.1(Y_train,
               X_train,
               maxrank = 8,
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
    Sres_r4 =   (Y_true - Ytest3)^2 
    PE_r4 = norm(Y_true - Ytest3,"2")/norm(Y_true,"2") 
    Sys.sleep(runif(5, min=0, max=1))
    #cat(Sres_r4, "\n", file = file_path_Sres_r4, append = TRUE)
    data.table::fwrite(as.list(Sres_r4), file = file_path_Sres_r4, append=TRUE, row.names=FALSE, col.names=FALSE)
    data.table::fwrite(as.list(PE_r4), file = file_path_PE_r4, append=TRUE, row.names=FALSE, col.names=FALSE)
    data.table::fwrite(as.list(timecost3), file = file_path_time_r4, append=TRUE, row.names=FALSE, col.names=FALSE)
    
    
  }else if(t>(3*rol) & t<(4*rol+1)){
    # RCGL 
    file_name_Sres_rcgl <- paste("Sres_rcgl.txt", sep = "_")
    file_path_Sres_rcgl <- paste(path_result, file_name_Sres_rcgl, sep = "")
    
    file_name_PE_rcgl <- paste("PE_rcgl.txt", sep = "_")
    file_path_PE_rcgl <- paste(path_result, file_name_PE_rcgl, sep = "")
    
    file_name_time_rcgl <- paste("time_rcgl.txt", sep = "_")
    file_path_time_rcgl <- paste(path_result, file_name_time_rcgl, sep = "")
    
    source("d:/Downloads/2022.0181/src/RCGL.R")
    set.seed(t)
    start6 = Sys.time()
    re6 = RCGL(Y_train,X_train,r=8) 
    stop6 = Sys.time()
    rank.est6 = re6 $ r
    B.est6 = re6$coef 
    timecost6 = stop6 - start6 # time record
    
    Ytest6 = X_test %*% B.est6
    Sres_rcgl =   (Y_true - Ytest6)^2  
    PE_rcgl = norm(Y_true - Ytest6,"2")/norm(Y_true,"2")
    Sys.sleep(runif(5, min=0, max=1))
    #cat(Sres_rcgl, "\n", file = file_path_Sres_rcgl, append = TRUE)
    data.table::fwrite(as.list(Sres_rcgl), file = file_path_Sres_rcgl, append=TRUE, row.names=FALSE, col.names=FALSE)
    data.table::fwrite(as.list(PE_rcgl), file = file_path_PE_rcgl, append=TRUE, row.names=FALSE, col.names=FALSE)
    data.table::fwrite(as.list(timecost6), file = file_path_time_rcgl, append=TRUE, row.names=FALSE, col.names=FALSE)
    
    
    
  }else if(t>(4*rol) & t<(5*rol+1)){
    # Secure 
    file_name_Sres_secure <- paste("Sres_secure.txt", sep = "_")
    file_path_Sres_secure <- paste(path_result, file_name_Sres_secure, sep = "")
    
    file_name_PE_secure <- paste("PE_secure.txt", sep = "_")
    file_path_PE_secure <- paste(path_result, file_name_PE_secure, sep = "")
    
    file_name_time_secure <- paste("time_secure.txt", sep = "_")
    file_path_time_secure <- paste(path_result, file_name_time_secure, sep = "")
    
    source("d:/Downloads/2022.0181/src/secure.R")
    set.seed(t)
    start7 = Sys.time()
    re7 = secure1 (Y_train,X_train, rank =8) 
    B.est7 <- re7$coef
    rank.est7 = re7$rank  
    stop7 = Sys.time()
    timecost7 = stop7 - start7 # time record
    
    Ytest7 = X_test %*% B.est7
    Sres_secure =   (Y_true - Ytest7)^2  
    PE_secure = norm(Y_true - Ytest7,"2")/norm(Y_true,"2")
    Sys.sleep(runif(5, min=0, max=1))
    #cat(Sres_secure, "\n", file = file_path_Sres_secure, append = TRUE)
    data.table::fwrite(as.list(Sres_secure), file = file_path_Sres_secure, append=TRUE, row.names=FALSE, col.names=FALSE)
    data.table::fwrite(as.list(PE_secure), file = file_path_PE_secure, append=TRUE, row.names=FALSE, col.names=FALSE)
    data.table::fwrite(as.list(timecost7), file = file_path_time_secure, append=TRUE, row.names=FALSE, col.names=FALSE)
    
  }else{
    break
  }
  
}


# snowfall
sfInit(parallel = TRUE, cpus = detectCores() - 1)
rol = nrow(Yused)-120
sfExport("path_result","Xused","Yused","h","rol")
sfLapply(seq_len(5*rol), once_function)
sfStop()



### read tables  PE 
file_name_PE_r4 <- paste("PE_r4.txt", sep = "_")
file_path_PE_r4 <- paste(path_result, file_name_PE_r4, sep = "")
PE_r4 = data.table::fread(file_path_PE_r4)


file_name_PE_rrr <- paste("PE_rrr3.txt", sep = "_")
file_path_PE_rrr <- paste(path_result, file_name_PE_rrr, sep = "")
PE_rrr = data.table::fread(file_path_PE_rrr)

file_name_PE_sss <- paste("PE_sss.txt", sep = "_")
file_path_PE_sss <- paste(path_result, file_name_PE_sss, sep = "")
PE_sss = data.table::fread(file_path_PE_sss)

file_name_PE_rcgl <- paste("PE_rcgl.txt", sep = "_")
file_path_PE_rcgl <- paste(path_result, file_name_PE_rcgl, sep = "")
PE_rcgl = data.table::fread(file_path_PE_rcgl)

file_name_PE_secure <- paste("PE_secure.txt", sep = "_")
file_path_PE_secure <- paste(path_result, file_name_PE_secure, sep = "")
PE_secure = data.table::fread(file_path_PE_secure)
PE_secure = PE_secure 



round(apply(PE_r4, 2, mean)  ,3)
round(apply(PE_rrr, 2, mean),3)
round(apply(PE_sss, 2, mean),3)
round(apply(PE_rcgl, 2, mean),3)
round(apply(PE_secure, 2, mean),3)


