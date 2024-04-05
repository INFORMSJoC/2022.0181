## data describe
# correspond to simulation data set with secure(X,B,E) + r4 outlier

############ simulatipn 

rm(list = ls())
setwd("~/scripts/Simu_Reproduce")
# load package
library(purrr) # for cross3
library(ncvreg) # for Mest computing
library(snowfall) # parallel programming
library(quadprog)
library(parallel)

##############    update to 22-4-17 (FARMA-1st) ################################
# initialization
path_result <- "Results-Sim22_5_1_context/"
if (!file_test("-d", path_result)) {
  dir.create(path_result)
}

# function of mask, swamp, detection
MSD = function(C_est,n,nout){
  mask = 0
  swamp = 0
  detection = 0
  for (i in 1:n) {
    if(i<=nout)
    {
      if( norm(C_est[i,],"2") == 0 )
      {mask = mask + 1}
    }else{
      if( norm(C_est[i,],"2") != 0 )
      {swamp = swamp + 1}
    }
  }
  
  if(mask == 0 & swamp == 0) detection = 1
  
  Mmask = mask/nout #the fraction of undetected outliers
  Mswamp = swamp/(n-nout) #the fraction of good points labeled as outliers
  Mdetection = detection #with no masking and no swamping.
  
  return(list(mask = Mmask, swamp = Mswamp, detection = Mdetection))
}





### step1 simulate 100 data set
runtimes = 100

cal_loss_once <- function( ii) 
{
  ## set 1
  set = "setting1"
  n = 400
  p = 500
  q = 200
  nrank = 3
  rho_X = 0.5
  rho_E = 0.5
  snr = 0.75
  nout = 5
  vout = NULL
  voutsd = 2
  nlev = 0
  vlev = NULL
  vlevsd = NULL
  sparsity = c(8,9,9) #sparsity = c(16,18,18)
  
  ### set 2
  #set = "setting2"
  #n = 400
  #p = 500
  #q = 200
  #nrank = 3
  #rho_X = 0.5
  #rho_E = 0.5
  #snr = 0.75
  #nout = 5
  #vout = NULL
  #voutsd = 2
  #nlev = 0
  #vlev = NULL
  #vlevsd = NULL
  #sparsity =  c(16,18,18)
  
  ### set 3
  # set = "setting3"
  # n = 400
  # p = 500
  # q = 200
  # nrank = 3
  # rho_X = 0.5
  # rho_E = 0.5
  # snr = 0.75
  # nout = 10
  # vout = NULL
  # voutsd = 2
  # nlev = 0
  # vlev = NULL
  # vlevsd = NULL
  # sparsity = c(8,9,9)
  
  ### set 4
  # set = "setting4"
  # n = 400
  # p = 500
  # q = 200
  # nrank = 3
  # rho_X = 0.5
  # rho_E = 0.5
  # snr = 0.75
  # nout = 10
  # vout = NULL
  # voutsd = 2
  # nlev = 0
  # vlev = NULL
  # vlevsd = NULL
  # sparsity =  c(16,18,18)
  
  ### set 5 
  # set = "setting5"
  # n = 400
  # p = 500
  # q = 200
  # nrank = 3
  # rho_X = 0.5
  # rho_E = 0.5
  # snr = 0.75
  # nout = 10
  # vout = NULL
  # voutsd = 4
  # nlev = 0
  # vlev = NULL
  # vlevsd = NULL
  # sparsity =  c(8,9,9)
  
  # ### set 6 
  # set = "setting6"
  # n = 400
  # p = 500
  # q = 200
  # nrank = 3
  # rho_X = 0.5
  # rho_E = 0.5
  # snr = 0.75
  # nout = 5
  # vout = NULL
  # voutsd = 4
  # nlev = 0
  # vlev = NULL
  # vlevsd = NULL
  # sparsity =   c(8,9,9)
  
  #
  method_list <- c("ROP", "ROPr", "RRR","RCGL","SECURE","R4")
  
  # output
  B_methods <- array(0, dim = length(method_list))
  C_methods <- array(0, dim = length(method_list))
  BC_methods <- array(0, dim = length(method_list))
  rank_methods = array(0, dim = length(method_list))
  time_methods = array(0, dim = length(method_list))
  mask_methods = array(0, dim = length(method_list))
  swamp_methods = array(0, dim = length(method_list))
  detection_methods = array(0, dim = length(method_list))
  
  # output files
  file_name_B <- paste(set, "B.txt", sep = "_")
  file_path_B <- paste(path_result, file_name_B, sep = "")
  
  file_name_C <- paste(set, "C.txt", sep = "_")
  file_path_C <- paste(path_result, file_name_C, sep = "")
  
  file_name_BC <- paste(set, "BC.txt", sep = "_")
  file_path_BC <- paste(path_result, file_name_BC, sep = "")
  
  file_name_rank <- paste(set, "rank.txt", sep = "_")
  file_path_rank <- paste(path_result, file_name_rank, sep = "")
  
  file_name_time <- paste(set, "time.txt", sep = "_")
  file_path_time <- paste(path_result, file_name_time, sep = "")
  
  file_name_mask <- paste(set, "mask.txt", sep = "_")
  file_path_mask <- paste(path_result, file_name_mask, sep = "")
  
  file_name_swamp <- paste(set, "swamp.txt", sep = "_")
  file_path_swamp <- paste(path_result, file_name_swamp, sep = "")
  
  file_name_detection <- paste(set, "detection.txt", sep = "_")
  file_path_detection <- paste(path_result, file_name_detection, sep = "")
  
  # debug
  file_name_debug <- paste(set, "debug.txt", sep = "_")
  file_path_debug <- paste(path_result, file_name_debug, sep = "")
  cat(ii, "\t", file = file_path_debug, append = TRUE)
  
  # library or source some useful packages and funcitons
  source("RS1.R")
  library("rrpack")
  source("setup_secure.R")
  source("r4.1.R")  # r4 something wrong with its original code
  source("RCGL.R")
  source("secure.R")
  
  # import data
  seed = 2022+ii
  simdata = setsecure(n = n,
                      p = p,
                      q = q,
                      s = sparsity,
                      nrank = nrank,
                      rho_X = rho_X,
                      rho_E = rho_E,
                      snr = snr,
                      nout = nout,
                      vout = vout,
                      voutsd = voutsd,
                      nlev = nlev,
                      vlev = vlev,
                      vlevsd = vlevsd) 
  
  Y <- simdata$Y
  X <- simdata$X
  B = simdata$B
  C = simdata$C
  
  trueEST = rbind(B,C)
  r = qr(trueEST)$rank
  Y.mean = X %*% B + C
  

  
  # ROP
  start1 = Sys.time()
  re1 = RS1(Y,X,r = (r+4),outlier = TRUE)
  stop1 = Sys.time()
  timecost1 = stop1 - start1
  rank.est1 = re1 $rank
  B.est1 = re1 $B
  C.est1 = re1 $C
  BC.est1 = re1 $coef
  temp1 = MSD(C.est1,n,nout)
  
  BC_methods[1] = norm(Y.mean - X %*% B.est1 - C.est1,"F")/norm(Y.mean,"F")
  B_methods[1] = norm(B - B.est1,"F")/norm(B,"F")
  C_methods[1] = norm(C - C.est1,"F")/norm(C,"F")
  rank_methods[1] = rank.est1
  time_methods[1] = timecost1
  mask_methods[1] = temp1$mask
  swamp_methods[1] = temp1$swamp
  detection_methods[1] = temp1$detection
  
  # ROPr  
  start2 = Sys.time()
  Y1 = Y - C.est1
  re2 = RS1(Y1,X,rank.est1,outlier = FALSE)
  stop2 = Sys.time()
  timecost2 = stop2 - start2 + timecost1
  B.est2 = re2$B
  rank.est2 = re2$rank
  
  BC_methods[2] = norm(Y.mean - X %*% B.est2 - C.est1,"F")/norm(Y.mean,"F")
  B_methods[2] = norm(B - B.est2,"F")/norm(B,"F")
  C_methods[2] = norm(C - C.est1,"F")/norm(C,"F")
  rank_methods[2] = rank.est2
  time_methods[2] = timecost2
  mask_methods[2] = temp1$mask
  swamp_methods[2] = temp1$swamp
  detection_methods[2] = temp1$detection
  
  # RRR
  start3 = Sys.time()
  re3 <- cv.rrr(Y, X, nfold = 10, maxrank = (r+4))
  stop3 = Sys.time()
  timecost3 = stop3 - start3
  B.est3 <- coef(re3)
  rank.est3 = re3$rank
  C.est3 = matrix(rep(0,n*q),n)
  
  BC_methods[3] = norm(Y.mean - X %*% B.est3 - C.est3,"F")/norm(Y.mean,"F")
  B_methods[3] = norm(B - B.est3,"F")/norm(B,"F")
  C_methods[3] = norm(C - C.est3,"F")/norm(C,"F")
  rank_methods[3] = rank.est3
  time_methods[3] = timecost3
  mask_methods[3] = 1
  swamp_methods[3] = 0
  detection_methods[3] = 0
  
  # RCGL
  start4 = Sys.time()
  re4 = RCGL(Y,X,rank=(r+4))
  stop4 = Sys.time()
  timecost4 = stop4 - start4
  rank.est4 = re4$r
  B.est4 = re4$coef
  C.est4 = matrix(rep(0,n*q),n)
  
  BC_methods[4] = norm(Y.mean - X %*% B.est4 - C.est4,"F")/norm(Y.mean,"F")
  B_methods[4] = norm(B - B.est4,"F")/norm(B,"F")
  C_methods[4] = norm(C - C.est4,"F")/norm(C,"F")
  rank_methods[4] = rank.est4
  time_methods[4] = timecost4
  mask_methods[4] = 1
  swamp_methods[4] = 0
  detection_methods[4] = 0
  
  # Secure
  start5 = Sys.time()
  re5 <- secure1(Y, X, rank = (r+4))
  stop5 = Sys.time()
  timecost5 = stop5 - start5
  B.est5 <- re5$coef
  rank.est5 = re5$rank
  C.est5 = matrix(rep(0,n*q),n)
  
  BC_methods[5] = norm(Y.mean - X %*% B.est5 - C.est5,"F")/norm(Y.mean,"F")
  B_methods[5] = norm(B - B.est5,"F")/norm(B,"F")
  C_methods[5] = norm(C - C.est5,"F")/norm(C,"F")
  rank_methods[5] = rank.est5
  time_methods[5] = timecost5
  mask_methods[5] = 1
  swamp_methods[5] = 0
  detection_methods[5] = 0
  
  # R4
  start6 = Sys.time()
  re6 = r4.1(Y,
             X,
             maxrank = (r+4),
             method = "rowl0",
             Gamma = NULL,
             ic.type ="PIC",
             modstr = list(),
             control = list())
  stop6 = Sys.time()
  timecost6 = stop6 - start6
  rank.est6 = re6$rank
  B.est6 =  coef(re6)
  C.est6 = re6$s
  temp6 = MSD(C.est6,n,nout)
  
  BC_methods[6] = norm(Y.mean - X %*% B.est6 - C.est6,"F")/norm(Y.mean,"F")
  B_methods[6] = norm(B - B.est6,"F")/norm(B,"F")
  C_methods[6] = norm(C - C.est6,"F")/norm(C,"F")
  rank_methods[6] = rank.est6
  time_methods[6] = timecost6
  mask_methods[6] = temp6$mask
  swamp_methods[6] = temp6$swamp
  detection_methods[6] = temp6$detection
  
  # output results
  Sys.sleep(runif(1, min=0, max=1))
  cat(BC_methods, "\n", file = file_path_BC, append = TRUE)
  cat(B_methods, "\n", file = file_path_B, append = TRUE)
  cat(C_methods, "\n", file = file_path_C, append = TRUE)
  cat(rank_methods, "\n", file = file_path_rank, append = TRUE)
  cat(time_methods, "\n", file = file_path_time, append = TRUE)
  cat(mask_methods, "\n", file = file_path_mask, append = TRUE)
  cat(swamp_methods, "\n", file = file_path_swamp, append = TRUE)
  cat(detection_methods, "\n", file = file_path_detection, append = TRUE)
  
  # debug
  cat("done\t", ii, "\n", file = file_path_debug, append = TRUE)
  # end debug
}

# snowfall
sfInit(parallel = TRUE, cpus = detectCores() - 1)
sfExport("path_result", "MSD")
sfLapply(1:100, cal_loss_once)
sfStop()


# draw table
set = "setting1"  #"setting2"
method_list <- c("ROP", "ROPr", "RRR","RCGL","SECURE","R4")
file_name_BC <- paste(set, "BC.txt", sep = "_")
file_path_BC <- paste(path_result, file_name_BC, sep = "")
rt1 = read.table(file_path_BC)
#index = which(rt1[,2]<rt1[,6])
index = 1:100
errorBC = as.numeric(round(apply(rt1[index,],2,mean), 3))
sdBC = as.numeric(round(apply(rt1[index,],2,sd), 3))

file_name_B <- paste(set, "B.txt", sep = "_")
file_path_B <- paste(path_result, file_name_B, sep = "")
rt = read.table(file_path_B)
errorB = as.numeric(round(apply(rt[index,],2,mean), 3))
sdB = as.numeric(round(apply(rt[index,],2,sd), 3))

file_name_C <- paste(set, "C.txt", sep = "_")
file_path_C <- paste(path_result, file_name_C, sep = "")
rt = read.table(file_path_C)
errorC = as.numeric(round(apply(rt[index,],2,mean), 3))
sdC = as.numeric(round(apply(rt[index,],2,sd), 3))

file_name_rank <- paste(set, "rank.txt", sep = "_")
file_path_rank <- paste(path_result, file_name_rank, sep = "")
rt = read.table(file_path_rank)
errorRank = as.numeric(round(apply(rt[index,],2,mean), 3))

file_name_time <- paste(set, "time.txt", sep = "_")
file_path_time <- paste(path_result, file_name_time, sep = "")
rt = read.table(file_path_time)
errorTime = as.numeric(round(apply(rt[index,],2,mean), 3))

file_name_mask <- paste(set, "mask.txt", sep = "_")
file_path_mask <- paste(path_result, file_name_mask, sep = "")
rt = read.table(file_path_mask)
errorMask = as.numeric(round(apply(rt[index,],2,mean), 3))


file_name_swamp <- paste(set, "swamp.txt", sep = "_")
file_path_swamp <- paste(path_result, file_name_swamp, sep = "")
rt = read.table(file_path_swamp)
errorSwamp = as.numeric(round(apply(rt[index,],2,mean), 3))

file_name_detection <- paste(set, "detection.txt", sep = "_")
file_path_detection <- paste(path_result, file_name_detection, sep = "")
rt = read.table(file_path_detection)
errorDetection = as.numeric(round(apply(rt[index,],2,sum)/length(index), 3))

RE = cbind(errorB,sdB,errorC,sdC, errorRank, errorMask, errorSwamp, errorDetection, errorTime)
rownames(RE)= method_list
print(RE)
